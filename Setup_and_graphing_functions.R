#Setup and Graphing functions
library(openxlsx)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggpubr)

roundToInt <- function(x, m) {
  if(x > 0) {
    r <- ceiling(x + m - x %% m)
  } else {
    x %<>% abs()
    r <- ceiling(x + m - x %% m)
    r <- r*-1
  }
  
  return(r)
}

#Merge lists into data.frames where names of list components form a new column (with name from newColName)
pcrMerge <- function(object, newColName = "Experiment", renewIds = TRUE, saveAllCols = TRUE, oblColumns = NULL, newId = TRUE) {
  #Merges list of data.frames into one data.frame:
  #object - named list where names have meaning;
  #newColName - name for a new column that will be created from list names;
  #renewIds - if Id column is present in components - make a new one;
  #saveAllCols - preserve all columns from all data.frames;
  #oblColumns - character vector of all column names to preserve - the rest are removed.
  
  if(is.null(names(object)) == TRUE) {
    stop("Merging impossible - object lacks names!")
  } else {
    #Manage columns
    if(is.null(oblColumns) == TRUE) {
      oblColumns <- lapply(object, function(experiment) colnames(experiment))
      
      #Keep all columns or only the common ones
      if(saveAllCols == TRUE) {
        oblColumns %<>% unlist() %>% unique()
      } else {
        allColNames <- oblColumns %>% unlist() %>% unique()
        
        namesToRemove <- lapply(oblColumns, function(experiment) allColNames %>% .[!. %in% experiment]) %>% unlist()
        
        oblColumns %<>% unlist() %>% unique() %>% .[!. %in% namesToRemove]
      }
    } else {
      #Remove unwanted columns
      object %<>% lapply(., function(experiment) experiment %<>% .[, colnames(.) %in% oblColumns])
    }
    
    #Check if any data.frames lack columns
    namesToAdd <- lapply(object, function(experiment) oblColumns %>% .[!. %in% colnames(experiment)]) %>% unlist() %>% length()
    
    if(namesToAdd > 0) {
      object %<>% lapply(., function(experiment) {
        absentColumns <- sum(!oblColumns %in% colnames(experiment))
        
        if(absentColumns > 0) {
          for(i in 1:absentColumns) {
            experiment <- cbind(experiment, NA)
          }
          colnames(experiment) <- c(colnames(experiment)[!colnames(experiment) %in% "NA"], 
                                    oblColumns[!oblColumns %in% colnames(experiment)])
        }
        
        return(experiment)
      })
    }
    
    #Extract list component names
    experimentNames <- names(object)
    
    #Add component names to the data.frames
    object <- lapply(seq_along(object), function(index) {
      object[[index]] %<>% dplyr::mutate(., NewCol = names(object)[index]) 
    })
    object %<>% do.call("rbind", .)
    colnames(object)[colnames(object) == "NewCol"] <- newColName
    
    #Put new column at the front
    object %<>% .[,c(newColName, colnames(object)[c(1:ncol(object)-1)])]
    
    if("Id" %in% colnames(object)) {
      if(newId == TRUE) {
        object$Id <- c(1:nrow(object))
      }
      
      allCols <- colnames(object) %>% .[!. %in% "Id"]
      object %<>% .[,c("Id", allCols)]
    }
    
    return(object)
  }
}

#Multiple regression analysis
aprilReg <- function(cntrlVec, effVec) {
  # Wild type
  SHAM <- array("SHAM", dim(cntrlVec)[1])
  VNS <- array("VNS", dim(effVec)[1])
  WT <- data.frame(c(cntrlVec$Time, effVec$Time), c(cntrlVec[[2]], effVec[[2]]), c(SHAM, VNS))
  colnames(WT) <- c("Time", "Neutrofill", "Type")
  
  res_WT <- lm(-log(Neutrofill) ~ Time + Time:Type, data = WT)
  summary(res_WT)
  
  # The log-linear plots
  plot(cntrlVec$Time, -log(cntrlVec[[2]]), col="red", ylim =c(-4,4), xlim = c(7,53),
       ylab = "neg log neutrofill",  xlab="time (h)", main = "Wild type")
  legend(45, -2,   # Coordinates (x also accepts keywords)
         legend = c("SHAM", "VNS"), # Vector with the name of each group
         fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
         col = par("red", "blue"), # Color of lines or symbols  
  )
  t <- c(10:50)
  lines(t,res_WT$coefficients[1]+res_WT$coefficients[2]*t, col="red")
  points(effVec$Time, -log(effVec[[2]]), col="blue")
  lines(t,res_WT$coefficients[1]+(res_WT$coefficients[2]+res_WT$coefficients[3])*t, col="blue")
  
  #MSE
  predWT <- data.frame(pNeutrofill = predict(res_WT, WT), Time = WT$Time)
  predWT$Type <- WT$Type
  
  WT %<>% dplyr::mutate(logN = -log(.$Neutrofill))
  
  mseWT <- WT %>%
    dplyr::mutate(., pLogNeutrofill = predWT[predWT$Time == .$Time & predWT$Type == .$Type, "pNeutrofill"])
  mseWT %<>% dplyr::mutate(., Distance = .$logN - .$pLogNeutrofill)
  mseWT %<>% dplyr::mutate(., sqDistance = .$Distance ^ 2)
  
  sumWT <- mseWT %>% dplyr::group_by(Type, Time) %>%
    dplyr::summarise(MSE = mean(sqDistance, na.rm = TRUE))
  sumWT %<>% dplyr::ungroup()
  sumWT %<>% dplyr::mutate(., Mean = unique(predWT$pNeutrofill))
  
  WT$Type %<>% factor(., levels = unique(.))
  WT %<>% .[.$Time != 4,]
  
  dataList <- list(Data = WT,
                   Model = res_WT,
                   Prediction = predWT,
                   Summary = sumWT)
  
  return(dataList)
}
aprilReg2 <- function(dataDf) {
  colnames(dataDf) <- c("Time", "Neutrofill")
  
  dataDf %<>% .[!is.na(.$Neutrofill),]
  
  resModel <- lm(-log(Neutrofill) ~ Time, data = dataDf)
  
  #MSE
  predData <- data.frame(pNeutrofill = predict(resModel, dataDf), Time = dataDf$Time)
  
  dataDf %<>% dplyr::mutate(logN = -log(.$Neutrofill))
  
  mse <- dataDf %>%
    dplyr::mutate(., pLogNeutrofill = predData[predData$Time == .$Time, "pNeutrofill"])
  mse %<>% dplyr::mutate(., Distance = .$logN - .$pLogNeutrofill)
  mse %<>% dplyr::mutate(., sqDistance = .$Distance ^ 2)
  
  sumDf <- mse %>% dplyr::group_by(Time) %>%
    dplyr::summarise(MSE = mean(sqDistance, na.rm = TRUE))
  sumDf %<>% dplyr::ungroup()
  sumDf %<>% dplyr::mutate(., Mean = unique(predData$pNeutrofill))
  
  dataDf %<>% .[.$Time != 4,]
  
  dataList <- list(Data = dataDf,
                   Model = resModel,
                   Prediction = predData,
                   Summary = sumDf)
  
  return(dataList)
}

#Regression line graph
regPlot <- function(dataList, yLims = NULL, xJitter = 1.2) {
  if(is.null(yLims) == TRUE) {
    yLims <- c(roundToInt(min(dataList$Data$logN), 0.5), roundToInt(max(dataList$Data$logN), 0.5))
  }
  
  ggplot2::ggplot() +
    geom_point(data = dataList$Data, mapping = aes(x = Time, y = -log(Neutrofill), color = Type), 
               position = position_jitter(width = xJitter)) +
    geom_line(data = dataList$Prediction, mapping = aes(x = Time, y = pNeutrofill, color = Type), size = 0.5) +
    geom_ribbon(data = dataList$Summary, mapping = aes(x = Time, y = Mean, ymax = Mean + MSE, ymin = Mean - MSE, fill = Type),
                size = 0.5, alpha = 0.3) +
    scale_y_continuous(limits = yLims, breaks = c(seq(yLims[1], yLims[2], by = 1)), expand = c(0,0)) +
    scale_x_continuous(breaks = c(12, 24, 36, 48)) +
    labs(y = expression(paste(-log[e], " neutrofill count"))) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
          axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
          legend.text = element_text(size = 10), legend.title = element_blank(),
          legend.key.size = unit(0.6, "lines"), legend.position = "top")
}
regPlot2 <- function(dataList, yLims = NULL, xJitter = 1.2, expIndex = 1) {
  if(is.null(yLims) == TRUE) {
    yLims <- c(roundToInt(min(dataList$Data$logN), 0.5), roundToInt(max(dataList$Data$logN), 0.5))
  }
  
  ggplot2::ggplot() +
    geom_point(data = dataList$Data, mapping = aes(x = Time, y = -log(Neutrofill), color = "tomato"), 
               position = position_jitter(width = xJitter)) +
    geom_line(data = dataList$Prediction, mapping = aes(x = Time, y = pNeutrofill), size = 0.5, color = "tomato") +
    geom_ribbon(data = dataList$Summary, mapping = aes(x = Time, y = Mean, ymax = Mean + MSE, ymin = Mean - MSE, fill = "tomato"),
                size = 0.5, alpha = 0.3) +
    scale_y_continuous(limits = yLims, breaks = c(seq(yLims[1], yLims[2], by = 1)), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(12, max(dataList$Summary), by = 12)) +
    labs(title = str_c("Experiment", expIndex, sep = " "), 
         y = expression(paste(-log[e], " neutrofill count"))) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
          axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
          legend.text = element_text(size = 10), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15), 
          legend.key.size = unit(0.6, "lines"), legend.position = "none") 
}

#Slope bargraph
plotSlope <- function(dataList) {
  sumDf <- summary(dataList$Model)$coefficients
  
  dataDf <- data.frame(CondGroup = factor(unique(dataList$Data$Type), levels = unique(dataList$Data$Type)),
                       Slope = c(sumDf[2, 1], sumDf[2, 1] + sumDf[3, 1]),
                       StdError = sumDf[c(2:3), 2]
  )
  
  path <- data.frame(CondGroup = rep(dataDf$CondGroup, each = 2),
                     y = c(max(dataDf$Slope) + max(dataDf$Slope)*0.2,
                           max(dataDf$Slope) + max(dataDf$Slope)*0.25,
                           max(dataDf$Slope) + max(dataDf$Slope)*0.25,
                           max(dataDf$Slope) + max(dataDf$Slope)*0.2
                     ),
                     Group = 1
  )
  
  
  ggplot(data = dataDf, mapping = aes(x = CondGroup, y = Slope, fill = CondGroup)) +
    geom_bar(stat = "identity", alpha = 0.5, width = 0.5) +
    geom_errorbar(data = dataDf, mapping = aes(y = Slope, x = CondGroup, ymin = Slope - StdError, ymax = Slope + StdError),
                  size = 0.5, width = 0.25) +
    geom_path(data = path, mapping = aes(y = y, x = CondGroup, group = Group), size = 0.6) +
    geom_text(mapping = aes(x = 1.5, y = max(dataDf$Slope)*1.3), size = 3.5,
              label = str_c("p = ", format(summary(dataList$Model)$coefficients[3, "Pr(>|t|)"], digits = 2), sep = ""))  +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(dataDf$Slope)*1.4)) +
    labs(y = "Slope") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
          axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
          legend.text = element_text(size = 12), legend.title = element_blank(),
          legend.key.size = unit(0.6, "lines"), legend.position = "top")
}

#Plot curvies
plotCurvies <- function(pointsA, pointsB, linesTimeA, linesNeuA, linesTimeB, linesNeuB, nameA = "Sham", nameB = "VNS", 
                        xJitter = 1.25, xBreak = 2.5, yLimits = NULL, addMeans = FALSE, 
                        addLines = FALSE, lineList = NULL, lineColor = "gray75") {
  #tt - time, x axis
  #nf - neutrofil count, y axis
  
  #Make points dataframe
  pointsDf <- list(pointsA, pointsB)
  names(pointsDf) <- c(nameA, nameB)
  
  pointsDf %<>% pcrMerge(., "Treatment")
  colnames(pointsDf) <- c("Treatment", "Time", "NF")
  
  #Make lines dataframe
  linesDf <- list(A = data.frame(`Time` = linesTimeA,
                                 `NF` = linesNeuA
  ),
  B = data.frame(`Time` = linesTimeB,
                 `NF` = linesNeuB
  ))
  names(linesDf) <- c(nameA, nameB)
  
  linesDf %<>% pcrMerge(., "Treatment")
  
  dataList <- list(Points = pointsDf,
                   Lines = linesDf
  )
  
  if(addMeans == TRUE) {
    #Calculates Means and SEMS
    sumDf <- dataList$Points %>%
      dplyr::group_by_at(., c("Treatment", "Time")) %>%
      dplyr::summarize(Mean = mean(NF, na.rm = TRUE),
                       SD = sd(NF, na.rm = TRUE),
                       Var = var(NF, na.rm = TRUE),
                       Count = sum(complete.cases(NF))
      ) %>% dplyr::ungroup()
    
    sumDf %<>% dplyr::mutate(., SEM = .$SD/sqrt(.$Count))
    
    if(is.null(yLimits) == TRUE) {
      yLimits <- c(0, roundToInt(max(sumDf$Mean + sumDf$SEM), xBreak))
    }
    
    Plot <-  Plot <- ggplot2::ggplot() +
      geom_point(data = sumDf, mapping = aes(x = Time, y = Mean, color = Treatment), size = 2) +
      geom_errorbar(data = sumDf, mapping = aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM, color = Treatment), width = 2) +
      geom_line(data = dataList$Lines, mapping = aes(x = Time, y = NF, color = Treatment), size = 0.5) +
      scale_y_continuous(expand = c(0,0), limits = yLimits,
                         breaks = seq(0, yLimits[2], by = xBreak)) +
      scale_x_continuous(expand = c(0,0), breaks = c(0, 4, 12, 24, 36, 48)) +
      labs(y = "Neutrofill count", x = "Time (h)") +
      theme_classic() +
      theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
            axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
            legend.text = element_text(size = 10), legend.title = element_blank(),
            legend.key.size = unit(0.6, "lines"), legend.position = "top")
  } else {
    if(is.null(yLimits) == TRUE) {
      yLimits <- c(0, roundToInt(max(dataList$Points$NF), xBreak))
    }
    
    Plot <- ggplot2::ggplot() +
      geom_point(data = dataList$Points, mapping = aes(x = Time, y = NF, color = Treatment), 
                 position = position_jitter(width = xJitter)) +
      geom_line(data = dataList$Lines, mapping = aes(x = Time, y = NF, color = Treatment), size = 0.5) +
      scale_y_continuous(expand = c(0,0), limits = yLimits,
                         breaks = seq(0, yLimits[2], by = xBreak)) +
      scale_x_continuous(expand = c(0,0), breaks = c(0, 4, 12, 24, 36, 48)) +
      labs(y = "Neutrofill count", x = "Time (h)") +
      theme_classic() +
      theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
            axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
            legend.text = element_text(size = 10), legend.title = element_blank(),
            legend.key.size = unit(0.6, "lines"), legend.position = "top")
  }
  
  if(addLines == TRUE) {
    if(is.null(lineList) == TRUE) {
      stop("List of data to make lines is absent!")
    } else {
      lineDf <- lineList %>% do.call("rbind", .) %>% tibble::as_tibble()
      lineDf$Treatment <- names(lineList)
      
      #Draw max line for control
      lineDf %<>% rbind(.,      
                        tibble::tibble(
                          Val_T_max = 0,
                          T_max = unname(unlist(lineDf[lineDf$Treatment == nameA, "T_max"])),
                          `Val50%T_max` = 0,
                          `Time_below_50%_Val` = 0,
                          Resolution_time = 0,
                          Treatment = nameA
                        )
      )
      
      lineDf <- lineDf[lineDf$Treatment == nameA,]
      
      Plot <- Plot + 
        geom_line(data = lineDf[lineDf$Treatment == nameA,], mapping = aes(y = Val_T_max, x = T_max, group = Treatment),
                  color = lineColor, linetype = "dashed")
      
      #Draw arrow for V max:
      lineDf <- lineList[names(lineList) == nameA] %>%
        do.call("rbind", .) %>%
        tibble::as_tibble()
      lineDf$Treatment <- nameA
      
      lineDf %<>% rbind(.,      
                        tibble::tibble(
                          Val_T_max = unname(unlist(lineDf[lineDf$Treatment == nameA, "Val_T_max"])),
                          T_max = 0.4,
                          `Val50%T_max` = 0,
                          `Time_below_50%_Val` = 0,
                          Resolution_time = 0,
                          Treatment = nameA
                        )
      )
      
      Plot <- Plot + 
        geom_line(data = lineDf[lineDf$Treatment == nameA,], mapping = aes(y = Val_T_max, x = T_max, group = Treatment),
                  color = lineColor, size = 1.5, arrow = arrow(type = "closed", ends = "first", length = unit(0.1, "inches")))
      
      #Draw half-peak lines:
      #Horizontal:
      lineDf <- lineList[names(lineList) == nameA] %>%
        do.call("rbind", .) %>%
        tibble::as_tibble()
      lineDf$Treatment <- nameA
      
      lineDf %<>% rbind(.,      
                        tibble::tibble(
                          Val_T_max = 0,
                          T_max = 0,
                          `Val50%T_max` = unname(unlist(lineDf[lineDf$Treatment == nameA, "Val50%T_max"])),
                          `Time_below_50%_Val` = 0,
                          Resolution_time = 0,
                          Treatment = nameA
                        )
      )
      
      Plot <- Plot + 
        geom_line(data = lineDf[lineDf$Treatment == nameA,], mapping = aes(y = `Val50%T_max`, x = `Time_below_50%_Val`, group = Treatment),
                  color = lineColor, linetype = "dashed")
      
      #Vertical
      lineDf <- lineList %>%
        do.call("rbind", .) %>%
        tibble::as_tibble()
      lineDf$Treatment <- names(lineList)
      
      #Adjust to match SHAM
      lineDf[lineDf$Treatment == "VNS", "Val50%T_max"] <- lineDf[lineDf$Treatment == "Sham", "Val50%T_max"]
      
      lineDf %<>% rbind(.,      
                        tibble::tibble(
                          Val_T_max = c(0, 0),
                          T_max = c(0, 0),
                          `Val50%T_max` = c(0, 0),
                          `Time_below_50%_Val` = unname(unlist(lineDf[, "Time_below_50%_Val"])),
                          Resolution_time = c(0, 0),
                          Treatment = names(lineList)
                        )
      )
      
      Plot <- Plot + 
        geom_line(data = lineDf, mapping = aes(y = `Val50%T_max`, x = `Time_below_50%_Val`, group = Treatment),
                  color = lineColor, linetype = "dashed")
      
      #Arrows for resolution index
      lineDf <- lineList %>%
        do.call("rbind", .) %>%
        tibble::as_tibble()
      lineDf$Treatment <- names(lineList)
      
      lineDf %<>% rbind(.,      
                        tibble::tibble(
                          Val_T_max = 0,
                          T_max = 0,
                          `Val50%T_max` = 0,
                          `Time_below_50%_Val` = unname(unlist(lineDf[lineDf$Treatment == nameA, "T_max"])),
                          Resolution_time = 0,
                          Treatment = names(lineList)
                        )
      ) %>% 
        dplyr::mutate(y = rep(seq(from = max(lineDf[lineDf$Treatment == nameA, "Val50%T_max"])/5,
                                  to = max(lineDf[lineDf$Treatment == nameA, "Val50%T_max"])/2, 
                                  length.out = length(lineList)
        ), times = length(lineList)
        )
        )
      
      Plot <- Plot + 
        geom_line(data = lineDf, mapping = aes(y = y, x = `Time_below_50%_Val`, group = Treatment),
                  color = lineColor, arrow = arrow(type = "closed", ends = "both", length = unit(0.075, "inches")))
    }
  }
  
  return(Plot)
}
