df <- neutrofil_count_Hildur
df <- df[5:nrow(df),4:8]
colnames(df) <- c("Four","Twelwe","Twentyfour","Fourtyeights", "Seventytwo")
# Organize data
n <- 22 #number of experiments
Times <- c(4,12,24,48,72)
nobs <- c(3,3,3,3,3,2,1,2,3,4,4,6,4,4,4,4,4,4,6,4,4,4)
start <- c(1,1+cumsum(nobs))
tt <- 0
nf <- 0

experiment <- data.frame(tt,nf)
for(i in 1:n){
  for(t in 2:5){
    for(j in start[i]:(start[i+1]-1)){
      tmp <- data.frame(Times[t],df[j,t])
      colnames(tmp) <- c("tt", "nf")
      experiment <- rbind(experiment, tmp)
    }
  }
}
experiment <- experiment[2:313,]
sstart <- c(1,1+cumsum(4*nobs))


pvals <- c()
ResTimes <- c()
t <- c(10:74)
plot(experiment$tt,-log(experiment$nf), ylim = c(-4,4), xlim=c(10,74))
for(i in 1:n){
    res <- lm(-log(nf) ~ tt, experiment[sstart[i]:(sstart[i+1]-1),])
    lines(t,res$coefficients[1]+res$coefficients[2]*t)
    fval <- summary(res)$fstatistic
    pvals <- c(pvals,pf(fval[1],fval[2],fval[3],lower.tail=F))
    ResTimes <- c(ResTimes, log(2)/res$coefficients[2])
}
pvals
ResTimes


# Individual experiments

i = 22
ex <- experiment[sstart[i]:(sstart[i+1]-1),]
res <- lm(-log(nf) ~ tt, ex)
plot(ex$tt, -log(ex$nf))
lines(t,res$coefficients[1]+res$coefficients[2]*t)
fval <- summary(res)$fstatistic
pf(fval[1],fval[2],fval[3],lower.tail=F)
log(2)/res$coefficients[2]


plot(ex$tt,ex$nf)
lines(t,exp(-res$coefficients[1]-res$coefficients[2]*t))


# Multiple regression analysis
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
regPlot3 <- function(dataList, yLims = NULL, xJitter = 1.2) {
  dataList <- lapply(seq_along(dataList), function(index) {
    dataList[[index]]$Data %<>% dplyr::mutate(., Experiment = str_c("Experiment", index, sep = " "))
    dataList[[index]]$Prediction %<>% dplyr::mutate(., Experiment = str_c("Experiment", index, sep = " "))
    dataList[[index]]$Summary %<>% dplyr::mutate(., Experiment = str_c("Experiment", index, sep = " "))
    
    return(dataList[[index]])
  })
  
  dataDf <- lapply(seq_along(dataList), function(index) dataList[[index]]$Data) %>% do.call("rbind", .)
  dataDf$Experiment %<>% factor(., levels = unique(.))
  
  predDf <- lapply(seq_along(dataList), function(index) dataList[[index]]$Prediction) %>% do.call("rbind", .)
  predDf$Experiment %<>% factor(., levels = unique(.))
  
  sumDf <- lapply(seq_along(dataList), function(index) dataList[[index]]$Summary) %>% do.call("rbind", .)
  sumDf$Experiment %<>% factor(., levels = unique(.))
  
  if(is.null(yLims) == TRUE) {
    yLims <- c(roundToInt(min(dataList$Data$logN), 0.5), roundToInt(max(dataList$Data$logN), 0.5))
  }
  
  ggplot2::ggplot() +
    geom_point(data = dataDf, mapping = aes(x = Time, y = -log(Neutrofill), color = Experiment), 
               position = position_jitter(width = xJitter)) +
    geom_line(data = predDf, mapping = aes(x = Time, y = pNeutrofill, color = Experiment), size = 0.5) +
    geom_ribbon(data = sumDf, mapping = aes(x = Time, y = Mean, ymax = Mean + MSE, ymin = Mean - MSE, fill = Experiment),
                size = 0.5, alpha = 0.3) +
    scale_y_continuous(limits = yLims, breaks = c(seq(yLims[1], yLims[2], by = 1)), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(12, max(sumDf$Time), by = 12)) +
    labs(y = expression(paste(-log[e], " neutrofill count"))) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
          axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
          legend.text = element_text(size = 10), legend.title = element_blank(),
          legend.direction = "vertical",
          legend.key.size = unit(0.6, "lines"), legend.position = "right") +
    guides(color = guide_legend(ncol = 1))
}

#Plot bargraph
plotSlope2 <- function(dataList) {
  dataDf <- lapply(dataList, function(resList) {
    sumDf <- summary(resList$Model)$coefficients
    
    dataDf <- data.frame(Slope = sumDf[2, 1],
                         StdError = sumDf[2, 2]
                         )
    }) %>% do.call("rbind", .) %>%
    dplyr::mutate(., Experiment = str_c("Experiment", seq_along(dataList), sep = " "))
  dataDf$Experiment %<>% factor(., levels = unique(.))
  
  ggplot(data = dataDf, mapping = aes(x = Experiment, y = Slope, fill = Experiment)) +
    geom_bar(stat = "identity", alpha = 0.5, width = 0.5) +
    geom_errorbar(data = dataDf, mapping = aes(y = Slope, x = Experiment, ymin = Slope - StdError, ymax = Slope + StdError),
                  size = 0.5, width = 0.25) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(dataDf$Slope)*1.4)) +
    labs(y = "Slope") +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
          axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
          legend.text = element_text(size = 12), legend.title = element_blank(),
          legend.key.size = unit(0.6, "lines"), legend.position = "top")
}

#Calculate mean of all slopes and plot as a bar graph
plotSlope3

res2 <- lapply(1:n, function(index) {
  experiment[sstart[index]:(sstart[index+1]-1),] %>%
    aprilReg2(.)
})

regPlots <- lapply(seq_along(res2), function(index) {
  print(index) 
  regPlot2(res2[[index]], yLims = c(-4, 4), expIndex = index)
})

pdf(file.path(getwd(), "hall.pdf"), width = 20, height = 20)
gridExtra::grid.arrange(grobs = regPlots)
dev.off()

pdf(file.path(getwd(), "hall2.pdf"), width = 6, height = 5)
regPlot3(res2,  yLims = c(-4, 4))
dev.off()

pdf(file.path(getwd(), "hall-Slope.pdf"), width = 10, height = 7)
plotSlope2(res2)
dev.off()

#Experiments to keep by Hildurs call: 9, 11,13,14,15,17,18,19,21,22
res2 %<>% .[c(9, 11, 13, 14, 15, 17, 18, 19, 21, 22)]

regPlots <- lapply(seq_along(res2), function(index) {
  print(index) 
  regPlot2(res2[[index]], yLims = c(-4, 4), expIndex = index)
})

pdf(file.path(getwd(), "newHall2.pdf"), width = 20, height = 10)
gridExtra::grid.arrange(grobs = regPlots, ncol = 5, nrow = 2)
dev.off()