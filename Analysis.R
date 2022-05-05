setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/April/Regression/Regression")

#Begin by reading the file "April_nov_2021_read.xlsx" into R
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

# Multiple regression analysis
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


df <- read.xlsx("April_nov_2021_read.xlsx")


#Wild type - loading data
WT_SH <- data.frame(df$Time...2, df$SHAM_WT)
WT_SH <- WT_SH[5:22,]
colnames(WT_SH) <- c("Time", "SHAM_WT")

#Alt
WT_SH <- data.frame(df[[2]], df$SHAM_WT)
WT_SH <- WT_SH[5:22,]
colnames(WT_SH) <- c("Time", "SHAM_WT")

#Henrik
WT_VNS <- data.frame(df$Time...4, df$VNS_WT)
WT_VNS <- WT_VNS[4:21,]
colnames(WT_VNS) <- c("Time", "VNS_WT")

#Alt
WT_VNS <- data.frame(df[[4]], df$VNS_WT)
WT_VNS <- WT_VNS[4:21,]
colnames(WT_VNS) <- c("Time", "VNS_WT")
 

#Alpha deficient - loading data
Alpha_SH <- data.frame(df$Time...6, df$SHAM_alpha)
Alpha_SH <- Alpha_SH[3:17,]
colnames(Alpha_SH) <- c("Time", "SHAM_alpha")

#Alt
Alpha_SH <- data.frame(df[[6]], df$SHAM_alpha)
Alpha_SH <- Alpha_SH[3:17,]
colnames(Alpha_SH) <- c("Time", "SHAM_alpha")


#Henrik
Alpha_VNS <- data.frame(df$Time...8, df$VNS_alpha)
Alpha_VNS <- Alpha_VNS[3:17,]
colnames(Alpha_VNS) <- c("Time", "VNS_alpha")

#Alt
Alpha_VNS <- data.frame(df[[8]], df$VNS_alpha)
Alpha_VNS <- Alpha_VNS[3:17,]
colnames(Alpha_VNS) <- c("Time", "VNS_alpha")


#Alox 15 deficient - loading data
Alox_SH <- data.frame(df$Time...8, df$SHAM_alox)
Alox_SH <- Alox_SH[3:16,]
colnames(Alox_SH) <- c("Time", "SHAM_alox")

#Alt
Alox_SH <- data.frame(df[[8]], df$SHAM_alox)
Alox_SH <- Alox_SH[3:16,]
colnames(Alox_SH) <- c("Time", "SHAM_alox")


#Henrik
Alox_VNS <- data.frame(df$Time...10, df$VNS_alox)
Alox_VNS <- Alox_VNS[3:16,]
colnames(Alox_VNS) <- c("Time", "VNS_alox")

#Alt
Alox_VNS <- data.frame(df[[10]], df$VNS_alox)
Alox_VNS <- Alox_VNS[3:16,]
colnames(Alox_VNS) <- c("Time", "VNS_alox")



# Multiple regression analysis

# Wild type
SHAM <- array("SHAM", dim(WT_SH)[1])
VNS <- array("VNS", dim(WT_VNS)[1])
WT <- data.frame(c(WT_SH$Time, WT_VNS$Time), c(WT_SH$SHAM_WT, WT_VNS$VNS_WT), c(SHAM, VNS))
colnames(WT) <- c("Time", "Neutrofill", "Type")
res_WT <- lm(-log(Neutrofill) ~ Time*Type, data = WT)
summary(res_WT)

res_WT <- lm(-log(Neutrofill) ~ Time + Time:Type, data = WT)
summary(res_WT)


# The log-linear plots
plot(WT_SH$Time, -log(WT_SH$SHAM_WT), col="red", ylim =c(-4,4), xlim = c(7,53),ylab = "neg log neutrofill", xlab="time (h)", main = "Wild type")
legend(45, -2,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)
t <- c(10:50)
lines(t,res_WT$coefficients[1]+res_WT$coefficients[2]*t, col="red")
points(WT_VNS$Time, -log(WT_VNS$VNS_WT), col="blue")
lines(t,res_WT$coefficients[1]+(res_WT$coefficients[2]+res_WT$coefficients[3])*t, col="blue")


#Get mean and sem per time point
sumWT <- WT %>% dplyr::group_by(Type, Time) %>%
  dplyr::summarise(Mean = mean(-log(Neutrofill)), na.rm = TRUE,
                   SD = sd(-log(Neutrofill)), na.rm = TRUE,
                   Count = sum(complete.cases(Neutrofill)))
sumWT %<>% dplyr::ungroup()
sumWT %<>% dplyr::mutate(SEM = .$SD/sqrt(.$Count))


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
sumWT %<>% dplyr::mutate(., Mean = unique(predWT$pNeutrofill))
sumWT %<>% dplyr::ungroup()

#ggplot
WT$Type %<>% factor(., levels = unique(.))
WT %<>% .[.$Time != 4,]

pdf(file.path(getwd(), "ggplot_WT.pdf"), width = 5, height = 5)
ggplot() +
  geom_point(data = WT, mapping = aes(x = Time, y = -log(Neutrofill), color = Type), 
             position = position_jitter(width = 0.6)) +
  geom_line(data = predWT, mapping = aes(x = Time, y = pNeutrofill, color = Type), size = 0.5) +
  geom_errorbar(data = sumWT, mapping = aes(x = Time, y = Mean, ymax = Mean + MSE, ymin = Mean - MSE, color = Type),
                width = 1, size = 0.5) +
  scale_y_continuous(limits = c(-4, 4)) +
  labs(y = expression(paste(-log[e], " neutrofill count"))) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
        axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.key.size = unit(0.6, "lines"), legend.position = "top")
dev.off()

pdf(file.path(getwd(), "ggplot_WT-area.pdf"), width = 5, height = 5)
ggplot() +
  geom_point(data = WT, mapping = aes(x = Time, y = -log(Neutrofill), color = Type), 
           position = position_jitter(width = 0.6)) +
  geom_line(data = predWT, mapping = aes(x = Time, y = pNeutrofill, color = Type), size = 0.5) +
  geom_ribbon(data = sumWT, mapping = aes(x = Time, y = Mean, ymax = Mean + MSE, ymin = Mean - MSE, fill = Type),
              size = 0.5, alpha = 0.3) +
  scale_y_continuous(limits = c(-4, 4)) +
  labs(y = expression(paste(-log[e], " neutrofill count"))) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11),
        axis.line = element_line(size = 0.6), axis.ticks = element_line(size = 0.6),
        legend.text = element_text(size = 10), legend.title = element_blank(),
        legend.key.size = unit(0.6, "lines"), legend.position = "top")
dev.off()

# Computing the resolution time
T_WT_SHAM <- log(2)/res_WT$coefficients[2] 
T_WT_VNS <- log(2)/(res_WT$coefficients[2] + res_WT$coefficients[3])
T_WT_SHAM
T_WT_VNS

# The exponential decay plots
plot(df$Time...2, df$SHAM_WT, col="red", ylim = c(0,25), xlim = c(0,55),  xlab = "Time (h)", ylab = "Neutrofill", main = "Wild type")
lines(t,exp(-(res_WT$coefficients[1]+res_WT$coefficients[2]*t)), col="red")
points(df$Time...4, df$VNS_WT, col="blue")
lines(t,exp(-(res_WT$coefficients[1]+(res_WT$coefficients[2]+res_WT$coefficients[3])*t)), col="blue")
legend(45, 25,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)

# Alpha deficient
SHAM <- array("SHAM", dim(Alpha_SH)[1])
VNS <- array("VNS", dim(Alpha_VNS)[1])
Alpha <- data.frame(c(Alpha_SH$Time, Alpha_VNS$Time), c(Alpha_SH$SHAM_alpha, Alpha_VNS$VNS_alpha), c(SHAM, VNS))
colnames(Alpha) <- c("Time", "Neutrofill", "Type")
res_Alpha <- lm(-log(Neutrofill) ~ Time*Type, data = Alpha)
summary(res_Alpha)

res_Alpha <- lm(-log(Neutrofill) ~ Time + Time:Type, data = Alpha)
summary(res_Alpha)

# The log-linear plots
plot(Alpha_SH$Time, -log(Alpha_SH$SHAM_alpha), col="red", ylim =c(-4,4), xlim = c(7,53), ylab = "neg log neutrofill", xlab="Time (h)", main = "Alpha7nAChR-deficient")
legend(45, -2,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)
t <- c(10:50)
lines(t,res_Alpha$coefficients[1]+res_Alpha$coefficients[2]*t, col="red")
points(Alpha_VNS$Time, -log(Alpha_VNS$VNS_alpha), col="blue")
lines(t,res_Alpha$coefficients[1]+(res_Alpha$coefficients[2]+res_Alpha$coefficients[3])*t, col="blue")

# Computing the resolution time
T_Alpha_SHAM <- log(2)/res_Alpha$coefficients[2] 
T_Alpha_VNS <- log(2)/(res_Alpha$coefficients[2] + res_Alpha$coefficients[3])
T_Alpha_SHAM
T_Alpha_VNS

# The exponential decay plots
plot(df$Time...6, df$SHAM_alpha, col="red", xlim = c(0,55), ylim = c(0,25), xlab = "Time (h)", ylab = "Neutrofill", main = "Alpha7nAChR-deficient")
lines(t,exp(-(res_Alpha$coefficients[1]+res_Alpha$coefficients[2]*t)), col="red")
points(df$Time...8, df$VNS_alpha, col="blue")
lines(t,exp(-(res_Alpha$coefficients[1]+(res_Alpha$coefficients[2]+res_Alpha$coefficients[3])*t)), col="blue")
legend(45, 25,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)


# Alox15 deficient
SHAM <- array("SHAM", dim(Alox_SH)[1])
VNS <- array("VNS", dim(Alox_VNS)[1])
Alox <- data.frame(c(Alox_SH$Time, Alox_VNS$Time), c(Alox_SH$SHAM_alox, Alox_VNS$VNS_alox), c(SHAM, VNS))
colnames(Alox) <- c("Time", "Neutrofill", "Type")
res_Alox <- lm(-log(Neutrofill) ~ Time*Type, data = Alox)
summary(res_Alox)

res_Alox <- lm(-log(Neutrofill) ~ Time+ Time:Type, data = Alox)
summary(res_Alox)

# The log-linear plots
plot(Alox_SH$Time, -log(Alox_SH$SHAM_alox), col="red", ylim =c(-4,4), xlim = c(7,53), ylab = "neg log neutrofill", xlab="time (h)", main = "Alox15-deficient")
legend(45, -2,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)
t <- c(10:50)
lines(t,res_Alox$coefficients[1]+res_Alox$coefficients[2]*t, col="red")
points(Alox_VNS$Time, -log(Alox_VNS$VNS_alox), col="blue")
lines(t,res_Alox$coefficients[1]+(res_Alox$coefficients[2]+res_Alox$coefficients[3])*t, col="blue")

# Computing the resolution time
T_Alox_SHAM <- log(2)/res_Alox$coefficients[2] 
T_Alox_VNS <- log(2)/(res_Alox$coefficients[2] + res_Alox$coefficients[3])
T_Alox_SHAM
T_Alox_VNS

# The log-linear plots
plot(df$Time...10, df$SHAM_alox, col="red", ylim = c(0,10), xlim = c(0,55),  xlab = "Time (h)", ylab = "Neutrofill", main = "Alox15-deficient")
lines(t,exp(-(res_Alox$coefficients[1]+res_Alox$coefficients[2]*t)), col="red")
points(df$Time...12, df$VNS_alox, col="blue")
lines(t,exp(-(res_Alox$coefficients[1]+(res_Alox$coefficients[2]+res_Alox$coefficients[3])*t)), col="blue")
legend(45, 10,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)

# All the slopes in the same graph

plot(WT_SH$Time, -log(WT_SH$SHAM_WT), col="red", ylim =c(-4,4), xlim = c(7,53),ylab = "neg log neutrofill", xlab="time (h)", main = "Wild type")
legend(45, -2,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)
t <- c(10:50)
lines(t,res_WT$coefficients[1]+res_WT$coefficients[2]*t, col="red")
points(WT_VNS$Time, -log(WT_VNS$VNS_WT), col="blue")
lines(t,res_WT$coefficients[1]+(res_WT$coefficients[2]+res_WT$coefficients[3])*t, col="blue")


points(Alpha_SH$Time, -log(Alpha_SH$SHAM_alpha), col="magenta")
#legend(45, -2,   # Coordinates (x also accepts keywords)
  #     legend = c("SHAM", "VNS"), # Vector with the name of each group
   #    fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
    #   col = par("red", "blue"), # Color of lines or symbols  
#)
t <- c(10:50)
lines(t,res_Alpha$coefficients[1]+res_Alpha$coefficients[2]*t, col="magenta")
points(Alpha_VNS$Time, -log(Alpha_VNS$VNS_alpha), col="green")
lines(t,res_Alpha$coefficients[1]+(res_Alpha$coefficients[2]+res_Alpha$coefficients[3])*t, col="green")

points(Alox_SH$Time, -log(Alox_SH$SHAM_alox), col="purple")
lines(t,res_Alox$coefficients[1]+res_Alox$coefficients[2]*t, col="purple")
points(Alox_VNS$Time, -log(Alox_VNS$VNS_alox), col="yellow")
lines(t,res_Alox$coefficients[1]+(res_Alox$coefficients[2]+res_Alox$coefficients[3])*t, col="yellow")






#WT
df <- read.xlsx("April_nov_2021_read.xlsx")

WT_SH <- data.frame(df[[2]], df$SHAM_WT)
WT_SH <- WT_SH[5:22,]
colnames(WT_SH) <- c("Time", "SHAM_WT")

WT_VNS <- data.frame(df[[4]], df$VNS_WT)
WT_VNS <- WT_VNS[4:21,]
colnames(WT_VNS) <- c("Time", "VNS_WT")

wt <- aprilReg(WT_SH, WT_VNS)

dev.off()

pdf(file.path(getwd(), "ggplot_WT.pdf"), width = 5, height = 5)
regPlot(wt)
dev.off()

pdf(file.path(getwd(), "slopes_WT.pdf"), width = 3, height = 5)
plotSlope(wt)
dev.off()




#Alox 15 deficient - loading data
Alox_SH <- data.frame(df[[10]], df$SHAM_alox)
Alox_SH <- Alox_SH[3:16,]
colnames(Alox_SH) <- c("Time", "SHAM_alox")

Alox_VNS <- data.frame(df[[12]], df$VNS_alox)
Alox_VNS <- Alox_VNS[3:16,]
colnames(Alox_VNS) <- c("Time", "VNS_alox")


alox <- aprilReg(Alox_SH, Alox_VNS)

pdf(file.path(getwd(), "ggplot_Alox.pdf"), width = 3, height = 5)
regPlot(alox, yLims = c(-3, 3), xJitter = 1.2)
dev.off()

pdf(file.path(getwd(), "slopes_Alox.pdf"), width = 3, height = 5)
plotSlope(alox)
dev.off()


#Alpha7
Alpha_SH <- data.frame(df[[6]], df$SHAM_alpha)
Alpha_SH <- Alpha_SH[3:17,]
colnames(Alpha_SH) <- c("Time", "SHAM_alpha")

Alpha_VNS <- data.frame(df[[8]], df$VNS_alpha)
Alpha_VNS <- Alpha_VNS[3:17,]
colnames(Alpha_VNS) <- c("Time", "VNS_alpha")


a7 <- aprilReg(Alpha_SH, Alpha_VNS)

dev.off()

pdf(file.path(getwd(), "ggplot_a7.pdf"), width = 5, height = 5)
regPlot(a7)
dev.off()

pdf(file.path(getwd(), "slopes_a7.pdf"), width = 3, height = 5)
plotSlope(a7)
dev.off()

gWT <- gridExtra::arrangeGrob(grobs = list(regPlot(wt), plotSlope(wt)), 
                              ncol = 2, nrow = 1, left = grid::textGrob("Wild Type", rot = 90))
gAlox <- gridExtra::arrangeGrob(grobs = list(regPlot(alox, yLims = c(-3, 3)), plotSlope(alox)), 
                                ncol = 2, nrow = 1, left = grid::textGrob("Alox15 KO", rot = 90))
gA7 <- gridExtra::arrangeGrob(grobs = list(regPlot(a7), plotSlope(a7)), 
                              ncol = 2, nrow = 1, left = grid::textGrob("Chrna7 KO", rot = 90))

pdf(file.path(getwd(), "all.pdf"), width = 6, height = 12)
gridExtra::grid.arrange(gWT, gAlox, gA7,  nrow = 3, ncol = 1)
dev.off()