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



#Regression plots for WT, ALOX15 KO and ALPHA7 KO

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



#ODE models
# Preliminaries
#Read data using inbuilt function
df <- April_nov_2021_read

SHAM_WT <- data.frame(df$Time...2,df$SHAM_WT) 
VNS_WT <- data.frame(df$Time...4, df$VNS_WT)
VNS_WT <- VNS_WT[1:21,]
SHAM_alpha <- data.frame(df$Time...6,df$SHAM_alpha) 
SHAM_alpha <- SHAM_alpha[1:17,]
VNS_alpha <- data.frame(df$Time...8,df$VNS_alpha) 
VNS_alpha <- VNS_alpha[1:17,]
SHAM_alox <- data.frame(df$Time...10,df$SHAM_alox) 
SHAM_alox <- SHAM_alox[1:16,]
VNS_alox <- data.frame(df$Time...12,df$VNS_alox) 
VNS_alox <- VNS_alox[1:16,]
colnames(SHAM_WT) <- c("tt", "nf")
colnames(VNS_WT) <- c("tt", "nf")
colnames(SHAM_alpha) <- c("tt", "nf")
colnames(VNS_alpha) <- c("tt", "nf")
colnames(SHAM_alox) <- c("tt", "nf")
colnames(VNS_alox) <- c("tt", "nf")

# define functions
de <- diffeqr::diffeq_setup()

# p = c(alpha, mu_max, gamma, K) parameters
mufun <- function(s, mu_max, K){
  return(mu_max * s/(s + K))
}

f <- function(u,p,t) {
  du1 = (mufun(u[2],p[2],p[4])-p[3])*u[1]
  du2 = -p[1]* mufun(u[2],p[2],p[4])*u[2]
  return(c(du1,du2))
}


tspan <- list(1,48)
saveat = c(4,12,24,48)


# Analysis 

DD <- SHAM_WT
plot(DD$tt, -log(DD$nf), ylim = c(-4,4))

#p <- c(0.699,2.55,0.0726,30)
#prob <- de$ODEProblem(f, u0, tspan, p)
#sol <- de$solve(prob, saveat=saveat)

loglike <- function(p, DD){
  prob <- de$ODEProblem(f, u0, tspan, p)
  sol <- de$solve(prob, saveat=saveat)
  mat <- sapply(sol$u,identity)
  udf <- as.data.frame(t(mat))
  logl <- c(0)
  for(i in 1:length(DD$tt)){
    if(DD$tt[i] == 4)
      ref <- -log(udf$V1[1])
    if(DD$tt[i] == 12)
      ref <- -log(udf$V1[2])
    if(DD$tt[i] == 24)
      ref <- -log(udf$V1[3])
    if(DD$tt[i] == 48)
      ref <- -log(udf$V1[4])
    if(DD$tt[i] == 72)
      ref <- -log(udf$V1[5])
    logl <- logl + (-log(DD$nf[i])-ref)^2
  }
  return(logl)
}

p <- c(0.699,2.55,0.0726,30)
u0 <- c(1,1)
res_SHAM_WT <- optim(p,loglike, DD = SHAM_WT)
p_SHAM_WT <- res_SHAM_WT$par
res_VNS_WT <- optim(p,loglike, DD = VNS_WT)
p_VNS_WT <- res_VNS_WT$par

u0 <- c(0.75,1)
res_SHAM_alpha <- optim(p,loglike, DD = SHAM_alpha)
p_SHAM_alpha <- res_SHAM_alpha$par
res_VNS_alpha <- optim(p,loglike, DD = VNS_alpha)
p_VNS_alpha <- res_VNS_alpha$par

u0 <- c(0.1,1)
res_SHAM_alox <- optim(p,loglike, DD = SHAM_alox)
p_SHAM_alox <- res_SHAM_alox$par
res_VNS_alox <- optim(p,loglike, DD = VNS_alox)
p_VNS_alox <- res_VNS_alox$par



# Plotting results

#ggplot
library(ggplot2)
plotCurvies <- function(pointsA, pointsB, linesTimeA, linesNeuA, linesTimeB, linesNeuB, 
                        nameA = "Sham", nameB = "VNS", xJitter = 1.25, xBreak = 2.5, yLimits = NULL, addMeans = FALSE) {
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
  
  return(Plot)
}

pdf("Curvies_WT.pdf", height = 2.75591, width = 2.75591)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()

# WT

tspan <- list(0,50)
saveat_fit = c(0:500)/10

p_fit <- p_SHAM_WT
u0 <- c(1,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_SHAM <- de$solve(prob, saveat=saveat_fit)
mat_SHAM <- sapply(sol_SHAM$u,identity)
udf_SHAM <- as.data.frame(t(mat_SHAM))

p_fit <- p_VNS_WT
u0 <- c(1,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_VNS <- de$solve(prob, saveat=saveat_fit)
mat_VNS <- sapply(sol_VNS$u,identity)
udf_VNS <- as.data.frame(t(mat_VNS))


plot(SHAM_WT$tt, SHAM_WT$nf, ylim = c(0,30), xlim = c(0,50), col ="red", ylab = "Neutrofill count", xlab = "Time (h)")
points(VNS_WT$tt, VNS_WT$nf, col = "blue")
lines(sol_SHAM$t, udf_SHAM$V1, col = "red")
lines(sol_VNS$t, udf_VNS$V1, col = "blue")
legend(40, 30,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)





plotWT <- plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))

pdf("Curvies_WT-mid.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()

pdf("Curvies_WT-big.pdf", height = 5, width = 5)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()


pdf("Curvies_WT-variants_mid2.pdf", height = 3.5, width = 10.5)
gridExtra::grid.arrange(
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()

pdf("Curvies_WT-variants_big2.pdf", height = 5, width = 15)
gridExtra::grid.arrange(
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()


# Alpha

tspan <- list(0,50)
saveat_fit = c(0:500)/10

p_fit <- p_SHAM_alpha
u0 <- c(0.75,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_SHAM <- de$solve(prob, saveat=saveat_fit)
mat_SHAM <- sapply(sol_SHAM$u,identity)
udf_SHAM <- as.data.frame(t(mat_SHAM))

p_fit <- p_VNS_WT
u0 <- c(0.75,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_VNS <- de$solve(prob, saveat=saveat_fit)
mat_VNS <- sapply(sol_VNS$u,identity)
udf_VNS <- as.data.frame(t(mat_VNS))


plot(SHAM_alpha$tt, SHAM_alpha$nf, ylim = c(0,30), xlim = c(0,50), col ="red", ylab = "Neutrofill count", xlab = "Time (h)")
points(VNS_alpha$tt, VNS_alpha$nf, col = "blue")
lines(sol_SHAM$t, udf_SHAM$V1, col = "red")
lines(sol_VNS$t, udf_VNS$V1, col = "blue")
legend(40, 30,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)

plotA7 <- plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))

pdf("Curvies_A7-mid.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()

pdf("Curvies_A7-big.pdf", height = 5, width = 5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()


pdf("Curvies_A7-variants_mid2.pdf", height = 3.5, width = 10.5)
gridExtra::grid.arrange(
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()

pdf("Curvies_A7-variants_big2.pdf", height = 5, width = 15)
gridExtra::grid.arrange(
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()


# Alox

tspan <- list(0,50)
saveat_fit = c(0:500)/10

p_fit <- p_SHAM_alox[1:4]
u0 <- c(0.1,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_SHAM <- de$solve(prob, saveat=saveat_fit)
mat_SHAM <- sapply(sol_SHAM$u,identity)
udf_SHAM <- as.data.frame(t(mat_SHAM))

p_fit <- p_VNS_alox[1:4]
u0 <- c(0.1,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_VNS <- de$solve(prob, saveat=saveat_fit)
mat_VNS <- sapply(sol_VNS$u,identity)
udf_VNS <- as.data.frame(t(mat_VNS))


plot(SHAM_alox$tt, SHAM_alox$nf, ylim = c(0,30), xlim = c(0,50), col ="red", ylab = "Neutrofill count", xlab = "Time (h)")
points(VNS_alox$tt, VNS_alox$nf, col = "blue")
lines(sol_SHAM$t, udf_SHAM$V1, col = "red")
lines(sol_VNS$t, udf_VNS$V1, col = "blue")
legend(40, 30,   # Coordinates (x also accepts keywords)
       legend = c("SHAM", "VNS"), # Vector with the name of each group
       fill = c("red", "blue"),   # Creates boxes in the legend with the specified colors
       col = par("red", "blue"), # Color of lines or symbols  
)


plotAlox <- plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))

pdf("Curvies_Alox-mid.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()

pdf("Curvies_Alox-big.pdf", height = 5, width = 5)
plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))
dev.off()


pdf("Curvies_Alox-variants_mid2.pdf", height = 3.5, width = 10.5)
gridExtra::grid.arrange(
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()

pdf("Curvies_Alox-variants_big2.pdf", height = 5, width = 15)
gridExtra::grid.arrange(
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35)),
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL),
  plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, yLimits = NULL, addMeans = TRUE),
  ncol = 3, nrow = 1
)
dev.off()



plotWT <- gridExtra::arrangeGrob(grobs = plotWT, left = grid::textGrob("WT", rot = 90))
plotAlox <- gridExtra::arrangeGrob(grobs = plotAlox, left = grid::textGrob("Alox15 KO", rot = 90))
plotA7 <- gridExtra::arrangeGrob(grobs = plotA7, left = grid::textGrob("Chrna7 KO", rot = 90))

pdf(file.path(getwd(), "Curvies-big.pdf"), width = 5.2, height = 15)
gridExtra::grid.arrange(plotWT, plotAlox, plotA7, nrow = 3, ncol = 1)
dev.off()


#Add additional lines
Plot <- plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = c(0, 35))

plotCurvies <- function(pointsA, pointsB, linesTimeA, linesNeuA, linesTimeB, linesNeuB, 
                        nameA = "Sham", nameB = "VNS", xJitter = 1.25, xBreak = 2.5, yLimits = NULL, addMeans = FALSE, addLines = FALSE, lineList = NULL) {
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
                  color = "black", linetype = "dashed")
      
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
                  color = "black", size = 1.5, arrow = arrow(type = "closed", ends = "first", length = unit(0.1, "inches")))
      
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
                  color = "black", linetype = "dashed")
      
      #Vertical
      lineDf <- lineList %>%
        do.call("rbind", .) %>%
        tibble::as_tibble()
      lineDf$Treatment <- names(lineList)
      
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
                  color = "black", linetype = "dashed")
      
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
                  color = "black", arrow = arrow(type = "closed", ends = "both", length = unit(0.075, "inches")))
    }
  }
  
  return(Plot)
}


SHAMnums <- c(
  Val_T_max = 12.69608,
  T_max = 5.3,
  `Val50%T_max` = 6.348795,
  `Time_below_50%_Val` = 20.9,
  Resolution_time = 15.6
)

VNSnums <- c(
  Val_T_max  = 11.69759,
  T_max = 4.3,
  `Val50%T_max` = 5.84804,
  `Time_below_50%_Val` = 15.6,
  Resolution_time = 12.3
)

lineList <- list(
  Sham = SHAMnums,
  VNS = VNSnums
)

pdf("Curvies_WT-variants_mid3.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_WT-variants_big3.pdf", height = 5, width = 5)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_WT-variants_verybig3.pdf", height = 7.5, width = 7.5)
plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()