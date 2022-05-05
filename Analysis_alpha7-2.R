library("ggplot2")



# Preliminaries

# read data
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
mufun <- function(s, mu_max,K){
  return(mu_max*s/(s+K))
}

f <- function(u,p,t) {
  du1 = (mufun(u[2],p[2],p[4])-p[3])*u[1]
  du2 = -p[1]* mufun(u[2],p[2],p[4])*u[2]
  return(c(du1,du2))
}


tspan <- list(1,48)
saveat = c(4,12,24,48)


# Analysis 

DD<-SHAM_WT
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

u0 <- c(0.1,1)
res_SHAM_alpha <- optim(p,loglike, DD = SHAM_alpha)
p_SHAM_alpha <- res_SHAM_alpha$par

u0 <- c(1,1)
res_VNS_alpha <- optim(p,loglike, DD = VNS_alpha)
p_VNS_alpha <- res_VNS_alpha$par

u0 <- c(0.1,1)
res_SHAM_alox <- optim(p,loglike, DD = SHAM_alox)
p_SHAM_alox <- res_SHAM_alox$par
res_VNS_alox <- optim(p,loglike, DD = VNS_alox)
p_VNS_alox <- res_VNS_alox$par



# Plotting results

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


# Confidence intervals


#data_tmp <- SHAM_WT$nf
#times_tmp <- SHAM_WT$tt
data_tmp <- VNS_WT$nf
times_tmp <- VNS_WT$tt



#type = "SHAM"
type = "VNS"

data4h <- c()
data12h <- c()
data24h <- c()
data48h <- c()

for(i in 1:length(times_tmp)){
  if(times_tmp[i] == 4)
    data4h <- c(data4h, data_tmp[i])
  else if(times_tmp[i] == 12)
    data12h <- c(data12h, data_tmp[i])
  else if(times_tmp[i] == 24)
    data24h <- c(data24h, data_tmp[i])
  else if(times_tmp[i] == 48)
    data48h <- c(data48h, data_tmp[i])
}
t <- c(4,12,24,48)
m <- c(mean(data4h), mean(data12h), mean(data24h), mean(data48h))
lo <- m -qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
up <- m+qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
# Vertical arrow
if(type == "SHAM"){colfig <- "red"}
if(type == "VNS"){colfig <- "blue"}

for(i in 1:4){
  arrows(x0=t[i], y0=lo[i], x1=t[i], y1=up[i], code=3, angle = 90, length = 0.2, col=colfig, lwd=2)
}

#stats
max(udf_SHAM$V1)
sol_SHAM$t[which.max(udf_SHAM$V1)]
half <- max(udf_SHAM$V1)/2
rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]
res_time <- rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]-sol_SHAM$t[which.max(udf_SHAM$V1)]
res_time

max(udf_VNS$V1)
sol_VNS$t[which.max(udf_VNS$V1)]
half <- max(udf_VNS$V1)/2
rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]
res_time <- rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]-sol_VNS$t[which.max(udf_VNS$V1)]
res_time


# Alpha

tspan <- list(0,50)
saveat_fit = c(0:500)/10

p_fit <- p_SHAM_alpha
u0 <- c(0.1,1)
prob <- de$ODEProblem(f, u0, tspan, p_fit)
sol_SHAM <- de$solve(prob, saveat=saveat_fit)
mat_SHAM <- sapply(sol_SHAM$u,identity)
udf_SHAM <- as.data.frame(t(mat_SHAM))

#p_fit <- p_VNS_alpha
p_fit <- p_SHAM_alpha
p_fit[3] <- 0.3
p_fit[2] <- 18
p_fit[1] <- 0.4
u0 <- c(0.1,1)
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


# Confidence intervals


#data_tmp <- SHAM_alpha$nf
#times_tmp <- SHAM_alpha$tt
#type = "SHAM"

data_tmp <- VNS_alpha$nf
times_tmp <- VNS_alpha$tt
type = "VNS"





data4h <- c()
data12h <- c()
data24h <- c()
data48h <- c()

for(i in 1:length(times_tmp)){
  if(times_tmp[i] == 4)
    data4h <- c(data4h, data_tmp[i])
  else if(times_tmp[i] == 12)
    data12h <- c(data12h, data_tmp[i])
  else if(times_tmp[i] == 24)
    data24h <- c(data24h, data_tmp[i])
  else if(times_tmp[i] == 48)
    data48h <- c(data48h, data_tmp[i])
}
t <- c(4,12,24,48)
m <- c(mean(data4h), mean(data12h), mean(data24h), mean(data48h))
lo <- m -qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
up <- m+qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
# Vertical arrow
if(type == "SHAM"){colfig <- "red"}
if(type == "VNS"){colfig <- "blue"}

for(i in 1:4){
  arrows(x0=t[i], y0=lo[i], x1=t[i], y1=up[i], code=3, angle = 90, length = 0.1, col=colfig, lwd=2)
}


#stats
max(udf_SHAM$V1)
sol_SHAM$t[which.max(udf_SHAM$V1)]
half <- max(udf_SHAM$V1)/2
half
rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]
res_time <- rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]-sol_SHAM$t[which.max(udf_SHAM$V1)]
res_time

max(udf_VNS$V1)
sol_VNS$t[which.max(udf_VNS$V1)]
half <- max(udf_VNS$V1)/2
half
rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]
res_time <- rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]-sol_VNS$t[which.max(udf_VNS$V1)]
res_time



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




# Confidence intervals


#data_tmp <- SHAM_alox$nf
#times_tmp <- SHAM_alox$tt
#type = "SHAM"

data_tmp <- VNS_alox$nf
times_tmp <- VNS_alox$tt
type = "VNS"





data4h <- c()
data12h <- c()
data24h <- c()
data48h <- c()

for(i in 1:length(times_tmp)){
  if(times_tmp[i] == 4)
    data4h <- c(data4h, data_tmp[i])
  else if(times_tmp[i] == 12)
    data12h <- c(data12h, data_tmp[i])
  else if(times_tmp[i] == 24)
    data24h <- c(data24h, data_tmp[i])
  else if(times_tmp[i] == 48)
    data48h <- c(data48h, data_tmp[i])
}
t <- c(4,12,24,48)
m <- c(mean(data4h), mean(data12h), mean(data24h), mean(data48h))
lo <- m -qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
up <- m+qnorm(0.975)*c(sd(data4h)/sqrt(length(data4h)), sd(data12h)/sqrt(length(data12h)), sd(data24h)/sqrt(length(data24h)), sd(data48h)/sqrt(length(data48h)))
# Vertical arrow
if(type == "SHAM"){colfig <- "red"}
if(type == "VNS"){colfig <- "blue"}

for(i in 1:4){
  arrows(x0=t[i], y0=lo[i], x1=t[i], y1=up[i], code=3, angle = 90, length = 0.1, col=colfig, lwd=2)
}


#stats
max(udf_SHAM$V1)
sol_SHAM$t[which.max(udf_SHAM$V1)]
half <- max(udf_SHAM$V1)/2
half
rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]
res_time <- rev(sol_SHAM$t)[which.max(rev(udf_SHAM$V1) > half)]-sol_SHAM$t[which.max(udf_SHAM$V1)]
res_time

max(udf_VNS$V1)
sol_VNS$t[which.max(udf_VNS$V1)]
half <- max(udf_VNS$V1)/2
half
rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]
res_time <- rev(sol_VNS$t)[which.max(rev(udf_VNS$V1) > half)]-sol_VNS$t[which.max(udf_VNS$V1)]
res_time


#Add additional lines
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


#alpha
SHAMnums <- c(
  Val_T_max = 7.92775,
  T_max = 9,
  `Val50%T_max` = 3.963875,
  `Time_below_50%_Val` = 17.9,
  Resolution_time = 8.9
)

VNSnums <- c(
  Val_T_max  = 9.784077,
  T_max = 8.1,
  `Val50%T_max` = 4.892039,
  `Time_below_50%_Val` = 17.3,
  Resolution_time = 7.8
)

lineList <- list(
  Sham = SHAMnums,
  VNS = VNSnums
)

plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, 
            yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)

pdf("Curvies_alpha-variants_mid3.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_alpha-variants_big3.pdf", height = 5, width = 5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_alpha-variants_verybig3.pdf", height = 7.5, width = 7.5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()


pdf("Curvies_alpha-variants_mid4.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = FALSE, lineList = NULL)
dev.off()

pdf("Curvies_alpha-variants_big4.pdf", height = 5, width = 5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = FALSE, lineList = NULL)
dev.off()

pdf("Curvies_alpha-variants_verybig4.pdf", height = 7.5, width = 7.5)
plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, yLimits = NULL, addMeans = TRUE, addLines = FALSE, lineList = NULL)
dev.off()

#alox
SHAMnums <- c(
  Val_T_max = 3.495074,
  T_max = 10.4,
  `Val50%T_max` = 1.747537,
  `Time_below_50%_Val` = 30.2,
  Resolution_time = 19.8
)

VNSnums <- c(
  Val_T_max  = 2.30537,
  T_max = 12.3,
  `Val50%T_max` = 1.152685,
  `Time_below_50%_Val` = 21.4,
  Resolution_time = 16.0
)

lineList <- list(
  Sham = SHAMnums,
  VNS = VNSnums
)

plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, 
            yLimits = c(0, 6), addMeans = TRUE, addLines = TRUE, lineList = lineList)

pdf("Curvies_alox-variants_mid3.pdf", height = 3.5, width = 3.5)
plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, 
            yLimits = c(0, 6), addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_alox-variants_big3.pdf", height = 5, width = 5)
plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, 
            yLimits = c(0, 6), addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()

pdf("Curvies_alox-variants_verybig3.pdf", height = 7.5, width = 7.5)
plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, 
            yLimits = c(0, 6), addMeans = TRUE, addLines = TRUE, lineList = lineList)
dev.off()


#New line color
wtPlot <- plotCurvies(SHAM_WT, VNS_WT, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, 
                      yLimits = c(0, 15), addMeans = TRUE, addLines = TRUE, lineList = lineList, lineColor = "gray60")
aloxPlot <- plotCurvies(SHAM_alox, VNS_alox, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 1, 
                        yLimits = c(0, 6), addMeans = TRUE, addLines = TRUE, lineList = lineList, lineColor = "gray60")
alphaPlot <- plotCurvies(SHAM_alpha, VNS_alpha, sol_SHAM$t, udf_SHAM$V1, sol_VNS$t, udf_VNS$V1, xBreak = 2.5, 
                         yLimits = NULL, addMeans = TRUE, addLines = TRUE, lineList = lineList, lineColor = "gray60")

pdf("Curvies_All_grayLines-1.pdf", height = 5, width = 15)
gridExtra::grid.arrange(grobs = list(wtPlot, aloxPlot, alphaPlot), nrow = 1, ncol = 3)
dev.off()

pdf("Curvies_All_grayLines-1_mid.pdf", height = 3.5, width = 10.5)
gridExtra::grid.arrange(grobs = list(wtPlot, aloxPlot, alphaPlot), nrow = 1, ncol = 3)
dev.off()