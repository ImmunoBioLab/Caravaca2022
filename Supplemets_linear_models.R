#Linear models of wild type mice
#Import data using Import Dataset
df <- neutrofil_count_Hildur
df <- df[5:nrow(df),4:8]
colnames(df) <- c("Four","Twelwe","Twentyfour","Fourtyeights", "Seventytwo")

#Organize data
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

#Regression
res2 <- lapply(1:n, function(index) {
  experiment[sstart[index]:(sstart[index+1]-1),] %>%
    aprilReg2(.)
})

#Experiments to keep by: 9, 11, 13, 14, 15, 17, 18, 19, 21, 22
res2 %<>% .[c(9, 11, 13, 14, 15, 17, 18, 19, 21, 22)]

#Regression line graph
regPlots <- lapply(seq_along(res2), function(index) {
  print(index) 
  regPlot2(res2[[index]], yLims = c(-4, 4), expIndex = index)
})

pdf(file.path(getwd(), "hall.pdf"), width = 20, height = 20)
gridExtra::grid.arrange(grobs = regPlots)
dev.off()

pdf(file.path(getwd(), "newHall2.pdf"), width = 20, height = 10)
gridExtra::grid.arrange(grobs = regPlots, ncol = 5, nrow = 2)
dev.off()