#Linear regression models of wild type (WT), Alox15 knock-out (alox) and CHRNA7 knock-out (a7) mice
#Load data
df <- read.xlsx("April_nov_2021_read.xlsx")

#WT
WT_SH <- data.frame(df[[2]], df$SHAM_WT)
WT_SH <- WT_SH[5:22,]
colnames(WT_SH) <- c("Time", "SHAM_WT")

WT_VNS <- data.frame(df[[4]], df$VNS_WT)
WT_VNS <- WT_VNS[4:21,]
colnames(WT_VNS) <- c("Time", "VNS_WT")

#Perform regression analysis
wt <- aprilReg(WT_SH, WT_VNS)

#Plot
pdf(file.path(getwd(), "ggplot_WT.pdf"), width = 5, height = 5)
regPlot(wt)
dev.off()

pdf(file.path(getwd(), "slopes_WT.pdf"), width = 3, height = 5)
plotSlope(wt)
dev.off()


#Alox 15 deficient
Alox_SH <- data.frame(df[[10]], df$SHAM_alox)
Alox_SH <- Alox_SH[3:16,]
colnames(Alox_SH) <- c("Time", "SHAM_alox")

Alox_VNS <- data.frame(df[[12]], df$VNS_alox)
Alox_VNS <- Alox_VNS[3:16,]
colnames(Alox_VNS) <- c("Time", "VNS_alox")

#Perform regression analysis
alox <- aprilReg(Alox_SH, Alox_VNS)

#Plot
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

#Perform regression analysis
a7 <- aprilReg(Alpha_SH, Alpha_VNS)

#Plot
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