#------------------------------------------------------------
# Homa Papoli
# Created on: May 2022
# Description: Script to create plots of Tajima's D distributions
# for Italian sparrow population based on thetas Analysis from ANGSD
# Note: Before running, set the working directory to the root of the project
#------------------------------------------------------------
rm(list = ls())
rho <- read.table("data/rho/pos_rho.txt", header=T)
cm <- read.table("data/geneticmap/pos_cm_smooth.txt", header=T)

plot(rho$positions, rho$rho_rates_kb, type="s", col=rgb(0,0,0.5), xlab = "Position (Mb)", ylab = "4Ner/Kb", las = 1)
par(new=TRUE)
plot(cm$positions, cm$sex_averaged_cM, col = "darkorchid1", xlab="", ylab="", ylim = c(0, 65), pch = 20, axes = FALSE)
lines(cm$positions, cm$sex_averaged_smooth, col = "darkorchid1")
mtext("cM",side=4,col="red",line=4) 
axis(4, ylim=c(0,95),las=1)

legend("topleft",legend=c("4Ner","cM"),
       text.col=c("black","darkorchid1"),pch=c(20),lty=1:2,col=c("black","darkorchid1"), box.lty=0)
