#------------------------------------------------------------
# Homa Papoli
# Created on: May 2022
# Description: Script to create plots of Tajima's D distributions
# for Italian sparrow population based on thetas Analysis from ANGSD
# Note: Before running, set the working directory to the root of the project
#------------------------------------------------------------
rm(list = ls())
male_boundary_coverage <- read.table("data/coverage/male_boundary_coverage.txt", header = F)
female_boundary_coverage <- read.table("data/coverage/female_boundary_coverage.txt", header = F)

plot(((max(male_boundary_coverage$V2)-male_boundary_coverage$V2)+3510000)/10^6, rowMeans(male_boundary_coverage[,3:7]), pch = 20, col = "blue", 
     ylab = "Mean Depth", xlab = "Position in Mb (superscaffold36)", bty="l", cex = 0.8)
points(((max(female_boundary_coverage$V2)-female_boundary_coverage$V2)+3510000)/10^6, rowMeans(female_boundary_coverage[,3:7]), pch = 20, col = "red", cex = 0.8)
abline(v = ((max(male_boundary_coverage$V2)-3516672)+3510000)/10^6, lty = 2)
abline(v = ((max(male_boundary_coverage$V2)-3524264)+3510000)/10^6, lty = 2)

legend(3.518, 50, legend=c("Male", "Female"), col=c("blue", "red"), pch = c(20, 20), cex=0.8, box.lty=0)
