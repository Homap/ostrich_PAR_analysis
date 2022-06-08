#*************************************************************************************************
#*Cleaning working directory
#*Loading the required libraries
#*************************************************************************************************
rm(list=ls())
library(dplyr)
library(reshape2) 
library(devtools)
library(ggplot2)
library(gridExtra)
library(data.table)
library(grid)
library(cowplot)
library(ggpubr)
library(tidyquant)
source('software_module/R/recombination_rate.R')
#*************************************************************************************************
#***********************************************************************************************
#***********************************************************************************************











#*************************************************************************************
#* Focus on superscaffold36
PAR_scaffold36 <- black_PAR_all[black_PAR_all$Scaffold=="superscaffold36",]
nonPAR_scaffold36 <- black_nonPAR_all[black_nonPAR_all$Scaffold=="superscaffold36",]
scaffold36 <- black_Z[black_Z$scaffold=="superscaffold36",]
plot((scaffold36$cstart+scaffold36$cend)/2*mb, scaffold36$pi/scaffold36$winsize)
scaf36 <- read.table("../sfs_measures/superscaffold36.10K.Pi.txt", header=T)
head(scaf36)
superscaffold26 =  25310599
superscaffold54 = superscaffold26 + 5000 + 16379243
superscaffold35 = superscaffold54 + 5000 + 4625539
scaf36$midpont = scaf36$midpont + superscaffold35


scaf_36_Pi = ggplot(scaf36, aes(y=rev(scaf36$pi/scaf36$win_length), x=scaf36$midpont/mb)) +
  geom_point(size = 1.3, col="blue", alpha = 0.5) + geom_smooth(colour="blue", size=0.3) + theme_classic() +
  ylab(expression(pi)) + xlim(50, 54) + ylim(0, 0.004) +
  theme(axis.title.x = element_blank()) + 
  geom_vline(xintercept = 52193205/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.001642126, linetype="dashed", color = "black")

scaf_rho <- read.table("/Users/homapapoli/Documents/projects/ostrich_Z/ldhat_dir/rho_results/superscaffold36.black.rho.10.1.1000", header = T)
scaf_rho_plot <- ggplot(scaf_rho, aes(y=rev(scaf_rho$rho_per_site), x=(scaf_rho$SNP1_pos+superscaffold35)/mb)) +
  geom_line(size = 0.5, col="red") + theme_classic() +
  xlab("Position (Mb)") + ylab(expression(rho))

ggplot(scaf_rho, aes(y=cumsum(scaf_rho$rho), x=(scaf_rho$SNP1_pos+superscaffold35)/mb)) +
  geom_line(size = 0.5, col="red") + theme_classic() + 
  xlab("Position (Mb)") + ylab(expression(rho))

grid.newpage()
grid.draw(rbind(ggplotGrob(scaf_36_Pi), ggplotGrob(scaf_rho_plot), size = "last"))


#*************************************************************************************
# Sex averaged cM/Mb
mean(black1Mb_stat$CM_per_bp*10^6, na.rm = T)
sd(black1Mb_stat$CM_per_bp*10^6, na.rm = T)
# PAR
mean(blackPAR1Mb$CM_per_bp*10^6, na.rm = T)
sd(blackPAR1Mb$CM_per_bp*10^6, na.rm = T)
# non-PAR
mean(blacknonPAR1Mb$CM_per_bp*10^6, na.rm = T)
sd(blacknonPAR1Mb$CM_per_bp*10^6, na.rm = T)

# superscaffold36


ggplot(s36_few_rho, aes((((SNP1_pos+SNP2_pos)/2)+superscaffold35 + 5869912)/mb)) + 
  geom_line(aes(y=cumsum(rho_per_site))) +
  geom_line(aes(y=cumsum(s36_few_rho_blue$rho_per_site)), color = "blue") + 
  geom_line(aes(y=cumsum(s36_few_rho_red$rho_per_site)), color = "red") +
  ylab("Cumulative recombination rate (rho)") + xlab ("Position (Mb)")+
  theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text.y = element_blank()) 


library(chromoMap)
chr_file_1 = "s36_scaf.bed"
anno_file_1 = "s36.PAR_nonPAR_genes"
# With segment_annotation = T we get the true length of a gene segment
chromoMap(chr_file_1,anno_file_1)

d <- read.table("Z.genes.chr.coordinates")
plot(d$V2, )

par(mfrow=c(2,1))
plot(c(0,80000000), c(0, 3), type = "n")
rect(d$V5, rep(1,849), d$V6, rep(2,849))

f <- read.table("chicken.Z.coordinates")
plot(c(0,80000000), c(0, 3), type = "n")
rect(f$V2, rep(1,1439), f$V3, rep(2,1439))

setwd("~/Documents/projects/ostrich_Z/results/manuscript_tables/")
s26 <- read.table("superscaffold26.50000.50000.rho.out", header=T)
s54 <- read.table("superscaffold54.50000.50000.rho.out", header=T)
s35 <- read.table("superscaffold35.50000.50000.rho.out", header=T)
s36 <- read.table("superscaffold36.50000.50000.rho.out", header=T)

Ne = 4*58598.67
head(s36)
(sum(s36$rho) + sum(s35$rho) + sum(s54$rho) + sum(s26$rho))/Ne

plot(c(s26$rho_per_site, s54$rho_per_site, s35$rho_per_site, s36$rho_per_site), pch = 20)

blackZ <- read.table("../sfs_measures/black.Z.500Kb.sfs.txt", header=T)
head(blackZ)
plot(blackZ$cstart, blackZ$pi/blackZ$winsize, pch = 20)

#*************************************************************************************
