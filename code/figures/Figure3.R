#------------------------------------------------------------
# Homa Papoli
# Created on: May 2022
# Description: Script to create plots of Tajima's D distributions
# for Italian sparrow population based on thetas Analysis from ANGSD
# Note: Before running, set the working directory to the root of the project
#------------------------------------------------------------
rm(list = ls())
library(cowplot)
# Z chromosome coordinates
black_Z <- read.table("data/diversity/black.Z.coordinates.sfs.txt", header=F)
colnames(x = black_Z) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")
td = ggplot(black_Z, aes(y=Td, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  xlab("Position (Mb)") + ylab("Tajima's D") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") 

mfFST <- read.table("data/fst/black_male_female_Z_100Kb.windowed.weir.Z.coord.fst", header=F)
mfFST$V8[mfFST$V8 < 0] <- 0
par_mfFST <- mfFST[(mfFST$V5<52193205),]
nonpar_mfFST <- mfFST[(mfFST$V5>52193205),]
mean(par_mfFST$V8)
mean(nonpar_mfFST$V8)
t.test(par_mfFST$V8, nonpar_mfFST$V8)
hist(par_mfFST$V8)
hist(nonpar_mfFST$V8)


fst = ggplot(mfFST, aes(y=mfFST$V8, x=(mfFST$V5+mfFST$V6)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  xlab("Position (Mb)") + ylab("FST") +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") 

ggdraw(ylim = c(0,0.9)) +
  draw_plot(td, x = 0.02, y = 0.5, width = 0.98, height = 0.4) +
  draw_plot(fst, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 12,
                  x = c(0, 0), y = c(0.9,0.5))
