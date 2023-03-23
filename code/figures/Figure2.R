#!/usr/bin/R

# Clean the workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(tidyquant)
library(patchwork)

# Figure 1 contain cM/Mb, rho and LD

mb = 10^6
kb = 10^3

#------------
# rho
#------------
rho_data <- read.table("../../data/rho/ldhat_rho/z/200Kb50Kb_rho_Z.chr.coord.txt", header=F)
colnames(rho_data) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")
rho_dataset <- data.frame(pos = (rho_data$Window_start+rho_data$Window_end)/(2*mb), rho_per_window=rho_data$rho_per_window/kb)

rho_chr4 <- read.table("../../data/rho/ldhat_rho/chr4/superscaffold11.200kb200Kb.rho.txt", header = T)
rho_chr5 <- read.table("../../data/rho/ldhat_rho/chr5/superscaffold8.200kb200Kb.rho.txt", header = T)
rho_a <- rbind(rho_chr4, rho_chr5)
mean_rho_a <- mean(rho_a$rho_per_window/kb) # 0.1295445
# Average rho for PAR

rho_plot <- ggplot(rho_dataset , aes(x = (rho_data$Window_start+rho_data$Window_end)/(2*mb), y = rho_per_window)) +
  geom_point(size = 0.8, alpha = 0.9, color = "grey")+ geom_ma(ma_fun = SMA, n = 5, size = 0.8, linetype = 6, na.rm = T, color = "black") + theme_classic()+
  ylab(expression(over(rho,Kb))) + geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=mean_rho_a, linetype="dashed", color = "#E69F00") +
  theme(legend.position = "none") +scale_x_continuous(breaks=seq(0,82,10))+
  theme(axis.text.y=element_text(size=12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) + ggtitle('A')

#------------
# LD
#------------
ld_data <- read.table("../../data/ld/ld_window/z/Z.LD.05-500.200kbBin.150.tab")
colnames(ld_data) <- c("chr", "Window_start", "Window_end", "Scaffold", "pos1", "pos2", "Dprime", "LOD", "r2", "CIlow", "CIhi", "meanDist", "overlap")
ld_dataset <- data.frame(pos = (pos = ld_data$Window_start+rho_data$Window_end)/(2*mb), r2=ld_data$r2)

ld_a <- read.table("../../data/ld/ld_window/autosome/chr4_5_pairwise_LD.tab", header = T)
mean_ld_a <- mean(ld_a$r2) # 0.1496411

ld_plot <- ggplot(ld_dataset , aes(x = pos, y = r2)) +
  geom_point(size = 0.8, alpha = 0.9, color = "grey")+ geom_ma(ma_fun = SMA, n = 5, size = 0.8, linetype = 6, na.rm = T, color = "black") + theme_classic()+
  ylab(expression(r^2)) + ylim(0.08, 0.3) + scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept = mean_ld_a, linetype="dashed", color = "#E69F00") +
  theme(axis.text.y=element_text(size=12),  axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5, size = 12), legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,-5,-20)) + ggtitle('B')


#------------
# Genetic map
#------------
female_map <- read.table("../../data/geneticmap/female.1MB.window.Z.txt", header = T)
male_map <- read.table("../../data/geneticmap/male.1MB.window.Z.txt", header = T)
sex_averaged_map <- read.table("../../data/geneticmap/sex_averaged.1MB.window.Z.txt", header = T)

map <- data.frame(
  position=c(rep((female_map$Window_start+female_map$Window_end)/(2*mb),3)),
  rate=c(female_map$CM_per_bp*mb, male_map$CM_per_bp*mb, sex_averaged_map$CM_per_bp*mb),
  sex=c(rep('Female', 80), rep('Male', 80), rep('Sex averaged', 80))
)

map_plot <- ggplot(map , aes(x = position, y = rate, group = sex)) +
  geom_point(aes(colour=sex), size = 0.8, alpha = 0.7)+ geom_ma(aes(color = sex), ma_fun = SMA, n = 5, na.rm = T, size = 0.8, linetype = 6) + theme_classic()+
  xlab("Position (Mb)") + ylab(expression(over(cM, Mb))) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  scale_colour_manual(values=c("red", "blue","black"), name = NULL) +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), legend.text = element_text(size=12), legend.position = "right", legend.justification="top", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) + scale_x_continuous(breaks=seq(0,82,10)) + ggtitle('C')


#-----------
# rho vs. cM/Mb
#-----------
rho_mb <- read.table("../../data/rho/ldhat_rho/z/1Mb1Mb.rho.chrom.Z.txt", header=T)
female_map <- read.table("../../data/geneticmap/female.1MB.window.Z.txt", header=T)
male_map <- read.table("../../data/geneticmap/male.1MB.window.Z.txt", header=T)
sex_averaged_map <- read.table("../../data/geneticmap/sex_averaged.1MB.window.Z.txt", header=T)

rho_mb_par <- rho_mb[rho_mb$Window_start<52000000,]
#rho_mb_nonpar <- rho_mb[rho_mb$Window_start>52000000,]

female_map_par <- female_map[female_map$Window_start<52000000,]
male_map_par <- male_map[male_map$Window_start<52000000,]
sex_averaged_map_par <- sex_averaged_map[sex_averaged_map$Window_start<52000000,]
#sex_averaged_map_nonpar <- sex_averaged_map[sex_averaged_map$Window_start>52000000,]

rho_cm_mb <- data.frame(cm = c(female_map_par$CM_per_bp[2:52]*mb, male_map_par$CM_per_bp[2:52]*mb, sex_averaged_map_par$CM_per_bp[2:52]*mb),
                        rho = c(rho_mb_par$rho_per_window[2:52]/kb, rho_mb_par$rho_per_window[2:52]/kb, rho_mb_par$rho_per_window[2:52]/kb),#, rho_mb_nonpar$rho_per_window[1:26]/kb),
                        chr = c(rep("Female", 51), rep("Male", 51), rep("Sex averaged", 51)))#, rep("nonPAR", 26)))

#rho_cm_mb <- data.frame(cm = c(female_map_par$CM_per_bp[2:52]*mb, sex_averaged_map_par$CM_per_bp[2:52]*mb),
#                        rho = c(rho_mb_par$rho_per_window[2:52]/kb, rho_mb_par$rho_per_window[2:52]/kb),#, rho_mb_nonpar$rho_per_window[1:26]/kb),
#                        chr = c(rep("Female", 51), rep("Sex averaged", 51)))#, rep("nonPAR", 26)))

rho_cm_plot <- ggplot(rho_cm_mb, aes(x = cm, y = rho, color = chr)) + geom_point(size = 1.3, alpha = 0.5) + theme_classic() + #, aes(color = chr)) + theme_classic() +
  xlab("cM/Mb") + scale_colour_manual(values=c("red", "blue", "black"), name = NULL) +
  scale_y_continuous(breaks=seq(0,3,0.5)) + ylab(expression(over(rho,Kb))) +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), legend.text = element_text(size=12), legend.position = "right", legend.justification = "top", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12), legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,-5,-15)) + ggtitle('D') + geom_smooth(method='lm', aes(fill = chr), fullrange = TRUE, se=FALSE, show.legend = FALSE)


#-----------
# rho vs. LD
#-----------
rho_data_par <- rho_data[rho_data$Window_start < 52185292,]
rho_data_nonpar <- rho_data[rho_data$Window_start > 52185292,]

ld_data_par <- ld_data[ld_data$Window_start < 52185292,]
ld_data_nonpar <- ld_data[ld_data$Window_start > 52185292,]

ld_rho_dataset <- data.frame(rho <- c(rho_data_par$rho_per_window, rho_data_nonpar$rho_per_window),
                             ld <- c(ld_data_par$r2, ld_data_nonpar$r2),
                             chr <- c(rep("PAR", 346), rep("SLR", 188)))

rho_ld_plot <- ggplot(ld_rho_dataset, aes(x = rho/kb, y = ld, color = chr)) + ylab(expression(r^2)) +
  geom_point(aes(color = chr), size = 1.3, alpha = 0.5) + theme_classic() + xlab(expression(over(rho,Kb)))+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), legend.text = element_text(size=12), legend.position = "right", legend.justification = "top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,-5,-15), axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  scale_colour_manual(values=c("black", "grey"), name = NULL)+ ggtitle('E')

#----------
# LD decay
#----------
# Autosome
ld_decay_autosome <- read.table("../../data/ld/ld_decay/autosome/chr4_5.bin")
colnames(ld_decay_autosome) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# PAR
ld_decay_par <- read.table("../../data/ld/ld_decay/z/par/PAR.bin")
colnames(ld_decay_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# nonPAR
ld_decay_nonpar <- read.table("../../data/ld/ld_decay/z/nonpar/nonPAR.bin")
colnames(ld_decay_nonpar) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# start_PAR
ld_decay_start_par <- read.table("../../data/ld/ld_decay/z/par/start.par.5Mb.bin")
colnames(ld_decay_start_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# mid_PAR
ld_decay_mid_par <- read.table("../../data/ld/ld_decay/z/par/mid.par.5Mb.bin")
colnames(ld_decay_mid_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# end_PAR
ld_decay_end_par <- read.table("../../data/ld/ld_decay/z/par/end.par.5Mb.bin")
colnames(ld_decay_end_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")

ld_decay_dataset <- data.frame(dist = rep(ld_decay_autosome$Dist/kb,6),
  r2 = c(ld_decay_autosome$Mean_r2, ld_decay_par$Mean_r2, ld_decay_nonpar$Mean_r2, ld_decay_start_par$Mean_r2, ld_decay_mid_par$Mean_r2, ld_decay_end_par$Mean_r2),
  chr = c(rep("Autosome", 3050), rep("PAR", 3050), rep("SLR", 3050), rep("PAR-1-5Mb", 3050), rep("midPAR", 3050), rep("PAR-boundary", 3050)))

ld_decay_dataset$chr <-factor(ld_decay_dataset$chr, levels=c("Autosome","PAR","SLR","PAR-1-5Mb","midPAR", "PAR-boundary"))

ld_decay_plot <- ggplot(data = ld_decay_dataset, aes(x = dist, y = r2, color = chr)) + xlim(0, 50) + ylab(expression(r^2)) +
geom_line(size = 0.4, alpha = 0.5) + theme_classic() + theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), legend.text = element_text(size=12), legend.position = "right", legend.justification = "top", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12), legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,-5,-8)) + geom_ma(aes(color = chr), ma_fun = SMA, n = 5, na.rm = T, size = 0.8, linetype = 6) +
scale_colour_manual(breaks = c("SLR", "midPAR", "Autosome","PAR-boundary", "PAR", "PAR-1-5Mb"), values=c("#999999", "#D55E00", "#E69F00", "#CC79A7", "#000000", "#0072B2"), name = NULL) + scale_y_continuous(breaks=seq(0.1,0.5,0.1)) + xlab("Distance (Kb)") + ggtitle('F')

#------------
# Putting six figures together
#------------

figure1 <- (rho_plot / ld_plot / map_plot) /(rho_cm_plot | rho_ld_plot | ld_decay_plot)
ggsave("../../figures/Figure2.pdf", figure1, width = 10.55, height = 10.75, units = "in")









