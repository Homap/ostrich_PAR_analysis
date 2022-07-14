#!/usr/bin/R

# Clean the workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(tidyquant)
library(cowplot)

# Figure 1 contain cM/Mb, rho and LD

mb = 10^6
kb = 10^3
#------------
# Genetic map
#------------
female_map <- read.table("../../data/geneticmap/female.1MB.window.Z.txt", header = T)
male_map <- read.table("../../data/geneticmap/male.1MB.window.Z.txt", header = T)
sex_averaged_map <- read.table("../../data/geneticmap/sex_averaged.1MB.window.Z.txt", header = T)

map <- data.frame(
  position=c(rep((female_map$Window_start+female_map$Window_end)/(2*mb),3)),
  rate=c(female_map$CM_per_bp*mb, male_map$CM_per_bp*mb, sex_averaged_map$CM_per_bp*mb),
  sex=c(rep('female', 80), rep('male', 80), rep('sex_averaged', 80))
)

map_plot <- ggplot(map , aes(x = position, y = rate, group = sex)) +
  geom_point(aes(colour=sex), size = 0.8, alpha = 0.7)+ geom_ma(aes(color = sex), ma_fun = SMA, n = 5, na.rm = T, size = 0.8, linetype = 6) + theme_classic()+
  ylab(expression(over(cM, Mb))) + xlim(0,82) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  scale_colour_manual(values=c("red", "blue","purple")) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5))

#------------
# rho
#------------
rho_data <- read.table("../../data/rho/ldhat_rho/z/200Kb50Kb_rho_Z.chr.coord.txt", header=F)
colnames(rho_data) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")
rho_dataset <- data.frame(pos = (rho_data$Window_start+rho_data$Window_end)/(2*mb), rho_per_window=rho_data$rho_per_window/kb)

rho_plot <- ggplot(rho_dataset , aes(x = (rho_data$Window_start+rho_data$Window_end)/(2*mb), y = rho_per_window)) +
  geom_point(size = 0.8, alpha = 0.9, color = "grey")+ geom_ma(ma_fun = SMA, n = 5, size = 0.8, linetype = 6, na.rm = T, color = "purple") + theme_classic()+
  ylab(expression(over(4~N[e],Kb))) + geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  theme(legend.position = "none") +xlim(0,82) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1))

#------------
# LD
#------------
ld_data <- read.table("../../data/ld/ld_window/z/Z.LD.05-500.200kbBin.150.tab")
colnames(ld_data) <- c("chr", "Window_start", "Window_end", "Scaffold", "pos1", "pos2", "Dprime", "LOD", "r2", "CIlow", "CIhi", "meanDist", "overlap")
ld_dataset <- data.frame(pos = (pos = ld_data$Window_start+rho_data$Window_end)/(2*mb), r2=ld_data$r2)

ld_plot <- ggplot(ld_dataset , aes(x = pos, y = r2)) +
  geom_point(size = 0.8, alpha = 0.9, color = "grey")+ geom_ma(ma_fun = SMA, n = 5, size = 0.8, linetype = 6, na.rm = T, color = "purple") + theme_classic()+
  xlab("Position (Mb)") + ylab(expression(r^2)) + ylim(0.08, 0.3) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) + scale_x_continuous(breaks=seq(0,82,10))

#-----------
# rho vs. LD
#-----------
rho_data_par <- rho_data[rho_data$Window_start < 52185292,]
rho_data_nonpar <- rho_data[rho_data$Window_start > 52185292,]

ld_data_par <- ld_data[ld_data$Window_start < 52185292,]
ld_data_nonpar <- ld_data[ld_data$Window_start > 52185292,]

ld_rho_dataset <- data.frame(rho <- c(rho_data_par$rho_per_window, rho_data_nonpar$rho_per_window),
                             ld <- c(ld_data_par$r2, ld_data_nonpar$r2),
                             chr <- c(rep("PAR", 346), rep("nonPAR", 188)))

rho_ld_plot <- ggplot(ld_rho_dataset, aes(x = rho, y = ld, group = chr)) + ylab(expression(r^2)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression(4~N[e]/Kb))+
  theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_colour_manual(values=c("azure4", "black"))

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

ld_decay_dataset <- data.frame(dist = rep(ld_decay_autosome$Dist/kb,3),
  r2 = c(ld_decay_autosome$Mean_r2, ld_decay_par$Mean_r2, ld_decay_nonpar$Mean_r2),
  chr = c(rep("Autosome", 3010), rep("PAR", 3010), rep("SDR", 3010)))

ld_decay_plot <- ggplot(data = ld_decay_dataset, aes(x = dist, y = r2, group = chr)) + xlim(0, 50) + ylab(expression(r^2)) +
geom_line(aes(color = chr)) + theme_classic() + theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
scale_colour_manual(values=c("darkorange", "azure4","black")) + scale_y_continuous(breaks=seq(0.1,0.45,0.05)) + xlab("Distance (Kb)")

#----------
# LD decay - end PAR, mid PAR, PAR boundary
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

ld_decay_dataset <- data.frame(dist = rep(ld_decay_autosome$Dist/kb,3),
                               r2 = c(ld_decay_autosome$Mean_r2, ld_decay_par$Mean_r2, ld_decay_nonpar$Mean_r2),
                               chr = c(rep("Autosome", 3010), rep("PAR", 3010), rep("SDR", 3010)))

ld_decay_plot <- ggplot(data = ld_decay_dataset, aes(x = dist, y = r2, group = chr)) + xlim(0, 50) + ylab(expression(r^2)) +
  geom_line(aes(color = chr)) + theme_classic() + theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_colour_manual(values=c("darkorange", "azure4","black")) + scale_y_continuous(breaks=seq(0.1,0.45,0.05)) + xlab("Distance (Kb)")

#------------
# Putting three figures together
#------------
figure1 <- ggdraw(xlim = c(0, 1), ylim = c(0,3.2)) +
  draw_plot(map_plot, x = 0.05, y = 2.45, width = 0.90, height = .68) +
  draw_plot(rho_plot, x = 0.035, y = 1.79, width = 0.92, height = .68) +
  draw_plot(ld_plot, x = 0.056, y = 1, width = 0.89, height = .8) +
  draw_plot(rho_ld_plot, x = 0.07, y = 0, width = 0.41, height = 1) +
  draw_plot(ld_decay_plot, x = 0.52, y = 0.015, width = 0.41, height = 1) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 12,
                 x = c(0.15, 0.15, 0.148, 0.16, 0.61), y = c(3.15, 2.49, 1.82, 1.05, 1.05))
ggsave("../../figures/Figure1.pdf", figure1, width = 5.89, height = 6.72, units = "in")









