#!/usr/bin/R

# Clean the workspace
rm(list = ls())

library(ggplot2)
library(tidyquant)
library(patchwork)

mb = 10^6
kb = 10^3

#-----------------------
# pi
#-----------------------
Z_pi <- read.table("../../data/diversity/genetic_variation/z/Z.200000.chr.coord.sfs.txt", header=F)
colnames(x = Z_pi) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")

pi_plot <- ggplot(Z_pi, aes(y=pi/win_length, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "purple", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  xlab("Position (Mb)") + ylab(expression(pi)) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.0016, linetype="dashed", color = "black") + ggtitle('A')

#-----------------------
# TD
#-----------------------
td_plot <- ggplot(Z_pi, aes(y=Td, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "purple", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  xlab("Position (Mb)") + ylab(expression("T"["d"])) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") + ggtitle('B')

#-----------------------
# Fst
#-----------------------
fst <- read.table("../../data/fst/z/black_male_female_Z_200Kb.windowed.weir.Z.coord.fst")
colnames(fst) <- c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start",
                   "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst$WEIGHTED_FST[fst$WEIGHTED_FST<0] <- 0

fst_plot <- ggplot(fst, aes(y=WEIGHTED_FST, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "purple", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  xlab("Position (Mb)") + ylab(expression('M-F F'['st'])) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") + ggtitle('C')


#-----------------------
# pi vs. rho
#-----------------------
rho <- read.table("../../data/rho/ldhat_rho/z/200Kb200Kb_rho_Z.chr.coord.txt")
colnames(rho) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")
GC <- read.table("../../data/genomic_features/z_scaf.200000.intergene.overlap.density.sorted.coord.txt", header=F)
colnames(GC) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                 "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
CDS <- read.table("../../data/genomic_features/z_scaf.200000.CDS.overlap.density.sorted.coord.txt", header=F)
colnames(CDS) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                 "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

par_rho <- rho[rho$Window_start < 52000000,]
nonpar_rho <- rho[rho$Window_start > 52000000,]

par_pi <- Z_pi[Z_pi$Chr_start < 52000000,]
nonpar_pi <- Z_pi[Z_pi$Chr_start > 52000000,]

par_gc <- GC[GC$Chr_start < 52000000,]
nonpar_gc <- GC[GC$Chr_start > 52000000,]

par_cds <- CDS[CDS$Chr_start < 52000000,]
nonpar_cds <- CDS[CDS$Chr_start > 52000000,]

cds_gc_rho_pi_dataset <- data.frame( cds = c(par_cds$Feat_Base_count/par_cds$Window_Base_count, nonpar_cds$Feat_Base_count/nonpar_cds$Window_Base_count),
                              gc = c(par_gc$Window_GC_count/par_gc$Window_Base_count, nonpar_gc$Window_GC_count/nonpar_gc$Window_Base_count),
                              rho = c(par_rho$rho_per_window, nonpar_rho$rho_per_window),
                              pi = c(par_pi$pi/par_pi$win_length, nonpar_pi$pi/nonpar_pi$win_length),
                              chr = c(rep("PAR", 258), rep("nonPAR", 139)))

rho_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = rho/kb, y = pi, group = chr)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression(4~N[e]/Kb))+
  theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_colour_manual(values=c("azure4", "black")) + ggtitle('D')


#-----------------------
# pi vs. gc
#-----------------------
gc_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = gc*100, y = pi, group = chr)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression("GC %"))+
  theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_colour_manual(values=c("azure4", "black")) + ggtitle('E')

#-----------------------
# pi vs. gene density
#-----------------------
cds_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = cds*100, y = pi, group = chr)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression("CDS %"))+
  theme(legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  scale_colour_manual(values=c("azure4", "black")) + ggtitle('F')

figure2 <- (pi_plot / td_plot / fst_plot) /(rho_pi_plot | gc_pi_plot | cds_pi_plot)
ggsave("../../figures/Figure2.pdf", figure2, width = 6.85, height = 7.85, units = "in")
