#!/usr/bin/R

# Clean the workspace
rm(list = ls())

library(ggplot2)
library(tidyquant)
library(patchwork)
library(nlme)

mb = 10^6
kb = 10^3

#-----------------------
# pi
#-----------------------
Z_pi <- read.table("../../data/diversity/genetic_variation/z/Z.200000.chr.coord.sfs.txt", header=F)
colnames(x = Z_pi) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")

pi_plot <- ggplot(Z_pi, aes(y=pi/win_length, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.text.y=element_text(size=12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  xlab("Position (Mb)") + ylab(expression(pi)) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.0016, linetype="dashed", color = "black") + ggtitle('A')

#-----------------------
# GC
#-----------------------
GC <- read.table("../../data/genomic_features/z_scaf.200000.intergene.overlap.density.sorted.coord.txt", header=F)
colnames(GC) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                 "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

GC_chr4 <- read.table("../../data/genomic_features/chr4.scaf.200000.intergene.overlap.density.sorted.txt", header=T)
GC_chr5 <- read.table("../../data/genomic_features/chr5.scaf.200000.intergene.overlap.density.sorted.txt", header=T)

GC_a <- rbind(GC_chr4, GC_chr5)
GC_a_mean <- mean(GC_a$Window_GC_count/GC_a$Window_Base_count)

gc_plot <- ggplot(GC, aes(y=(GC$Window_GC_count/GC$Window_Base_count)*100, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.text.y=element_text(size=12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  xlab("Position (Mb)") + ylab("GC %") +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=GC_a_mean*100, linetype="dashed", color = "black") + ggtitle('B')

#-----------------------
# CDS
#-----------------------
CDS <- read.table("../../data/genomic_features/z_scaf.200000.CDS.overlap.density.sorted.coord.txt", header=F)
colnames(CDS) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                  "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

CDS_chr4 <- read.table("../../data/genomic_features/chr4.scaf.200000.CDS.overlap.density.sorted.txt", header=T)
CDS_chr5 <- read.table("../../data/genomic_features/chr5.scaf.200000.CDS.overlap.density.sorted.txt", header=T)
CDS_a <- rbind(CDS_chr4, CDS_chr5)

CDS_a_mean <- mean(CDS_a$Feat_Base_count/CDS_a$Window_Base_count)

cds_plot <- ggplot(CDS, aes(y=(CDS$Feat_Base_count/CDS$Window_Base_count)*100, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.text.y=element_text(size=12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  xlab("Position (Mb)") + ylab("CDS %") +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept = CDS_a_mean*100, linetype="dashed", color = "black") + ggtitle('C')

#-----------------------
# TD
#-----------------------
chr4_pi <- read.table("../../data/diversity/genetic_variation/chr4/chr4.200000.sfs.txt", header=T)
chr5_pi <- read.table("../../data/diversity/genetic_variation/chr5/chr5.200000.sfs.txt", header=T)

pi_a <- rbind(chr4_pi, chr5_pi)
td_a_mean <- mean(pi_a$Td)

td_plot <- ggplot(Z_pi, aes(y=Td, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.text.y=element_text(size=12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  xlab("Position (Mb)") + ylab(expression("T"["d"])) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept = td_a_mean, linetype="dashed", color = "black") + ggtitle('D')

#-----------------------
# Fst
#-----------------------
fst <- read.table("../../data/fst/z/black_male_female_Z_200Kb.windowed.weir.Z.coord.fst")
colnames(fst) <- c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start",
                   "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst$WEIGHTED_FST[fst$WEIGHTED_FST<0] <- 0

chr4_fst <- read.table("../../data/fst/autosome/chr4/chr4.200Kb.windowed.weir.fst", header = F)
chr5_fst <- read.table("../../data/fst/autosome/chr5/chr5.200Kb.windowed.weir.fst", header = F)

fst_a <- rbind(chr4_fst, chr5_fst)
colnames(fst_a) <- c("Scaffold", "Window_start", "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst_a$WEIGHTED_FST[fst_a$WEIGHTED_FST<0] <- 0

fst_a_mean <- mean(fst_a$WEIGHTED_FST)

fst_plot <- ggplot(fst, aes(y=WEIGHTED_FST, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.text.y=element_text(size=12), axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  xlab("Position (Mb)") + ylab(expression('F'['FM'])) +  scale_x_continuous(breaks=seq(0,82,10)) +
  geom_vline(xintercept = (52185292)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept = fst_a_mean, linetype="dashed", color = "black") + ggtitle('E')


#-----------------------
# pi vs. rho
#-----------------------
rho <- read.table("../../data/rho/ldhat_rho/z/200Kb200Kb_rho_Z.chr.coord.txt")
colnames(rho) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")

par_rho <- rho[rho$Window_start < 52000000,]
nonpar_rho <- rho[rho$Window_start > 52000000,]

par_pi <- Z_pi[Z_pi$Chr_start < 52000000,]
nonpar_pi <- Z_pi[Z_pi$Chr_start > 52000000,]

par_gc <- GC[GC$Chr_start < 52000000,]
nonpar_gc <- GC[GC$Chr_start > 52000000,]

par_cds <- CDS[CDS$Chr_start < 52000000,]
nonpar_cds <- CDS[CDS$Chr_start > 52000000,]

cds_gc_rho_pi_dataset <- data.frame( cds_data = c(par_cds$Feat_Base_count/par_cds$Window_Base_count),
                              gc_data = c(par_gc$Window_GC_count/par_gc$Window_Base_count),
                              rho_data = c(par_rho$rho_per_window),
                              pi_data = c(par_pi$pi/par_pi$win_length), 
                              chr = c(rep("PAR", 258)))

cds_gc_rho_pi_dataset$obs = as.numeric(rownames(cds_gc_rho_pi_dataset))

fitgls_rho_pi = gls( pi_data ~  +  rho_data, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = cds_gc_rho_pi_dataset )
cds_gc_rho_pi_dataset$predgls_rho_pi = predict(fitgls_rho_pi)

fitgls_gc_pi = gls( pi_data ~  +  gc_data, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = cds_gc_rho_pi_dataset )
cds_gc_rho_pi_dataset$predgls_gc_pi = predict(fitgls_gc_pi)

fitgls_cds_pi = gls( pi_data ~  +  cds_data, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = cds_gc_rho_pi_dataset )
cds_gc_cds_pi_dataset$predgls_cds_pi = predict(fitgls_cds_pi)

rho_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = rho_data/kb, y = pi_data)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression(over(rho,Kb)))+
  theme(axis.text.y=element_text(size=12), legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  scale_colour_manual(values=c("black"))  + ggtitle('F') + geom_line(aes (y = predgls_rho_pi), linewidth =1, col = "grey" )



#-----------------------
# pi vs. gc
#-----------------------
gc_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = gc_data*100, y = pi_data)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression("GC %"))+
  theme(axis.text.y=element_text(size=12), legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  scale_colour_manual(values=c("black")) + ggtitle('G') + geom_line(aes (y = predgls_gc_pi), linewidth =1, col = "grey" )

#-----------------------
# pi vs. gene density
#-----------------------
cds_pi_plot <- ggplot(cds_gc_rho_pi_dataset, aes(x = cds_data*100, y = pi_data)) + ylab(expression(pi)) +
  geom_point(aes(color = chr), size = 0.8, alpha = 0.6) + theme_classic() + xlab(expression("CDS %"))+
  theme(axis.text.y=element_text(size=12), legend.position = "none", axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)) +
  scale_colour_manual(values=c("black")) + ggtitle('H')+ geom_line(aes (y = predgls_cds_pi), linewidth =1, col = "grey" )

figure2 <- (pi_plot / gc_plot / cds_plot / td_plot / fst_plot) /(rho_pi_plot | gc_pi_plot | cds_pi_plot)
ggsave("../../figures/Figure3.pdf", figure2, width = 7.33, height = 10.46, units = "in")
