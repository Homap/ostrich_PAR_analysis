#------------------------------------------------------------
# Homa Papoli
# Created on: May 2022
# Description: Script to create plots of Tajima's D distributions
# for Italian sparrow population based on thetas Analysis from ANGSD
# Note: Before running, set the working directory to the root of the project
#------------------------------------------------------------
rm(list = ls())
library(tidyverse)
# SFS
black_SFS <- read.table("data/diversity/black_A.foldedSFS.txt", header = T)
black_par_SFS <- read.table("data/diversity/black_PAR.foldedSFS.txt", header = T)
black_nonpar_SFS <- read.table("data/diversity/black_nonPAR.foldedSFS.txt", header = T)

# SFS
folded_20 = c()
i = seq(1, 10)
n = 20
for( h in i){
  if(h == 10){
    folded_20[h] =((1/h)+(1/(n-h)))/(1+1)
  }else{
    folded_20[h] = ((1/h)+(1/(n-h)))
  }
}

SFS <- data.frame(
  species = rep("Black", 30),
  chr = c(rep("Expected", 10), rep("Autosome", 10), rep("PAR", 10)),
  sfs = c(folded_20/sum(folded_20), black_SFS$chr_norm, black_par_SFS$chr_norm),
  count = c(black_SFS$count, black_SFS$count, black_par_SFS$count)
)

SFS <- within(SFS, chr <- factor(chr, levels=c("Expected", "Autosome", "PAR")))

plotSFS = ggplot(SFS, aes(x = factor(count), y=sfs, fill = chr)) +
  geom_bar(stat = "identity",position=position_dodge()) + ylab("Proportion") + xlab("No. minor alleles") +
  facet_grid(species ~ .) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  strip.text.y = element_blank()) + 
  scale_fill_manual("", values = c("Expected" = "grey", "Autosome" = "blue", "PAR" = "orange")) + theme(legend.position = "none")

# Pi distribution
# For all sequences regardless of functional categories
black_autosome_all <- read.table("data/diversity/black.autosome.sfs.txt", header=T)
black_PAR_all <- read.table("data/diversity/black.PAR.sfs.txt", header=T)
black_nonPAR_all <- read.table("data/diversity/black.nonPAR.sfs.txt", header=T)


black_autosome_all <- black_autosome_all[which(black_autosome_all$win_length>80000),]
black_PAR_all <- black_PAR_all[which(black_PAR_all$win_length>80000),]
black_nonPAR_all <- black_nonPAR_all[which(black_nonPAR_all$win_length>80000),]

pi_dist = data.frame(
  Pi=c(black_autosome_all$pi[1:12451]/black_autosome_all$win_length[1:12451],black_PAR_all$pi/black_PAR_all$win_length, black_nonPAR_all$pi/black_nonPAR_all$win_length),
  Species=rep('black', 13210),
  Chromosome=c(rep("Autosome", 12451), rep("PAR", 534), rep("nonPAR", 225))
)

pdat_pi <- pi_dist %>%
  group_by(Species, Chromosome) %>%
  do(data.frame(loc = density(.$Pi, na.rm = T)$x,
                dens = density(.$Pi, na.rm = T)$y))

pdat_pi$dens <- ifelse(pdat_pi$Chromosome == 'Autosome', pdat_pi$dens * -1, pdat_pi$dens)

plotPi = ggplot(pdat_pi, aes(x = loc, y = dens, fill = Chromosome, group = interaction(Chromosome, Species))) + 
  geom_polygon(colour = "black") + 
  xlab(expression(pi)) + xlim(c(0, 0.004)) + ylab ("Density") + ylim(c(0, 2000)) +
  theme_classic()  + scale_fill_manual(values = alpha(c("blue", "green", "orange"), .6)) + theme(legend.position = "none")

# Z chromosome coordinates
black_Z <- read.table("data/diversity/black.Z.coordinates.sfs.txt", header=F)
colnames(x = black_Z) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")

mb = 10^6
Pi = ggplot(black_Z, aes(y=pi/black_Z$win_length, x=(black_Z$Chr_start+black_Z$Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Position (Mb)") + ylab(expression(pi)) + xlim(0,81) +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.001642126, linetype="dashed", color = "black") 

gc_density_500Kb <- read.table("data/gc_genedensity/black.Z.genedensity.500Kb.txt")
colnames(x = gc_density_100Kb) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_Base_count", 
                                   "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

GC = ggplot(gc_density_100Kb, aes(y=100*(gc_density_100Kb$Window_GC_count/gc_density_100Kb$Window_Base_count), x=((gc_density_100Kb$Chr_start)+(gc_density_100Kb$Chr_end))/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7)+ geom_ma(ma_fun = SMA, n =2, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Position (Mb)") + ylab("GC%") +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") + xlim(0,81)

gene_density = ggplot(gc_density_100Kb, aes(y=gc_density_100Kb$Feat_Base_count/gc_density_100Kb$Window_Base_count, x=((gc_density_100Kb$Chr_start)+(gc_density_100Kb$Chr_end))/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7)+ geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Position (Mb)") + ylab("CDS density") +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") +xlim(0,81)


mb = 10^6
female_map <- read.table("data/geneticmap/LGZ3.female.lifted.cleaned.bed.txt", sep = "")
colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
male_map <- read.table("data/geneticmap/LGZ3.male.lifted.cleaned.bed.txt", sep = "")
colnames(male_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")

sex_averaged_map <- inner_join(female_map, male_map, by = "Pos")
sex_averaged_cM <- (sex_averaged_map$cM.x + sex_averaged_map$cM.y)/2
Pos <- sex_averaged_map$Pos
female_cM <- sex_averaged_map$cM.x
male_cM <- sex_averaged_map$cM.y
sex_averaged_map <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM)
head(sex_averaged_map)

female_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 2, span = 0.3)
male_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 3, span = 0.3)
sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 4, span = 0.3)

rec <- data.frame(
  position=c((female_rate[[2]]$start+female_rate[[2]]$end)/(2*mb), 
             (male_rate[[2]]$start+male_rate[[2]]$end)/(2*mb),
             (sex_averaged_rate[[2]]$start+sex_averaged_rate[[2]]$end)/(2*mb)),
  rate=c(female_rate[[2]]$pair_cm_per_site, male_rate[[2]]$pair_cm_per_site, sex_averaged_rate[[2]]$pair_cm_per_site),
  sex=c(rep('female', 78), rep('male', 78), rep('sex_averaged', 78))
)

rec_rate <- ggplot(rec , aes(x = position, y = rate, group = sex)) + 
  geom_point(aes(colour=sex), size = 1.3, alpha = 0.7)+ geom_ma(aes(color = sex), ma_fun = SMA, n = 5, size = 0.8, linetype = 6) + theme_classic()+
  xlab("Position (Mb)") + ylab("cM/Mb") +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") + 
  scale_colour_manual(values=c("red", "blue","purple")) + theme(legend.position = "none") +xlim(0,81)


ggdraw(ylim = c(0,2.6)) +
  draw_plot(plotSFS, x = 0.02, y = 2.1, width = .48, height = .5) +
  draw_plot(plotPi, x = .5, y = 2.11, width = .48, height = .49) +
  draw_plot(Pi, x = 0, y = 1.6, width = 1, height = 0.5) +
  draw_plot(GC, x = 0.02, y = 1.1, width = 0.98, height = 0.5) +
  draw_plot(gene_density, x = 0.005, y = 0.6, width = 0.99, height = 0.5) +
  draw_plot(rec_rate, x = 0.03, y = 0, width = 0.97, height = 0.6) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 12,
                  x = c(0, 0.5, 0, 0, 0, 0), y = c(2.6, 2.6, 2.1,1.6,1.1,0.6))