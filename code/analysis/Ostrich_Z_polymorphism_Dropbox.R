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
library("cowplot")
library(ggpubr)
library(tidyquant)
#*************************************************************************************************
#*Set the working directory
#*************************************************************************************************
setwd("~/Documents/projects/ostrich_Z/results/manuscript_tables")

par(mfrow=c(1,1))

male_boundary_coverage <- read.table("/Users/homapapoli/Documents/projects/ostrich_Z/results/manuscript_tables/male_boundary_coverage.txt", header = F)
female_boundary_coverage <- read.table("/Users/homapapoli/Documents/projects/ostrich_Z/results/manuscript_tables/female_boundary_coverage.txt", header = F)


plot(((max(male_boundary_coverage$V2)-male_boundary_coverage$V2)+3510000)/10^6, rowMeans(male_boundary_coverage[,3:7]), pch = 20, col = "blue", 
     ylab = "Mean Depth", xlab = "Position in Mb (superscaffold36)", bty="l", cex = 0.8)
points(((max(female_boundary_coverage$V2)-female_boundary_coverage$V2)+3510000)/10^6, rowMeans(female_boundary_coverage[,3:7]), pch = 20, col = "red", cex = 0.8)
abline(v = ((max(male_boundary_coverage$V2)-3516672)+3510000)/10^6, lty = 2)
abline(v = ((max(male_boundary_coverage$V2)-3524264)+3510000)/10^6, lty = 2)

legend(3.518, 50, legend=c("Male", "Female"), col=c("blue", "red"), pch = c(20, 20), cex=0.8, box.lty=0)

male_female_het_residuals <- read.table("~/Documents/projects/sex_chr/ostrich_z/male_female_het_residuals.txt", header=T)

par(mfrow=c(1,1))
plot(male_female_het_residuals[male_female_het_residuals$scaf=="superscaffold67",2],
     male_female_het_residuals[male_female_het_residuals$scaf=="superscaffold67",5])
plot(male_female_het_residuals$pos, male_female_het_residuals$residual)

plot(male_female_het_residuals$male_het, male_female_het_residuals$female_het)
hist(male_female_het_residuals$male_het)
hist(male_female_het_residuals$female_het)

hist(male_female_het_residuals$male_het)
hist(male_female_het_residuals[male_female_het_residuals$scaf=="superscaffold36",4])

length(which(male_female_het_residuals$female_het==1))
length(which(male_female_het_residuals$female_het==0.8))
length(which(male_female_het_residuals$female_het==0.6))
length(which(male_female_het_residuals$female_het==0.4))
length(which(male_female_het_residuals$female_het==0.2))
length(which(male_female_het_residuals$female_het==0))

# Gametolog scaffold coverage
setwd("~/Documents/projects/sex_chr/ostrich_z")

male_female_median <- function(male_cov, female_cov){
  m <- read.table(male_cov)
  f <- read.table(female_cov)
  ratio <- median(rowMeans(m[,3:7]))/median(rowMeans(f[,3:7]))
  return(ratio)
}

male_female_median("male_gameto_coverage/scaffold102.male.coverage", "female_gameto_coverage/scaffold102.female.coverage")
male_female_median("male_gameto_coverage/scaffold1051.male.coverage", "female_gameto_coverage/scaffold1051.female.coverage")
male_female_median("male_gameto_coverage/scaffold1066.male.coverage", "female_gameto_coverage/scaffold1066.female.coverage")
male_female_median("male_gameto_coverage/scaffold1051.male.coverage", "female_gameto_coverage/scaffold1101.female.coverage")
male_female_median("male_gameto_coverage/scaffold125.male.coverage", "female_gameto_coverage/scaffold125.female.coverage")
male_female_median("male_gameto_coverage/scaffold178.male.coverage", "female_gameto_coverage/scaffold178.female.coverage")
male_female_median("male_gameto_coverage/scaffold203.male.coverage", "female_gameto_coverage/scaffold203.female.coverage")
male_female_median("male_gameto_coverage/scaffold207.male.coverage", "female_gameto_coverage/scaffold207.female.coverage")
male_female_median("male_gameto_coverage/scaffold240.male.coverage", "female_gameto_coverage/scaffold240.female.coverage")

male_female_median("male_gameto_coverage/scaffold246.male.coverage", "female_gameto_coverage/scaffold246.female.coverage")
male_female_median("male_gameto_coverage/scaffold264.male.coverage", "female_gameto_coverage/scaffold264.female.coverage")
male_female_median("male_gameto_coverage/scaffold283.male.coverage", "female_gameto_coverage/scaffold283.female.coverage")
male_female_median("male_gameto_coverage/scaffold291.male.coverage", "female_gameto_coverage/scaffold291.female.coverage")
male_female_median("male_gameto_coverage/scaffold292.male.coverage", "female_gameto_coverage/scaffold292.female.coverage")
male_female_median("male_gameto_coverage/scaffold314.male.coverage", "female_gameto_coverage/scaffold314.female.coverage")
male_female_median("male_gameto_coverage/scaffold330.male.coverage", "female_gameto_coverage/scaffold330.female.coverage")
male_female_median("male_gameto_coverage/scaffold414.male.coverage", "female_gameto_coverage/scaffold414.female.coverage")
male_female_median("male_gameto_coverage/scaffold432.male.coverage", "female_gameto_coverage/scaffold432.female.coverage")

male_female_median("male_gameto_coverage/scaffold463.male.coverage", "female_gameto_coverage/scaffold463.female.coverage")
male_female_median("male_gameto_coverage/scaffold466.male.coverage", "female_gameto_coverage/scaffold466.female.coverage")
male_female_median("male_gameto_coverage/scaffold502.male.coverage", "female_gameto_coverage/scaffold502.female.coverage")
male_female_median("male_gameto_coverage/scaffold525.male.coverage", "female_gameto_coverage/scaffold525.female.coverage")
male_female_median("male_gameto_coverage/scaffold545.male.coverage", "female_gameto_coverage/scaffold545.female.coverage")
male_female_median("male_gameto_coverage/scaffold589.male.coverage", "female_gameto_coverage/scaffold589.female.coverage")
male_female_median("male_gameto_coverage/scaffold615.male.coverage", "female_gameto_coverage/scaffold615.female.coverage")
male_female_median("male_gameto_coverage/scaffold624.male.coverage", "female_gameto_coverage/scaffold624.female.coverage")
male_female_median("male_gameto_coverage/scaffold67.male.coverage", "female_gameto_coverage/scaffold67.female.coverage")

male_female_median("male_gameto_coverage/scaffold687.male.coverage", "female_gameto_coverage/scaffold687.female.coverage")
male_female_median("male_gameto_coverage/scaffold91.male.coverage", "female_gameto_coverage/scaffold91.female.coverage")
male_female_median("male_gameto_coverage/scaffold916.male.coverage", "female_gameto_coverage/scaffold916.female.coverage")
male_female_median("male_gameto_coverage/superscaffold17.male.coverage", "female_gameto_coverage/superscaffold17.female.coverage")
male_female_median("male_gameto_coverage/superscaffold20.male.coverage", "female_gameto_coverage/superscaffold20.female.coverage")
male_female_median("male_gameto_coverage/superscaffold45.male.coverage", "female_gameto_coverage/superscaffold45.female.coverage")
male_female_median("male_gameto_coverage/superscaffold51.male.coverage", "female_gameto_coverage/superscaffold51.female.coverage")
male_female_median("male_gameto_coverage/superscaffold67.male.coverage", "female_gameto_coverage/superscaffold67.female.coverage")
male_female_median("male_gameto_coverage/superscaffold7.male.coverage", "female_gameto_coverage/superscaffold7.female.coverage")
male_female_median("male_gameto_coverage/superscaffold92.male.coverage", "female_gameto_coverage/superscaffold92.female.coverage")
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
# FIGURE 1
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
# Read gene coordinates
mRNA_coords <- read.table("../../../../sex_chr/ostrich_z/Z_mRNA_chro_coordinates.txt")
gameto_mRNA_coords <- read.table("/Users/homapapoli/Documents/projects/sex_chr/ostrich_z/gametolog_chro_coordinates.txt")

head(mRNA_coords)
plot(c(0, 83), c(0, 3), type= "n", xlab = "", ylab = "",  bty="l", axes=F)
mb = 10^6

for(i in seq(1,849)){
  rect(mRNA_coords[i,5]/mb, 0, mRNA_coords[i,6]/mb, 1, col = "grey", border = "transparent")
}

for(i in seq(1,42)){
  rect(gameto_mRNA_coords[i,5]/mb, 0, gameto_mRNA_coords[i,6]/mb, 1, col = "red", border = "transparent")
}

rect(76174922/mb, 0,  76230182/mb, 1, col = "blue", border = "black")
abline(v = 52193205/mb, col = "black", lty = 2)
rect(0, 0, 82, 1, border = "black")
axis(side = 1, at = seq(0, 82, 10))
text(x = 76174922/mb, y = 1.2, "DMRT1")
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
# FIGURE 2
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
setwd("/Users/homapapoli/Documents/projects/ostrich_Z/results/manuscript_tables/sfs_measures")
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
# SFS
black_SFS <- read.table("black_A.foldedSFS.txt", header = T)
black_par_SFS <- read.table("black_PAR.foldedSFS.txt", header = T)
black_nonpar_SFS <- read.table("black_nonPAR.foldedSFS.txt", header = T)

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

sfs_A <- SFS[SFS$chr=="Autosome",]$sfs
sfs_PAR <- SFS[SFS$chr=="PAR",]$sfs
sfs_expected <- SFS[SFS$chr=="Expected",]$sfs

wilcox.test(sfs_PAR , sfs_A) 
wilcox.test(sfs_expected , sfs_A) 
wilcox.test(sfs_expected , sfs_PAR) 


plotSFS = ggplot(SFS, aes(x = factor(count), y=sfs, fill = chr)) +
  geom_bar(stat = "identity",position=position_dodge()) + ylab("Proportion") + xlab("No. minor alleles") +
  facet_grid(species ~ .) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  strip.text.y = element_blank()) + 
  scale_fill_manual("", values = c("Expected" = "grey", "Autosome" = "blue", "PAR" = "orange")) + theme(legend.position = "none")

# Pi distribution
# For all sequences regardless of functional categories
black_autosome_all <- read.table("black.autosome.sfs.txt", header=T)
black_PAR_all <- read.table("black.PAR.sfs.txt", header=T)
black_nonPAR_all <- read.table("black.nonPAR.sfs.txt", header=T)


black_autosome_all <- black_autosome_all[which(black_autosome_all$win_length>80000),]
black_PAR_all <- black_PAR_all[which(black_PAR_all$win_length>80000),]
black_nonPAR_all <- black_nonPAR_all[which(black_nonPAR_all$win_length>80000),]
# pi
mean(black_autosome_all$pi/black_autosome_all$win_length, na.rm = T)
sd(black_autosome_all$pi/black_autosome_all$win_length, na.rm = T)
mean(black_PAR_all$pi/black_PAR_all$win_length, na.rm = T)
sd(black_PAR_all$pi/black_PAR_all$win_length, na.rm = T)
mean(black_nonPAR_all$pi/black_nonPAR_all$win_length, na.rm = T)
sd(black_nonPAR_all$pi/black_nonPAR_all$win_length, na.rm = T)
# theta
mean(black_autosome_all$theta/black_autosome_all$win_length, na.rm = T)
sd(black_autosome_all$theta/black_autosome_all$win_length, na.rm = T)
mean(black_PAR_all$theta/black_PAR_all$win_length, na.rm = T)
sd(black_PAR_all$theta/black_PAR_all$win_length, na.rm = T)
mean(black_nonPAR_all$theta/black_nonPAR_all$win_length, na.rm = T)
sd(black_nonPAR_all$theta/black_nonPAR_all$win_length, na.rm = T)
# Tajima's D
mean(black_autosome_all$Td, na.rm = T)
sd(black_autosome_all$Td, na.rm = T)
mean(black_PAR_all$Td, na.rm = T)
sd(black_PAR_all$Td, na.rm = T)
mean(black_nonPAR_all$Td, na.rm = T)
sd(black_nonPAR_all$Td, na.rm = T)

t.test(black_autosome_all$Td, black_PAR_all$Td)

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
black_Z <- read.table("black.Z.coordinates.sfs.txt", header=F)
colnames(x = black_Z) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")

mb = 10^6
Pi = ggplot(black_Z, aes(y=pi/black_Z$win_length, x=(black_Z$Chr_start+black_Z$Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Position (Mb)") + ylab(expression(pi)) + xlim(0,81) +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.001642126, linetype="dashed", color = "black") 

gc_density_100Kb <- read.table("../genomic_features/black.Z.genedensity.500Kb.txt")
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

working_dir = "~/Documents/projects/ostrich_Z/results/manuscript_tables"
mb = 10^6
female_map <- read.table(paste(working_dir, "/linkage_map/LGZ3.female.lifted.cleaned.bed.txt", sep = ""))
colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
male_map <- read.table(paste(working_dir, "/linkage_map/LGZ3.male.lifted.cleaned.bed.txt", sep = ""))
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

#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*Figure 3
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
black_LD <- read.table("../LD/LD_chromosome_plot/black.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out")
colnames(x = black_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                           "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                           "mean_SNP_distance", "NUM_snp_pairs")
# LD decay
scaf26 <- read.table("../LD/scaf.split/black.superscaffold26.LDdecay.bin")
scaf54<- read.table("../LD/scaf.split/black.superscaffold54.LDdecay.bin")
scaf35<- read.table("../LD/scaf.split/black.superscaffold35.LDdecay.bin")
scaf36<- read.table("../LD/scaf.split/black.superscaffold36.LDdecay.bin")
scaf62<- read.table("../LD/scaf.split/black.superscaffold62.LDdecay.bin")
scaf67<- read.table("../LD/scaf.split/black.superscaffold67.LDdecay.bin")
scaf69_1<- read.table("../LD/scaf.split/black.superscaffold69-1.LDdecay.bin")
scaf93<- read.table("../LD/scaf.split/black.superscaffold93.LDdecay.bin")
scaf63<- read.table("../LD/scaf.split/black.superscaffold63.LDdecay.bin")
scaf88<- read.table("../LD/scaf.split/black.superscaffold88.LDdecay.bin")
scaf83<- read.table("../LD/scaf.split/black.superscaffold83.LDdecay.bin")
scaf92<- read.table("../LD/scaf.split/black.superscaffold92.LDdecay.bin")

# LD for PAR and non-PAR
black_LDPAR <- black_LD[black_LD$Window_start < 53065700,]
black_LDnonPAR <- black_LD[black_LD$Window_start > 53065700,]


library(tidyquant)
library(LDheatmap)
mb = 10^6
plotld = ggplot(black_LD, aes(y=r_squared, x=(Window_start+Window_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  geom_vline(xintercept = 52193205/mb, linetype="dashed", color = "black") + 
  xlab("Position (Mb)") + ylab(expression("r"^{2})) 


LD <-as.matrix(read.table("black.superscaffold36.bothsexes.100Kb.boundary.ld"))
snpinfo <- read.table("boundary.z.coordinates.txt", header=F)
LDheat <- LDheatmap(LD, genetic.distances = snpinfo$V5, distances = "physical", color=heat.colors(20), title = NULL, flip = TRUE)

# LD decay at the boundary
head(snpinfo$V4-snpinfo$V3)
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*Figure 4
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
td = ggplot(black_Z, aes(y=Td, x=(Chr_start+Chr_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  xlab("Position (Mb)") + ylab("Tajima's D") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_vline(xintercept = (52193205)/mb, linetype="dashed", color = "black") 

mfFST <- read.table("black_male_female_Z_100Kb.windowed.weir.Z.coord.fst", header=F)
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

#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
#*Figure 5
#*************************************************************************************************
#*************************************************************************************************
#*************************************************************************************************
positions <- cum_positions[1:(length(cum_positions )-1)]/1000000
rho_rates_kb <- cum_rates[2:length(cum_rates)]/1000

pos_rho <- data.frame(positions, rho_rates_kb)
write.table(pos_rho, "~/Documents/projects/ostrich_Z/results/manuscript_tables/pos_rho.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

positions <- sex_averaged_map$Pos/mb
sex_averaged_cM <- sex_averaged_map$sex_averaged_cM
sex_averaged_smooth <- sex_average_map_smoothed25

pos_cm_smooth <- data.frame(positions, sex_averaged_cM, sex_averaged_smooth)
write.table(pos_cm_smooth, "~/Documents/projects/ostrich_Z/results/manuscript_tables/pos_cm_smooth.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


plot(cum_positions[1:(length(cum_positions )-1)]/1000000, y=cum_rates[2:length(cum_rates)]/1000, type="s", col=rgb(0,0,0.5), 
     xlab = "Position (Mb)", ylab = "4Ner/Kb", las = 1)

par(new=TRUE)

plot(sex_averaged_map$Pos/mb, sex_averaged_map$sex_averaged_cM, col = "darkorchid1", xlab="", ylab="", ylim = c(0, 95), pch = 20, axes = FALSE)
lines(sex_averaged_map$Pos/mb, sex_average_map_smoothed25, col = "darkorchid1")

mtext("cM",side=4,col="red",line=4) 
axis(4, ylim=c(0,95),las=1)

legend("topleft",legend=c("4Ner","cM"),
       text.col=c("black","darkorchid1"),pch=c(20),lty=1:2,col=c("black","darkorchid1"), box.lty=0)


rho <- read.table("~/Documents/projects/ostrich_Z/results/manuscript_tables/pos_rho.txt", header=T)
cm <- read.table("~/Documents/projects/ostrich_Z/results/manuscript_tables/pos_cm_smooth.txt", header=T)

plot(rho$positions, rho$rho_rates_kb, type="s", col=rgb(0,0,0.5), xlab = "Position (Mb)", ylab = "4Ner/Kb", las = 1)
par(new=TRUE)
plot(cm$positions, cm$sex_averaged_cM, col = "darkorchid1", xlab="", ylab="", ylim = c(0, 65), pch = 20, axes = FALSE)
lines(cm$positions, cm$sex_averaged_smooth, col = "darkorchid1")
mtext("cM",side=4,col="red",line=4) 
axis(4, ylim=c(0,95),las=1)

legend("topleft",legend=c("4Ner","cM"),
       text.col=c("black","darkorchid1"),pch=c(20),lty=1:2,col=c("black","darkorchid1"), box.lty=0)


#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
#*Measures of site frequency spectrum
#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
# Functional categories diversity values (Table 1)
# Intergenic
black_A_intergenic <- read.table("black.A.intergenic.sfs.txt", header = T)
mean(black_A_intergenic$pi/black_A_intergenic$win_length, na.rm = T)
black_PAR_intergenic <- read.table("black.PAR.intergenic.sfs.txt", header = T)
mean(black_PAR_intergenic$pi/black_PAR_intergenic$win_length, na.rm = T)
black_nonPAR_intergenic <- read.table("black.nonPAR.intergenic.sfs.txt", header = T)
mean(black_nonPAR_intergenic$pi/black_nonPAR_intergenic$win_length, na.rm = T)
# Intronic
black_A_intron <- read.table("black.A.intronic.sfs.txt", header = T)
mean(black_A_intron$pi/black_A_intron$win_length, na.rm = T)
black_PAR_intron <- read.table("black.PAR.intronic.sfs.txt", header = T)
mean(black_PAR_intron$pi/black_PAR_intron$win_length, na.rm = T)
black_nonPAR_intron <- read.table("black.nonPAR.intronic.sfs.txt", header = T)
mean(black_nonPAR_intron$pi/black_nonPAR_intron$win_length, na.rm = T)
# CDS: 4fold and 0fold
black_A_pnps <- read.table("black.A.pnps.txt", header = T)
mean(black_A_pnps$ps, na.rm = T)
mean(black_A_pnps$pn, na.rm = T)
mean(black_A_pnps$pnps, na.rm = T)
black_PAR_pnps <- read.table("black.PAR.pnps.txt", header = T)
mean(black_PAR_pnps$ps, na.rm = T)
mean(black_PAR_pnps$pn, na.rm = T)
mean(black_PAR_pnps$pnps, na.rm = T)
black_nonPAR_pnps <- read.table("black.nonPAR.pnps.txt", header = T)
mean(black_nonPAR_pnps$ps, na.rm = T)
mean(black_PAR_pnps$pn, na.rm = T)
mean(black_PAR_pnps$pnps, na.rm = T)

# Z coordinate
black_Z <- read.table("black.sfs.Z.txt", header =T)
black_Z_intergenic <- read.table("black.Z.intergenic.sfs.txt", header = T)

# SFS measures 0.5 and 1 MB
black_PAR_all_0.5Mb_PAR <- read.table("black.PAR.500Kb.sfs.txt", header=T)
black_nonPAR_all_0.5Mb_PAR <- read.table("black.nonPAR.500Kb.sfs.txt", header=T)

black_PAR_all_1Mb_PAR <- read.table("black.PAR.1000Kb.sfs.txt", header=T)
black_nonPAR_all_1Mb_PAR <- read.table("black.nonPAR.1000Kb.sfs.txt", header=T)

# Genomic features
gene_GC_0.1Mb_PAR <- read.table("../genomic_features/black.par_scaf.100Kb.CDS.overlap.density.sorted.txt")
gene_GC_0.5Mb_PAR <- read.table("../genomic_features/black.par_scaf.500Kb.CDS.overlap.density.sorted.txt")
gene_GC_1Mb_PAR <- read.table("../genomic_features/black.par_scaf.1000Kb.CDS.overlap.density.sorted.txt")
colnames(x = gene_GC_0.1Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                    "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_0.5Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                    "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_1Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                  "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
gene_GC_0.1Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.100Kb.CDS.overlap.density.sorted.txt")
gene_GC_0.5Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.500Kb.CDS.overlap.density.sorted.txt")
gene_GC_1Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.1000Kb.CDS.overlap.density.sorted.txt")
colnames(x = gene_GC_0.1Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                       "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_0.5Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                       "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_1Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                     "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

GC_PAR <- read.table("../genomic_features/black.Z.density.txt", header=T)
gc_density_100Kb <- read.table("../genomic_features/black.Z.genedensity.100Kb.txt")

gc_density_1000Kb <- read.table("../genomic_features/black.Z.genedensity.1000Kb.txt")

gc_density_500Kb <- read.table("../genomic_features/black.Z.genedensity.500Kb.txt")
colnames(x = gc_density_500Kb) = c("Chr", "Scaffold", "Window_start", "Window_end", "cstart", "cend", "Window_Base_count", 
                                   "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

Z_500Kb <- read.table("black.Z.500Kb.sfs.txt", header =T)

# Recombination rate in 0.5 Mb windows
start_chr <- sex_averaged_rate[[2]]$start
end_chr <- sex_averaged_rate[[2]]$end
smooth_rate <- sex_averaged_rate[[2]]$pair_cm_per_site
df1 <- data.frame(start_chr, end_chr, smooth_rate, smooth_rate/mb)

sex_averaged_map_1 <- rbind(c(0,0,0,0),sex_averaged_map, c(80543761, 80926604, NA, NA))
sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map_1, len_markers = 195, column = 4, span = 0.3)
test <-c(sex_averaged_rate[[4]],0,0,0,0,0,0,0,0)

#*************************************************************************************************
#*************************************************************************************************
# GLS
#*************************************************************************************************
#*************************************************************************************************
# PAR section 1
plot((gc_density_500Kb$Window_GC_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
cor.test((gc_density_500Kb$Window_GC_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
plot((gc_density_500Kb$Feat_Base_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
cor.test((gc_density_500Kb$Feat_Base_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]),  Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])

# PAR section 2
plot((gc_density_500Kb$Window_GC_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
cor.test((gc_density_500Kb$Window_GC_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
plot((gc_density_500Kb$Feat_Base_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
cor.test((gc_density_500Kb$Feat_Base_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])

# nonPAR
cor.test(gc_density_500Kb$Window_GC_count[106:167]/gc_density_500Kb$Window_Base_count[106:167], Z_500Kb$pi[106:167]/Z_500Kb$winsize[106:167])
cor.test(gc_density_500Kb$Feat_Base_count[106:167]/gc_density_500Kb$Window_Base_count[106:167], Z_500Kb$pi[106:167]/Z_500Kb$winsize[106:167])

# GLS
library(nlme)
scaled_Z <- scale(Z_500Kb$pi/Z_500Kb$winsize)
scaled_gc <- scale(gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count)
scaled_gene <- scale(gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count)
d <- data.frame(scaled_Z, scaled_gc, scaled_gene, scaled_rec)

cor_strut <- pcr_df 
d$obs = as.numeric(rownames(d))
full_PAR <- gls( scaled_Z ~  +  scaled_gc + scaled_gene + scaled_rec, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = d )

summary(full_PAR)

scaled_rec <- scale(rec)
diversity <- Z_500Kb$pi/Z_500Kb$winsize
gc <- gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count
gene_dense <- gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count
rec <- test

cor.test(diversity, gc)
cor.test(diversity, gene_dense)
cor.test(diversity, rec)

cor.test(gc, rec)

pcr_df <- data.frame(diversity, gc, gene_dense, rec)

library(pls)
set.seed(1)
model <- pcr(diversity~gc+gene_dense+rec, data=pcr_df, scale=TRUE, validation="CV")

summary(model)

validationplot(model, val.type="R2", cex.axis=0.7)
axis(side = 1, at = c(8), cex.axis=0.7)
abline(v = 8, col = "blue", lty = 3)

validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

# The model improves by adding three PCAs

train <- pcr_df[1:130, c("diversity", "gc", "gene_dense", "rec")]
y_test <- pcr_df[131:nrow(pcr_df), c("diversity")]
test <- pcr_df[131:nrow(pcr_df), c("gc", "gene_dense", "rec")]
model <- pcr(diversity~gc+gene_dense+rec, data=train, scale=TRUE, validation="CV")
pcr_pred <- predict(model, test, ncomp=3)
sqrt(mean((pcr_pred - y_test)^2))

library(factoextra)
d <- prcomp(~Z_500Kb$pi/Z_500Kb$winsize+gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count+
              gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count+test, scale = TRUE)
fviz_eig(d)
fviz_pca_ind(d,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#*************************************************************************************************
#* rho estimation
#*************************************************************************************************
r_Ne_table <- read.table("../rho_r_Ne_table.txt", header = T)
head(r_Ne_table)
attach(r_Ne_table)
cor.test(rho, smooth_rate_map_cM_Mb)


plot((r_Ne_table$start+r_Ne_table$end)/(2*mb), Ne_smooth_r, pch = 20)

range(r_Ne_table$Ne_smooth_r)
range(r_Ne_table$Ne_real_r, na.rm = T)

#*************************************************************************************************
#*************************************************************************************************
# Linkage Disequilibrium
#***********************************************************************************************
mean(black_LD$r_squared, na.rm = T)
sd(black_LD$r_squared, na.rm = T)

mean(black_LDPAR$r_squared)
sd(black_LDPAR$r_squared)

mean(black_LDnonPAR$r_squared, na.rm = T)
sd(black_LDnonPAR$r_squared, na.rm = T)

t.test(black_LDPAR$r_squared, black_LDnonPAR$r_squared)
#***********************************************************************************************
# LD decay
#***********************************************************************************************
plot(scaf26$V1, scaf26$V2)
PAR_ld <- c(scaf26$V2, scaf54$V2, scaf35$V2)
nonPAR_ld <- c(scaf62$V2, scaf63$V2, scaf67$V2, scaf69_1$V2, scaf83$V2, scaf88$V2, scaf92$V2, scaf93$V2)
# mean of LD decay PAR
mean(PAR_ld)
sd(PAR_ld)
# mean of LD decay nonPAR
mean(nonPAR_ld)
sd(nonPAR_ld)
#***********************************************************************************************
#***********************************************************************************************
# Dos and HKA
dnds <- read.table("../divergence/ostrich_dnds_pariwise_txt", header = T)
head(dnds)


boxplot(black_A_pnps$ps, black_PAR_pnps$ps, black_nonPAR_pnps$ps, outline = F, ylab = "Ps", names = c("Autosome Ps", "PAR Ps", "nonPAR Ps"))
boxplot(black_A_pnps$pn, black_PAR_pnps$pn, black_nonPAR_pnps$pn, outline = F, ylab = "Pn", names = c("Autosome Pn", "PAR Pn", "nonPAR Pn"))
boxplot(black_A_pnps$pnps, black_PAR_pnps$pnps, black_nonPAR_pnps$pnps, outline = F, ylab = "PnPs", names = c("Autosome PnPs", "PAR PnPs", "nonPAR PnPs"))
median(black_A_pnps$ps)
median(black_PAR_pnps$ps)
median(black_nonPAR_pnps$ps)

median(black_A_pnps$pn)
median(black_PAR_pnps$pn)
median(black_nonPAR_pnps$pn)

median(black_A_pnps$pnps)
median(black_PAR_pnps$pnps)
median(black_nonPAR_pnps$pnps)

# dnds
dnds <- read.table("../divergence/ostrich_dnds_pariwise_txt", header = T)
head(dnds)
PAR_scaffolds = c("superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36")
nonPAR_scaffolds = c("superscaffold69-1", "superscaffold67", "superscaffold62", "superscaffold63", "superscaffold88", "superscaffold93")
dnds_PAR <- dnds[dnds$Scaffold %in% PAR_scaffolds,]
dnds_nonPAR <- dnds[dnds$Scaffold %in% nonPAR_scaffolds,]

boxplot(dnds_PAR$dS, dnds_nonPAR$dS, outline = F, ylab = "ds", names = c("PAR ds", "nonPAR ds"))
boxplot(dnds_PAR$dN, dnds_nonPAR$dN, outline = F, ylab = "dn", names = c("PAR dn", "nonPAR dn"))
boxplot(black_A_pnps$pnps, black_PAR_pnps$pnps, black_nonPAR_pnps$pnps, dnds_PAR$dNdS, dnds_nonPAR$dNdS, outline = F, ylab = "dnds", names = c("PAR dnds", "nonPAR dnds"))

median(dnds_PAR$dS)
median(dnds_nonPAR$dS)
median(dnds_PAR$dN)
median(dnds_nonPAR$dN)
median(dnds_PAR$dNdS)
median(dnds_nonPAR$dNdS)


DoS = Dn/(Dn + Ds) - Pn/(Pn + Ps)

PAR <- read.table("black.PAR.pnps.number_syn_numer_nonsyn.txt", header = T)
nonPAR <- read.table("black.nonPAR.pnps.number_syn_numer_nonsyn.txt", header = T)
PAR[,2] <- sub("_XP_rna", "", PAR[,2])
dnds[,2] <- sub("XM_", "", dnds[,2])
nonPAR[,2] <- sub("_XP_rna", "", nonPAR[,2])

pnps_dnds_list <- dnds[dnds$geneID %in% PAR$gene,]
pnps_pnps_list <- PAR[PAR$gene %in% pnps_dnds_list$geneID,]

nonPAR_pnps_dnds_list <- dnds[dnds$geneID %in% nonPAR$gene,]
nonPAR_pnps_pnps_list <- nonPAR[nonPAR$gene %in% nonPAR_pnps_dnds_list$geneID,]

DoS <- function(Dn, Ds, Pn, Ps){
  DoS_cal = Dn/(Dn + Ds) - Pn/(Pn + Ps)
  return(DoS_cal)
}

alpha <- function(Dn, Ds, Pn, Ps){
  alpha_cal <- 1 - (Ds*Pn)/(Dn*Ps) 
  return(alpha_cal)
}

DoS_PAR <- c()
c = 0
Pn_t = 0
Ps_t = 0
Dn_t = 0
Ds_t = 0
positive_DoS = c()
for(gene in pnps_dnds_list$geneID){
  c = c + 1
  Pn = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NN
  Ps = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NS
  N = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$N
  S = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$S
  dN = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dN
  dS = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dS
  Dn = N * dN
  Ds = S * dS
  
  DoS_PAR[c] = DoS(Dn, Ds, Pn, Ps)
  if(DoS(Dn, Ds, Pn, Ps) > 0){
    print(c(gene, DoS(Dn, Ds, Pn, Ps)))
    positive_DoS[c] = gene
  }
  Pn_t = Pn_t + Pn
  Ps_t = Ps_t + Ps
  Dn_t = Dn_t + Dn
  Ds_t = Ds_t + Ds
  #print(c(Pn_t, Ps_t, Dn_t, Ds_t))
}

c = 0
gene_id = c()
L = c()
Seg = c()
n = c()
D = c()
startingtheta = c()
inheritance = c()
for(gene in pnps_dnds_list$geneID){
  c = c + 1
  gene_id[c] = gene
  L[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$whole_syn
  Seg[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NS
  n[c] = 20
  S = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$S
  dS = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dS
  D[c] = S * dS
  startingtheta[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$ps
  inheritance[c] = 1
}

HKA_dataframe <- data.frame(gene_id, L, Seg, n, D, startingtheta, inheritance)
write.table(HKA_dataframe, file="HKA_dataframe", sep = "\t", quote=FALSE, row.names = FALSE)


PAR[PAR$gene %in% positive_DoS,]

alpha <- function(Dn, Ds, Pn, Ps){
  alpha_cal <- 1 - (Ds*Pn)/(Dn*Ps) 
  return(alpha_cal)
}

alpha(4508.9460, 7627.0692, 94.5464, 60.6771)


alpha_PAR[c] = alpha(Dn, Ds, Pn, Ps)

hist(DoS_PAR, breaks = 10)

median(DoS_PAR)

length(which(DoS_PAR<0))/length(DoS_PAR)

DoS_nonPAR <- c()
c = 0
for(gene in nonPAR_pnps_dnds_list$geneID){
  c = c + 1
  Pn = nonPAR_pnps_pnps_list[nonPAR_pnps_pnps_list$gene==gene,]$NN
  Ps = nonPAR_pnps_pnps_list[nonPAR_pnps_pnps_list$gene==gene,]$NS
  N = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$N
  S = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$S
  dN = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$dN
  dS = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$dS
  Dn = N * dN
  Ds = S * dS
  print(c(Pn, Ps, Dn, Ds))
  print(DoS(Dn, Ds, Pn, Ps))
  DoS_nonPAR[c] = DoS(Dn, Ds, Pn, Ps)
}

median(DoS_nonPAR)

length(which(DoS_nonPAR<0))/length(DoS_nonPAR)

hist(DoS_nonPAR, breaks = 10)


#*************************************************************************************************
#*************************************************************************************************
# PSMC Ne
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Black_P1878_110_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Blue_P1878_121_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Red_P1878_128_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)


#********************************************************************************************************************************************
#********************************************************************************************************************************************
#********************************************************************************************************************************************
#********************************************************************************************************************************************











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
# SFS
black_SFS_4_A <- read.table("cds_sfs_measures/fourfold_A_SFS.txt", header = T)
black_SFS_0_A <- read.table("cds_sfs_measures/zerofold_A_SFS.txt", header = T)

black_SFS_4_PAR <- read.table("cds_sfs_measures/fourfold_PAR_SFS.txt", header = T)
black_SFS_0_PAR <- read.table("cds_sfs_measures/zerofold_PAR_SFS.txt", header = T)

black_SFS_4_nonPAR <- read.table("cds_sfs_measures/fourfold_nonPAR_SFS.txt", header = T)
black_SFS_0_nonPAR <- read.table("cds_sfs_measures/zerofold_nonPAR_SFS.txt", header = T)


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
  chr = c(rep("nonsyn", 10), rep("syn", 10), rep("Expected", 10)),
  sfs = c(black_SFS_0_A$chr_norm, black_SFS_4_A$chr_norm, folded_20/sum(folded_20)),
  count = c(black_SFS_0_A$count, black_SFS_4_A$count, black_SFS_4_A$count)
)

SFS <- within(SFS, chr <- factor(chr, levels=c("nonsyn", "syn", "Expected")))

SFS_PAR <- data.frame(
  species = rep("Black", 30),
  chr = c(rep("nonsyn", 10), rep("syn", 10), rep("Expected", 10)),
  sfs = c(black_SFS_0_PAR$chr_norm, black_SFS_4_PAR$chr_norm, folded_20/sum(folded_20)),
  count = c(black_SFS_0_PAR$count, black_SFS_4_PAR$count, black_SFS_4_PAR$count)
)

SFS_PAR <- within(SFS_PAR, chr <- factor(chr, levels=c("nonsyn", "syn", "Expected")))


plotSFS = ggplot(SFS, aes(x = factor(count), y=sfs, fill = chr)) +
  geom_bar(stat = "identity",position=position_dodge()) + ylab("Proportion") + xlab("No. minor alleles") +
  facet_grid(species ~ .) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  strip.text.y = element_blank()) + 
  scale_fill_manual("", values = c("nonsyn" = "red", "syn" = "blue", "Expected" = "grey")) + theme(legend.position = "none")

plotSFS_PAR = ggplot(SFS_PAR, aes(x = factor(count), y=sfs, fill = chr)) +
  geom_bar(stat = "identity",position=position_dodge()) + ylab("Proportion") + xlab("No. minor alleles") +
  facet_grid(species ~ .) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                  strip.text.y = element_blank()) + 
  scale_fill_manual("", values = c("nonsyn" = "red", "syn" = "blue", "Expected" = "grey")) + theme(legend.position = "none")

start_rev	end_rev	real_rate_map_cM_Mb_rev	smooth_rate_map_cM_Mb_rev	rho_rev	Ne_smooth_r_rev	Ne_real_r_rev
0	1e+06	1.59632630373809	1.072339914362	538.280208596404	12549.1973530767	8429.98400978426

chrom   chromStart      chromEnd
superscaffold26 0       25310599
superscaffold54 0       16379243
superscaffold35 0       4625539
superscaffold36 3524263 9394175
