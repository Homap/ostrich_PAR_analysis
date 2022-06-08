#------------------------------------------------------------
# Homa Papoli
# Created on: May 2022
# Description: Script to create plots of Tajima's D distributions
# for Italian sparrow population based on thetas Analysis from ANGSD
# Note: Before running, set the working directory to the root of the project
#------------------------------------------------------------
rm(list = ls())
library(tidyquant)
library(LDheatmap)
library(ggplot2)
#*************************************************************************************************
#*Chromosome plot
#*************************************************************************************************
# Read gene coordinates
mRNA_coords <- read.table("data/gene_gameto_coordinates/Z_mRNA_chro_coordinates.txt")
gameto_mRNA_coords <- read.table("data/gene_gameto_coordinates/gametolog_chro_coordinates.txt")

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
#*LD (r^2) across the chromosome Z
#*************************************************************************************************
black_LD <- read.table("data/ld/window_ld/black.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out")
colnames(x = black_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                           "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                           "mean_SNP_distance", "NUM_snp_pairs")


mb = 10^6
plotld = ggplot(black_LD, aes(y=r_squared, x=(Window_start+Window_end)/(2*mb))) +
  geom_point(size = 1.3, col="grey", alpha = 0.7) + geom_ma(ma_fun = SMA, n = 5, col = "black", size = 0.8, linetype = 6) + theme_classic()+
  geom_vline(xintercept = 52193205/mb, linetype="dashed", color = "black") + 
  xlab("Position (Mb)") + ylab(expression("r"^{2})) 


LD <-as.matrix(read.table("data/ld/window_ld/boundary/black.superscaffold36.bothsexes.100Kb.boundary.ld"))
snpinfo <- read.table("data/ld/window_ld/boundary/boundary.z.coordinates.txt", header=F)
LDheat <- LDheatmap(LD, genetic.distances = snpinfo$V5, distances = "physical", color=heat.colors(20), title = NULL, flip = TRUE)


#*************************************************************************************************
#*LD decay
#*************************************************************************************************
# LD decay
scaf26 <- read.table("data/ld/ld_decay/black.superscaffold26.LDdecay.bin")
scaf54<- read.table("data/ld/ld_decay/black.superscaffold54.LDdecay.bin")
scaf35<- read.table("data/ld/ld_decay/black.superscaffold35.LDdecay.bin")
scaf36<- read.table("data/ld/ld_decay/black.superscaffold36.LDdecay.bin")
scaf62<- read.table("data/ld/ld_decay/black.superscaffold62.LDdecay.bin")
scaf67<- read.table("data/ld/ld_decay/black.superscaffold67.LDdecay.bin")
scaf69_1<- read.table("data/ld/ld_decay/black.superscaffold69-1.LDdecay.bin")
scaf93<- read.table("data/ld/ld_decay/black.superscaffold93.LDdecay.bin")
scaf63<- read.table("data/ld/ld_decay/black.superscaffold63.LDdecay.bin")
scaf88<- read.table("data/ld/ld_decay/black.superscaffold88.LDdecay.bin")
scaf83<- read.table("data/ld/ld_decay/black.superscaffold83.LDdecay.bin")
scaf92<- read.table("data/ld/ld_decay/black.superscaffold92.LDdecay.bin")


# LD for PAR and non-PAR
black_LDPAR <- black_LD[black_LD$Window_start < 53065700,]
black_LDnonPAR <- black_LD[black_LD$Window_start > 53065700,]

# LD decay
d <- read.table("~/Downloads/scaffold25.LDdecay.stat", header=F)
head(d)
f <- read.table("~/Downloads/OUT.bin", header=F)
head(f)
plot(f$V1, f$V2, type = "l", xlim = c(0, 10000))
