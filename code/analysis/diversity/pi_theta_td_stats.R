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

# Z chromosome coordinates
black_Z <- read.table("black.Z.coordinates.sfs.txt", header=F)
colnames(x = black_Z) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_mid",
                          "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")
Z_500Kb <- read.table("black.Z.500Kb.sfs.txt", header =T)