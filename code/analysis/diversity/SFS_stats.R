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
