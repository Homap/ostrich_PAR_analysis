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


## Linkage disequilibrium
```{r}
data_dir = "../results/manuscript_tables"
LD_file = paste(data_dir,"/LD/LD_chromosome_plot/black.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out", sep="")
black_LD <- read.table(LD_file)
colnames(x = black_LD) = c("Chromosome", "Window_start", "Window_end", "Dprime_LD",
                           "LOD_Dprime", "r_squared", "CIlow_rsq", "CIHi_rsq",
                           "mean_SNP_distance", "NUM_snp_pairs")
# LD decay
scaf26 <- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold26.LDdecay.bin", sep=""))
scaf54<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold54.LDdecay.bin", sep=""))
scaf35<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold35.LDdecay.bin", sep=""))
scaf36<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold36.LDdecay.bin", sep=""))
scaf62<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold62.LDdecay.bin", sep=""))
scaf67<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold67.LDdecay.bin", sep=""))
scaf69_1<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold69-1.LDdecay.bin", sep=""))
scaf93<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold93.LDdecay.bin", sep=""))
scaf63<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold63.LDdecay.bin", sep=""))
scaf88<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold88.LDdecay.bin", sep=""))
scaf83<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold83.LDdecay.bin", sep=""))
scaf92<- read.table(paste(data_dir,"/LD/scaf.split/black.superscaffold92.LDdecay.bin", sep=""))
```

## Mean and sd for LD in PAR and nonPAR
```{r}
# LD for PAR and non-PAR
black_LDPAR <- black_LD[black_LD$Window_start < 53065700,]
black_LDnonPAR <- black_LD[black_LD$Window_start > 53065700,]

# Black mean and sd
mean(black_LD$r_squared, na.rm = T)
sd(black_LD$r_squared, na.rm = T)

mean(black_LDPAR$r_squared)
sd(black_LDPAR$r_squared)

mean(black_LDnonPAR$r_squared, na.rm = T)
sd(black_LDnonPAR$r_squared, na.rm = T)

t.test(black_LDPAR$r_squared, black_LDnonPAR$r_squared)
```

## Comparison of LD decay in the PAR and nonPAR
```{r}
plot(scaf26$V1, scaf26$V2)
PAR_ld <- c(scaf26$V2, scaf54$V2, scaf35$V2)
nonPAR_ld <- c(scaf62$V2, scaf63$V2, scaf67$V2, scaf69_1$V2, scaf83$V2, scaf88$V2, scaf92$V2, scaf93$V2)
# mean of LD decay PAR
mean(PAR_ld)
sd(PAR_ld)
# mean of LD decay nonPAR
mean(nonPAR_ld)
sd(nonPAR_ld)
```

## Measuring LD across the PAR-nonPAR boundary 
```{r}




```