# TODO: Add comment
# 
# Calculate mean LD (r-squared) by sliding window. 
# LD is calculated by PopLDdecay.
# Sliding window setting is similar to this paper: http://www.g3journal.org/content/4/1/121
# "Sliding window analysis (bin 100 kb, step 50 kb) of linkage (r2) for SNPs separated by 1−10 kb. 
# Windows with less than 100 values were excluded, which resulted in reduced representation in the centromeric regions, where SNP density was lower."
# SNP pairs separated by 1-10kb are used. This is done by awk: cat black.superscaffold88.LDdecay.LD | grep -v "^#" | awk '($9>=1000 && $9<=10000)' > black.superscaffold88.LDdecay.LD.1-10 
#
# Author: takikawakami
###############################################################################

#Clear the workspace
rm(list=ls())

#library(stringr)
#options(scipen=999)

comarg <- commandArgs() 
in.file <- comarg[6]              # in.file <- "black.superscaffold88.LDdecay.LD.1-10"  
scaf.file <- comarg[7]  			# scaf.file <- "scaf.size.txt"
bin.size <- as.numeric(comarg[8])  			# scaf.file <- "scaf.size.txt"
step.size <- as.numeric(comarg[9])  			# scaf.file <- "scaf.size.txt"

df <- read.table(in.file, header=FALSE)
colnames(df) <- c("chr", "Site1", "Site2","Dprime", "LOD","r2","CIlow","CIhi", "Dist")

# scaf name
scaf <- strsplit(in.file[1], "\\.")[[1]][5]

# scaf list
df.scaf <- read.table(scaf.file, header=TRUE)
colnames(df.scaf) <- c("chr", "pos1", "pos2")

# scaf size
scaf.size <- subset(df.scaf, chr==scaf)$pos2

# set up sliding windows
win.starts <- seq(0, scaf.size-step.size, by=step.size)
win.ends <- win.starts + bin.size
win.ends[length(win.ends)] <- scaf.size

# store results here
chr.out <- rep(scaf, length(win.starts))
nSNP.out <- vector(length=length(win.starts))
Dprime.out <- vector(length=length(win.starts))
LOD.out <- vector(length=length(win.starts))
r2.out <- vector(length=length(win.starts))
CIlow.out <- vector(length=length(win.starts))
CIhi.out <- vector(length=length(win.starts))
meanDist.out <- vector(length=length(win.starts))


# loop
for (i in 1:length(win.starts)) {
	start <- win.starts[i]
	end <- win.ends[i]
	
	# subset LD data
	df.sub <- subset(df, Site1>=start & Site2<=end)
	
	# calculate mean LD etc.
	nSNP.out[i] <- nrow(df.sub)
	Dprime.out[i] <- mean(df.sub$Dprime)
	LOD.out[i] <- mean(df.sub$LOD)
	r2.out[i] <- mean(df.sub$r2)
	CIlow.out[i] <- mean(df.sub$CIlow)
	CIhi.out[i] <- mean(df.sub$CIhi)
	meanDist.out[i] <- mean(df.sub$Dist)
	
}

## output
df.out <- data.frame(chr=chr.out, pos1=win.starts, pos2=win.ends, Dprime=Dprime.out, LOD=LOD.out, r2=r2.out, CIlow=CIlow.out, CIhi=CIhi.out, meanDist=meanDist.out, nSNP=nSNP.out)
write.table(df.out, file=paste(in.file, ".", bin.size/1000, "kbBin.", step.size/1000, ".kbStep.out", sep=""), sep='\t', quote=FALSE, col.names = FALSE, row.names = FALSE)

#### plot option

# plot relative rho
# pdf(paste(in.file, ".", bin.size/1000, "kbBin.", step.size/1000, ".kbStep.pdf", sep=""), width=6, height=6)
# par(mfrow=c(2,1), mar=c(2,4,1,0.5), oma=c(3,2,2,1))
# plot((df.out$pos1+df.out$pos2)/2000000, df.out$r2, xlab="Physical position (Mb)", ylab="LD (r^2)")
# plot((df.out$pos1+df.out$pos2)/2000000, df.out$Dprime, xlab="Physical position (Mb)", ylab="LD (Dprime)")
# dev.off()

