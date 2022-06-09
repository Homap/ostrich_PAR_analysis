# Calculate mean LD (r-squared) by sliding window. 
# LD is calculated by PopLDdecay.
# Sliding window setting is similar to this paper: http://www.g3journal.org/content/4/1/121
# "Sliding window analysis (bin 100 kb, step 50 kb) of linkage (r2) for SNPs separated by 1−10 kb. 
#
# Author: takikawakami modified by Homa Papli Yazdi
###############################################################################

#Clear the workspace
rm(list=ls())
#Remove scientific notation
options(scipen = 100)

comarg <- commandArgs(trailingOnly=TRUE) 
in.file <- comarg[1]              # in.file <- "black.superscaffold88.LDdecay.LD.1-10"  
scaf.file <- comarg[2]  			# scaf.file <- "scaf.size.txt"
bin.size <- as.numeric(comarg[3])  			# scaf.file <- "scaf.size.txt"
step.size <- as.numeric(comarg[4])  			# scaf.file <- "scaf.size.txt"


df <- read.table(in.file, header=FALSE)
colnames(df) <- c("chr", "Site1", "Site2","Dprime", "LOD","r2","CIlow","CIhi", "Dist")

# scaf name
scaf_file <- strsplit(in.file[1], "\\/")[[1]][8]
scaf <- strsplit(scaf_file, "\\.")[[1]][1]

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
write.table(df.out, file=paste(in.file, ".", bin.size/1000, "kbBin.", step.size/1000, ".kbStep.out", sep=""), sep='\t', quote=FALSE, row.names = FALSE)



