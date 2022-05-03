# TODO: Add comment
# 
# Calculate mean LD (r-squared) by sliding window. 
# LD is calculated by PopLDdecay.
# Sliding window setting is similar to this paper: http://www.g3journal.org/content/4/1/121
# "Sliding window analysis (bin 100 kb, step 50 kb) of linkage (r2) for SNPs separated by 1−10 kb in the 20 Tanzanian samples. Windows with less than 100 values were excluded, which resulted in reduced representation in the centromeric regions, where SNP density was lower."
# SNP pairs separated by 1-10kb are used. This is done by awk: cat black.superscaffold88.LDdecay.LD | grep -v "^#" | awk '($9>=1000 && $9<=10000)' > black.superscaffold88.LDdecay.LD.1-10 
#
# Author: takikawakami
###############################################################################

#Clear the workspace
rm(list=ls())

#library(stringr)
#options(scipen=999)

comarg <- commandArgs() 
in.file1 <- comarg[6]      
in.file2 <- comarg[7] 
in.file3 <- comarg[8]         # in.file <- "black.superscaffold88.LDdecay.LD.1-10"  
scaf.file <- comarg[9]  			# scaf.file <- "scaf.size.txt"
bin.size <- as.numeric(comarg[10])  			# scaf.file <- "scaf.size.txt"
step.size <- as.numeric(comarg[11])  			# scaf.file <- "scaf.size.txt"

		
# setwd("/Users/takikawakami/Documents/ostrich_Z/LD_decay")

mean_ld <- function(in.file = in.file, scaf.file = scaf.file, bin.size = bin.size, step.size = step.size){
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

	return(df.out)
}

# Generate the data frames
black <- mean_ld(in.file = in.file1, scaf.file = scaf.file, bin.size = bin.size, step.size = step.size)
blue <- mean_ld(in.file = in.file2, scaf.file = scaf.file, bin.size = bin.size, step.size = step.size)
red <- mean_ld(in.file = in.file3, scaf.file = scaf.file, bin.size = bin.size, step.size = step.size)

#### plot option

# plot relative rho
pdf(paste("allthree", ".", bin.size/1000, "kbBin.", step.size/1000, ".kbStep.pdf", sep=""), width=10, height=8)
par(mfrow=c(2,1), mar=c(2,4,1,0.5), oma=c(3,2,2,1))
plot((black$pos1+black$pos2)/2000000, black$r2, xlab="Physical position (Mb)", ylab="LD (r^2)", pch = 20, col = "black")
points((blue$pos1+blue$pos2)/2000000, blue$r2, pch = 20, col = "blue")
points((red$pos1+red$pos2)/2000000, red$r2, pch = 20, col = "red")
abline(v = 3524263/1000000, col="dimgrey", lty = 2)

plot((black$pos1+black$pos2)/2000000, black$Dprime, xlab="Physical position (Mb)", ylab="LD (Dprime)", pch = 20, col = "black")
points((blue$pos1+black$pos2)/2000000, blue$Dprime, pch = 20, col = "blue")
points((red$pos1+black$pos2)/2000000, red$Dprime, pch = 20, col = "red")
abline(v = 3524263/1000000, col="dimgrey", lty = 2)
dev.off()

