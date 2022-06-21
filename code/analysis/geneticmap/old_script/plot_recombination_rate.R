#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
df_male = read.table(args[1], header=FALSE)
df_female = read.table(args[2], header=FALSE)

# scaf_vector = unique(df$V1)

output=strsplit(strsplit(args[1], "\\_")[[1]][6], "\\.")[[1]][1]
print(output)

# print(scaf_vector)

# for(scaf in scaf_vector){
# 	print(scaf)
# 	df_scaf = subset(df, df$V1==scaf)
# 	print(df_scaf)
#  	pdf(paste(scaf, ".", output, ".pdf", sep=""))
#  	plot(df_scaf$V2, df_scaf$V3, type = "l", xlab = "Physical position (Mb)", ylab = "Recombination rate")
#  	dev.off()
# }

df_scaf_male = subset(df_male, df_male$V1=="superscaffold36")
df_scaf_female = subset(df_female, df_female$V1=="superscaffold36")

pdf(paste("superscaffold36", ".", output, ".pdf", sep=""), width = 12, height = 8)
plot(df_scaf_male$V2/1000000, df_scaf_male$V3, type = "l", xlab = "Physical position (Mb)", ylab = "Recombination rate", col = "blue")
lines(df_scaf_female$V2/1000000, df_scaf_female$V3, col = "red")
abline(v = 3524263/1000000, col="dimgrey", lty = 2)
dev.off()