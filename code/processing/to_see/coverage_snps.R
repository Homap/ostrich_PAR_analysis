#!/usr/bin/env Rscript

# Written by Homa Papoli

args = commandArgs(trailingOnly=TRUE)
df = read.table(args[1], header=TRUE)

mean_dp <- mean(df$MEAN_DEPTH)
sd_dp <- sd(df$MEAN_DEPTH)

print(mean_dp)
print(sd_dp)

low_threshold <- mean_dp/3
up_threshold <- mean_dp*2

print(low_threshold)
print(up_threshold)

low_cov_sites <- which(df$MEAN_DEPTH < low_threshold)
#print(length(low_cov_sites))
print(length(data.frame(df[-low_cov_sites, ])$MEAN_DEPTH))
high_cov_sites <- which(df$MEAN_DEPTH > up_threshold)
print(length(data.frame(df[-high_cov_sites, ])$MEAN_DEPTH))
#print(length(high_cov_sites))
sites_to_filter <- c(low_cov_sites, high_cov_sites)

filtered_df <- data.frame(df[-sites_to_filter, ])
print(length(filtered_df$MEAN_DEPTH))

#chrom = filtered_df$CHROM
#chromStart = filtered_df$POS-1
#chromEnd = filtered_df$POS
#filtered_df_bed <- data.frame(chrom, chromStart, chromEnd)
# filtered_df_out = data.frame(filtered_df)

write.table(filtered_df, args[2], sep = "\t", dec = ".",
            row.names = FALSE, quote=FALSE,
            eol = "\n")

pdf(args[3])
hist(filtered_df$MEAN_DEPTH, breaks = 50, main = args[4], xlab = "Depth per site")
dev.off()
