#!/usr/bin/R

par_frq <- read.table("../../data/vcf/par_vcf/stats_finalFilter/par.frq", skip = 1)
nonpar_frq <- read.table("../../data/vcf/nonpar_vcf/stats_finalFilter/nonpar.frq", skip = 1)
a_frq <- read.table("../../data/vcf/a_vcf/stats_finalFilter/a.frq", skip = 1)

par_ldepth <- read.table("../../data/vcf/par_vcf/stats_finalFilter/par.ldepth.mean", skip = 1)
nonpar_ldepth <- read.table("../../data/vcf/nonpar_vcf/stats_finalFilter/nonpar.ldepth.mean", skip = 1)
a_ldepth <- read.table("../../data/vcf/a_vcf/stats_finalFilter/a.ldepth.mean", skip = 1)

pdf("../../figures/FigureS2.pdf", w = 9, h = 6.5)
par(mfrow=c(2,3))
hist(par_frq$V6, xlab = "PAR alternatve allele freq.", main = NULL)
hist(a_frq$V6, xlab = "Autosome alternatve allele freq.", main = NULL)
hist(nonpar_frq$V6, xlab = "SLR alternatve allele freq.", main = NULL)

hist(par_ldepth$V3, xlab = "PAR depth per site", main = NULL)
hist(a_ldepth$V3, xlab = "Autosome depth per site", main = NULL)
hist(nonpar_ldepth$V3, xlab = "SLR depth per site", main = NULL)

dev.off()
