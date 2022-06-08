# Examining the vcf statistics produced by vcftools and generating figures

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

variant_statistics <- function(site_quality, site_depth, site_missing, site_frequency) {
    # Takes vcftools ouputs for site quality as input and outputs
    # graphs summarizing statistics.
    # Variant quality
    var_qual <- read_delim(site_quality, delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
    var_qual_fig <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    var_qual_fig + theme_light() + xlim(0, 10000)

    # Variant mean depth
    var_depth <- read_delim(site_depth, delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
    var_depth_fig <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    var_depth_fig + theme_light() + xlim(0, 100)
    var_depth_summary <- summary(var_depth$mean_depth)

    # Variant missingness
    var_miss <- read_delim(site_missing, delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
    var_miss_fig <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    var_miss_fig + theme_light()

    # Minor allele frequency
    var_freq <- read_delim(site_frequency, delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
    # find minor allele frequency
    var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
    var_freq_fig <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    var_freq_fig + theme_light()

    return(list(var_qual_fig, var_depth_fig, var_miss_fig, var_freq_fig, var_depth_summary))
}

individual_statistics <- function(ind_depth, ind_miss, ind_het){
    # Takes vcftools ouputs for individual quality as input and outputs
    # graphs summarizing statistics.
    # Mean depth per individual
    ind_depth <- read_delim(ind_depth, delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
    ind_depth_fig <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    ind_depth_fig + theme_light()

    # Proportion of missing data per individual
    ind_miss <- read_delim(ind_miss, delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
    ind_miss_fig <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    ind_miss_fig + theme_light()

    # Heterozygosity and inbreeding coefficient per individual
    ind_het <- read_delim(ind_het, delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
    ind_het_fig <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
    ind_het_fig + theme_light()

    return(list(ind_depth_fig, ind_miss_fig, ind_het_fig))
}