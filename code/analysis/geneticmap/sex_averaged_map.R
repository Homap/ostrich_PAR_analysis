#!/use/bin/R

# Homa Papoli Yazdi
# Script to obtain the sex-averaged genetic map and the smoothed
# map across the Z chromosome.

library(dplyr)

## Genetic map
female_map <- read.table("../../../data/geneticmap/LGZ3.female.cleaned.chr.bed", sep = "")
colnames(female_map) <- c("Chr", "chrPos", "chrPosPlus1", "scaffold", "scaffoldpos", "scaffoldposplus1", "cM")
male_map <- read.table("../../../data/geneticmap/LGZ3.male.cleaned.chr.bed", sep = "")
colnames(male_map) <- c("Chr", "chrPos", "chrPosPlus1", "scaffold", "scaffoldpos", "scaffoldposplus1", "cM")

both_sex_map <- inner_join(female_map, male_map, by = "chrPosPlus1")

start <- both_sex_map$chrPosPlus1[1:192]
end <- both_sex_map$chrPosPlus1[2:193]
female_cM <- round(abs(both_sex_map$cM.x[2:193] - both_sex_map$cM.x[1:192]), 3)
male_cM <- round(abs(both_sex_map$cM.y[2:193] - both_sex_map$cM.y[1:192]), 3)
Pos <- (start + end)/2
first_scaffold <- both_sex_map$scaffold.x[1:192]
second_scaffold <- both_sex_map$scaffold.x[2:193]
scaffold_start <- both_sex_map$scaffoldposplus1.x[1:192]
scaffold_end <- both_sex_map$scaffoldposplus1.x[2:193]

rate_data <- data.frame(start, end, Pos, female_cM, male_cM, first_scaffold, second_scaffold, scaffold_start, scaffold_end)

par <- rate_data[which(rate_data$start < 53050777),]
nonpar <- rate_data[which(rate_data$start >= 53050777),]

sex_averaged_cM <- c((par$female_cM + par$male_cM)*0.5, (nonpar$male_cM)*(2/3))

sex_averaged_map <- data.frame(start, end, Pos, female_cM, male_cM, sex_averaged_cM, first_scaffold, second_scaffold, scaffold_start, scaffold_end)

# LOESS for smoothing of the three maps
female_map_loessMod25 <- loess(female_cM ~ Pos, data=sex_averaged_map, span=0.25) # 25% smoothing span
male_map_loessMod25 <- loess(male_cM ~ Pos, data=sex_averaged_map, span=0.25) # 25% smoothing span
sex_averaged_map_loessMod25 <- loess(sex_averaged_cM ~ Pos, data=sex_averaged_map, span=0.25) # 25% smoothing span

female_map_smoothed25 <- predict(female_map_loessMod25)
male_map_smoothed25 <- predict(male_map_loessMod25)
sex_average_map_smoothed25 <- predict(sex_averaged_map_loessMod25)

female_smoothed25 <- round(female_map_smoothed25, 3)
male_smoothed25 <- round(male_map_smoothed25, 3)
sex_averaged_smoothed25 <- round(sex_average_map_smoothed25, 3)

chrom = rep("ChrZ", 192)

sex_averaged_map_loess <- data.frame(chrom, start, end, Pos, female_cM, male_cM,sex_averaged_cM,
                                    first_scaffold, second_scaffold, scaffold_start, scaffold_end,
                                    female_smoothed25, male_smoothed25, sex_averaged_smoothed25)

write.table(sex_averaged_map_loess, "../../../data/geneticmap/sex_averaged_map_loess.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)