#!/use/bin/R

# Homa Papoli Yazdi
# Script to obtain the sex-averaged genetic map and the smoothed
# map across the Z chromosome.

library(dplyr)

## Genetic map
female_map <- read.table("../../../data/geneticmap/LGZ3.female.lifted.cleaned.bed.txt", sep = "")
colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
male_map <- read.table("../../../data/geneticmap/LGZ3.male.lifted.cleaned.bed.txt", sep = "")
colnames(male_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")

both_sex_map <- inner_join(female_map, male_map, by = "Pos")

par <- both_sex[which(both_sex$Pos < 53065777),]
nonpar <- both_sex[which(both_sex$Pos >= 53065777),]

sex_averaged_cM <- c((par$cM.x + par$cM.y)*0.5, (nonpar$cM.y)*(2/3))
Pos <- both_sex$Pos
female_cM <- both_sex$cM.x
male_cM <- both_sex$cM.y

sex_averaged_map <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM)

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

sex_averaged_map_loess <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM, 
                                    female_smoothed25, male_smoothed25, sex_averaged_smoothed25)

write.table(sex_averaged_map_loess, "../../../data/geneticmap/sex_averaged_map_loess.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)