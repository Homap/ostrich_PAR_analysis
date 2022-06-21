#!/use/bin/R

library(dplyr)

## Genetic map
female_map <- read.table("../../../data/geneticmap/LGZ3.female.lifted.cleaned.bed.txt", sep = "")
colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
male_map <- read.table("../../../data/geneticmap/LGZ3.male.lifted.cleaned.bed.txt", sep = "")
colnames(male_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")

sex_averaged_map <- inner_join(female_map, male_map, by = "Pos")

par <- sex_averaged_map[which(sex_averaged_map$Pos < 53065777),]
nonpar <- sex_averaged_map[which(sex_averaged_map$Pos >= 53065777),]

sex_averaged_cM <- c((par$cM.x + par$cM.y)*0.5, (nonpar$cM.y)*(2/3))
Pos <- sex_averaged_map$Pos
female_cM <- sex_averaged_map$cM.x
male_cM <- sex_averaged_map$cM.y
sex_averaged_map <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM)

write.table(sex_averaged_map, "../../../data/geneticmap/sex_averaged_map.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)