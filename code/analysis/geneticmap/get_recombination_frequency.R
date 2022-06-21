## Recombination frequency

source("../../../software_module/R/recombination_rate.R")

sex_averaged_map <- read.table("../../../data/geneticmap/sex_averaged_map_loess.txt", header = T)

female_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 2, span = 0.3)
male_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 3, span = 0.3)
sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 4, span = 0.3)

kosambi_female <- female_rate[[1]]
kosambi_male <-  male_rate[[1]]
kosambi_sex_averaged <- sex_averaged_rate[[1]]

write.table(kosambi_female, "../../../data/geneticmap/kosambi_female.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(kosambi_male, "../../../data/geneticmap/kosambi_male.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(kosambi_sex_averaged, "../../../data/geneticmap/kosambi_sex_averaged.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

start_pos <- female_rate[[2]]$start
end_pos <- female_rate[[2]]$end
female_pair_cm_per_site <- female_rate[[2]]$pair_cm_per_site
male_pair_cm_per_site <- male_rate[[2]]$pair_cm_per_site
sex_averaged_pair_cm_per_site <- sex_averaged_rate[[2]]$pair_cm_per_site

female_smoothed <- female_rate[[4]]
male_smoothed <- male_rate[[4]]
sex_averaged_smoothed <- sex_averaged_rate[[4]]

recombination_rate <- data.frame(start_pos, end_pos, female_pair_cm_per_site, male_pair_cm_per_site, sex_averaged_pair_cm_per_site, 
                                female_smoothed, male_smoothed, sex_averaged_smoothed)

write.table(recombination_rate, "../../../data/geneticmap/recombination_rate.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

