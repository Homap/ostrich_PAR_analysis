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
female_pair_cm_per_mb <- round(female_rate[[2]]$pair_cm_per_site, 3)
male_pair_cm_per_mb <- round(male_rate[[2]]$pair_cm_per_site, 3)
sex_averaged_pair_cm_per_mb <- round(sex_averaged_rate[[2]]$pair_cm_per_site, 3)

female_smoothed <- round(female_rate[[3]], 3)
male_smoothed <- round(male_rate[[3]], 3)
sex_averaged_smoothed <- round(sex_averaged_rate[[3]], 3)

recombination_rate <- data.frame(start_pos, end_pos, female_pair_cm_per_mb, male_pair_cm_per_mb, sex_averaged_pair_cm_per_mb, 
                                female_smoothed, male_smoothed, sex_averaged_smoothed)

write.table(recombination_rate, "../../../data/geneticmap/recombination_rate.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Genetic length in PAR in females, males and sex-averaged where PAR is defined at 53065777
female_PAR_length = 80.628  
male_PAR_length = 42.641
sex_averaged_PAR_length = 61.6345

female_PAR_r = round((0.5 * ((exp(4*(female_PAR_length/100))-1)/(exp(4*(female_PAR_length/100))+1))), 3)
male_PAR_r = round((0.5 * ((exp(4*(male_PAR_length/100))-1)/(exp(4*(male_PAR_length/100))+1))), 3)
sex_averaged_PAR_r = round((0.5 * ((exp(4*(sex_averaged_PAR_length/100))-1)/(exp(4*(sex_averaged_PAR_length/100))+1))), 3)

par_map_rate <- data.frame(female_PAR_length, male_PAR_length, sex_averaged_PAR_length,
                           female_PAR_r, male_PAR_r, sex_averaged_PAR_r)

write.table(par_map_rate, "../../../data/geneticmap/par_map_rate.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


