mb = 10^6
female_map <- read.table(paste(working_dir, "/linkage_map/LGZ3.female.lifted.cleaned.bed.txt", sep = ""))
colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
male_map <- read.table(paste(working_dir, "/linkage_map/LGZ3.male.lifted.cleaned.bed.txt", sep = ""))
colnames(male_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")

sex_averaged_map <- inner_join(female_map, male_map, by = "Pos")
sex_averaged_cM <- (sex_averaged_map$cM.x + sex_averaged_map$cM.y)/2
Pos <- sex_averaged_map$Pos
female_cM <- sex_averaged_map$cM.x
male_cM <- sex_averaged_map$cM.y
sex_averaged_map <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM)
head(sex_averaged_map)

female_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 2, span = 0.3)
male_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 3, span = 0.3)
sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 4, span = 0.3)

rec <- data.frame(
  position=c((female_rate[[2]]$start+female_rate[[2]]$end)/(2*mb),
             (male_rate[[2]]$start+male_rate[[2]]$end)/(2*mb),
             (sex_averaged_rate[[2]]$start+sex_averaged_rate[[2]]$end)/(2*mb)),
  rate=c(female_rate[[2]]$pair_cm_per_site, male_rate[[2]]$pair_cm_per_site, sex_averaged_rate[[2]]$pair_cm_per_site),
  sex=c(rep('female', 78), rep('male', 78), rep('sex_averaged', 78))
)
# 


positions <- sex_averaged_map$Pos/mb
sex_averaged_cM <- sex_averaged_map$sex_averaged_cM
sex_averaged_smooth <- sex_average_map_smoothed25

pos_cm_smooth <- data.frame(positions, sex_averaged_cM, sex_averaged_smooth)
write.table(pos_cm_smooth, "~/Documents/projects/ostrich_Z/results/manuscript_tables/pos_cm_smooth.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Recombination rate in 0.5 Mb windows
start_chr <- sex_averaged_rate[[2]]$start
end_chr <- sex_averaged_rate[[2]]$end
smooth_rate <- sex_averaged_rate[[2]]$pair_cm_per_site
df1 <- data.frame(start_chr, end_chr, smooth_rate, smooth_rate/mb)

sex_averaged_map_1 <- rbind(c(0,0,0,0),sex_averaged_map, c(80543761, 80926604, NA, NA))
sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map_1, len_markers = 195, column = 4, span = 0.3)
test <-c(sex_averaged_rate[[4]],0,0,0,0,0,0,0,0)