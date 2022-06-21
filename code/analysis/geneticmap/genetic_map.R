
rm(list=ls())



## Genetic map
# mb = 10^6
# female_map <- read.table("../../../data/geneticmap/LGZ3.female.lifted.cleaned.bed.txt", sep = "")
# colnames(female_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")
# male_map <- read.table("../../../data/geneticmap/LGZ3.male.lifted.cleaned.bed.txt", sep = "")
# colnames(male_map) <- c("Chr", "Pos", "PosPlus1", "mapname", "cM", "scaffold")

# sex_averaged_map <- inner_join(female_map, male_map, by = "Pos")

# par <- sex_averaged_map[which(sex_averaged_map$Pos < 53065777),]
# nonpar <- sex_averaged_map[which(sex_averaged_map$Pos => 53065777),]

# sex_averaged_cM <- c((par$cM.x + par$cM.y)*0.5, (nonpar$cM.y)*(2/3))
# Pos <- sex_averaged_map$Pos
# female_cM <- sex_averaged_map$cM.x
# male_cM <- sex_averaged_map$cM.y
# sex_averaged_map <- data.frame(Pos, female_cM, male_cM,sex_averaged_cM)
# head(sex_averaged_map)

# write.table(sex_averaged_map, "../../../data/geneticmap/sex_averaged_map.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ## LOESS for smoothing of the three maps
# female_map_loessMod25 <- loess(cM ~ Pos, data=female_map, span=0.25) # 25% smoothing span
# male_map_loessMod25 <- loess(cM ~ Pos, data=male_map, span=0.25) # 25% smoothing span
# sex_averaged_map_loessMod25 <- loess(sex_averaged_cM ~ Pos, data=sex_averaged_map, span=0.25) # 25% smoothing span


# female_map_smoothed25 <- predict(female_map_loessMod25)
# male_map_smoothed25 <- predict(male_map_loessMod25)
# sex_average_map_smoothed25 <- predict(sex_averaged_map_loessMod25)

# ## Plotting genetic map 
# # Shade the area where PAR-nonPAR boundary is residing
# plot(NULL, xlim = c(0, 82), ylim = c(0, 95), ylab = "Genetic map position (cM)", xlab = "Physical position (Mb)")
# rect(51738942/mb, -5, 53065777/mb, 100, lty = 2, col = "lightgrey")

# # Plot the points for the three genetic maps, female, male and sex-averaged
# points(male_map$Pos/mb, male_map$cM, col = "blue")
# lines(male_map$Pos/mb, male_map_smoothed25, col = "blue")
# points(female_map$Pos/mb, female_map$cM, col = "red")
# lines(female_map$Pos/mb, female_map_smoothed25, col = "red")
# points(sex_averaged_map$Pos/mb, sex_averaged_map$sex_averaged_cM, col = "darkorchid1")
# lines(sex_averaged_map$Pos/mb, sex_average_map_smoothed25, col = "darkorchid1")
# legend(x = -2, y = 98, legend = c("male", "female", "sex-averaged"), fill = c("blue", "red", "darkorchid1"), box.lty = 0)


## Recombination frequency
# female_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 2, span = 0.3)
# male_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 3, span = 0.3)
# sex_averaged_rate <- rate_function(genetic_map = sex_averaged_map, len_markers = 193, column = 4, span = 0.3)

# head(female_rate[[1]])
# head(male_rate[[1]])
# head(sex_averaged_rate[[1]])

# plot(NULL, xlim = c(0, 82), ylim = c(0,7), ylab = "cM/Mb", xlab = "Physical position (Mb)")
# rect(51738942/mb, -5, 53065777/mb, 100, lty = 2, col = "lightgrey")

# points((female_rate[[2]]$start+female_rate[[2]]$end)/(2*mb), female_rate[[2]]$pair_cm_per_site, col = "red", pch = 20)
# lines((female_rate[[2]]$start+female_rate[[2]]$end)/(2*mb), female_rate[[4]], col = "red")

# points((male_rate[[2]]$start+male_rate[[2]]$end)/(2*mb), male_rate[[2]]$pair_cm_per_site, col = "blue", pch = 20)
# lines((male_rate[[2]]$start+male_rate[[2]]$end)/(2*mb), male_rate[[4]], col = "blue")

# points((sex_averaged_rate[[2]]$start+sex_averaged_rate[[2]]$end)/(2*mb), sex_averaged_rate[[2]]$pair_cm_per_site, col = "darkorchid1", pch = 20)
# lines((sex_averaged_rate[[2]]$start+sex_averaged_rate[[2]]$end)/(2*mb), sex_averaged_rate[[4]], col = "darkorchid1")



