# Recombination rate from the genetic map

We use linkage map data from Yazid and Ellegren 2018. 

## Sex-averaged genetic map and smoothed map using loess function
To obtain the sex-averaged genetic map, for the PAR is (male_map + female_map)x0.5
and for the nonPAR with only recombination in male is (male_map)x(2/3).

`Rscript sex_averaged_map.R`

## Kosambi-recombination frequency and the smoothed recombination rate using loess function

`Rscript get_recombination_frequency.R` 