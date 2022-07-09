## Recombination frequency per Mb

library(magrittr)
library(dplyr)
library(reshape2) 
library(devtools)

# Using interpolated data for males, females and the sex-averaged map, I use the 
# Kosambi map function to obtain an estimate of recombination fraction per meiosis.

rate_function <- function(genetic_map, len_markers, column, span, first_marker = 1113307, mb = 10^6){
  start <- genetic_map$start
  end <- genetic_map$end
  pair_cm <- genetic_map[,column]
  length_region <- end - start
  
  # Calculate recombination frequency for each pair of markers
  kosambi_r_length_region = c()
  c = 0
  for( i in pair_cm){
    c = c + 1
    kosambi_r_length_region[c] = 0.5 * ((exp(4*(i/100))-1)/(exp(4*(i/100))+1))
  }
  pair_cm_per_site <- pair_cm/(end-start)
  kosambi_r_per_site <- kosambi_r_length_region/(end-start)
  df <- data.frame(start, end, pair_cm, pair_cm_per_site, kosambi_r_length_region, kosambi_r_per_site, length_region)
  return(df)
}
