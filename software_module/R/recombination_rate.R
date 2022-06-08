rate_function <- function(genetic_map, len_markers, column, span, first_marker = 1113306, mb = 10^6){
  start <- genetic_map$Pos[1:len_markers-1]
  end <- genetic_map$Pos[2:len_markers]-1
  pair_cm <- abs(diff(genetic_map[,column]))
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
  df1 <- df %>% 
    rowwise() %>% 
    mutate(row = list(seq(from = start, to = end))) %>% 
    ungroup() %>% 
    tidyr::unnest(row) %>% 
    mutate(group = (row - first_marker) %/% mb)
  
  df2 <- df1 %>% 
    group_by(group) %>% 
    summarise(
      start = min(row), 
      end = max(row), 
      pair_cm_per_site = mean(pair_cm_per_site)*mb, 
      .groups = "drop"
    ) %>% 
    select(-group)
  rate_loess <- loess(pair_cm_per_site ~ start, data=df2, span=span, na.action = na.exclude)
  rate_smoothed <- predict(rate_loess)
  return(list(df, df2, rate_loess, rate_smoothed))
}