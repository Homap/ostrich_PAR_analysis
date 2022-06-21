# Recombination rate from the genetic map

We use linkage map data from Yazid and Ellegren 2018. 

## Sex-averaged genetic map and smoothed map using loess function
To obtain the sex-averaged genetic map, for the PAR is (male_map + female_map)x0.5
and for the nonPAR with only recombination in male is (male_map)x(2/3).

`Rscript sex_averaged_map.R`

- Table output of sex_averaged_map.R

| Pos     | female_cM | male_cM | sex_averaged_cM | female_smoothed25 | male_smoothed25 | sex_averaged_smoothed25 |
| ------- | --------- | ------- | --------------- | ----------------- | --------------- | ----------------------- | 
| 1113306 | 0         | 0       | 0               | -1.348 | 0.644 | -0.352 | 
| 3461664 | 4.768     | 3.096   | 3.932           | 6.853 | 3.646 | 5.249 | 
| 4790875 | 11.351    | 4.705   | 8.028           | 11.122 | 5.421 | 8.271 | 
| 5218699 | 13.249    | 9.764   | 11.5065         | 12.437 | 6.002 | 9.22 | 

## Kosambi-recombination frequency and the smoothed recombination rate using loess function

`Rscript get_recombination_frequency.R` 

- Recombination rate table

| start_pos | end_pos | female_pair_cm_per_site | male_pair_cm_per_site | sex_averaged_pair_cm_per_site | female_smoothed | male_smoothed | sex_averaged_smoothed |
| ------- | --------- | ------- | --------------- | ----------------- | --------------- | ----------------------- | ---|
| 1113306 | 2113305 | 2.03 |   1.318  | 1.674 |  2.747 |  1.721 |  2.234 |
| 2113306 | 3113305 | 2.03 |   1.318  | 1.674 |  2.855 |  1.653 |  2.254 |

- Map length and recombination rate for the PAR

| female_PAR_length | male_PAR_length | sex_averaged_PAR_length | female_PAR_r | male_PAR_r | sex_averaged_PAR_r |
| ----------------- | --------------- | ----------------------- | ------------ | ---------- | ------------------ |
| 80.628 | 42.641 | 61.6345 | 0.462 | 0.346 | 0.422 |