# Recombination rate from the genetic map

We use linkage map data from Yazid and Ellegren 2018 (Genome Biology and Evolution).  

## Convert to chromosome coordinates
`python scaffold_to_chr.py ../../../data/geneticmap/LGZ3.female.cleaned.bed Z > ../../../data/geneticmap/LGZ3.female.cleaned.chr.bed`
`python scaffold_to_chr.py ../../../data/geneticmap/LGZ3.male.cleaned.bed Z > ../../../data/geneticmap/LGZ3.male.cleaned.chr.bed`

## Sex-averaged genetic map and smoothed map using loess function
To obtain the sex-averaged genetic map, for the PAR is (male_map + female_map)x0.5
and for the nonPAR with only recombination in male is (male_map)x(2/3).

`Rscript sex_averaged_map.R`

- Table output of sex_averaged_map.R

| chrom | start | end | Pos | female_cM | male_cM | sex_averaged_cM | first_scaffold | second_scaffold | scaffold_start | scaffold_end | female_smoothed25 | male_smoothed25 | sex_averaged_smoothed25 |
| ----- | ----- | --- | --- | --------- | ------- | --------------- | -------------- | --------------- | -------------- | ------------ | ----------------- | --------------- | ----------------------- |
| ChrZ | 1113307 | 3461665 |2287486 | 4.768 | 3.096 | 3.932 | superscaffold26 | superscaffold26 | 1113307 | 3461665 | 6.535 | 3.707 | 5.121 |
| ChrZ | 3461665 | 4790876 |4126270.5 | 6.583 | 1.609 | 4.096 | superscaffold26 | superscaffold26 | 3461665 | 4790876 | 5.401 | 2.986 | 4.194 |

## Kosambi-recombination frequency and the smoothed recombination rate using loess function

`Rscript get_recombination_frequency.R` 

- Kosambi sex-averaged recombination frequency

| start |  end    | pair_cm |pair_cm_per_site     |   kosambi_r_length_region| kosambi_r_per_site   |   length_region|
| ----- | ------- | ------- | ------------------- | ------------------------ | -------------------- | -------------- |
| 1113306 |3461663| 3.932  | 1.674e-06  |  0.03924    |  1.6709e-08 |   2348357|
| 3461664 |4790874 |4.096  | 3.0815e-06  |  0.04087    |  3.0746e-08 |   1329210|
| 4790875 |5218698 |3.4785 | 8.1307e-06  |  0.03473    |  8.1177e-08  |  427823|

## Window-based recombination rate

```
module load bioinfo-tools BEDTools/2.29.2

awk 'BEGIN{OFS="\t"}NR==1{$(NF+1)="CHROM"} NR>1{$(NF+1)="ChrZ"} {chr=$NF; print chr,$1,$2,$3,$4,$5,$6,$7}' ../../../data/geneticmap/kosambi_sex_averaged.txt > temp && mv -f temp ../../../data/geneticmap/kosambi_sex_averaged.txt

# 1Mb window
bedtools intersect -a ../../../data/sliding_window/Z.coord.1Mb.windows.PAR.txt -b ../../../data/geneticmap/kosambi_sex_averaged.txt -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""pair_cm""\t""pair_cm_per_site""\t""kosambi_r_length_region""\t""kosambi_r_per_site""\t""length_region"}{print $0}' > ../../../data/geneticmap/kosambi_sex_averaged.bedtools.txt

python get_recombination_per_window.py ../../../data/geneticmap/kosambi_sex_averaged.bedtools.txt ../../../data/sliding_window/Z.coord.1Mb.windows.PAR.txt > ../../../data/geneticmap/kosambi_sex_averaged_1MB.window.PAR.txt

rm -f ../../../data/geneticmap/kosambi_sex_averaged.bedtools.txt 
```

- Map length and recombination rate for the PAR

| female_PAR_length | male_PAR_length | sex_averaged_PAR_length | female_PAR_r | male_PAR_r | sex_averaged_PAR_r |
| ----------------- | --------------- | ----------------------- | ------------ | ---------- | ------------------ |
| 80.628 | 42.641 | 61.6345 | 0.462 | 0.346 | 0.422 |