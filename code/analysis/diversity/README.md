# Measures of genetic variation

Using SNPs, we calculate different measures of genetic variation:

- Average number of pairwise differences (π)
- Average number of segregating sites (θ)
- Tajima's D which reflects the relationship between π and θ
- Folded site frequency spectrum

## Getting allele counts from scaffolds of chromosomes 4, 5, PAR and nonPAR

`bash 1_get_allele_count.sh`

## Calculate measures of genetic variation

`bash 2_get_genetic_variation.sh`

## Measuring genetic variation in synonymous and nonsynonymous sites

`bash 3_get_genetic_variation_CDS.sh`