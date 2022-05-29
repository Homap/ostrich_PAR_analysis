# Processing data

This document describes the steps for reproducing the filtered data used for analyses in this project.
`/proj/snic2020-16-269/private/cornwallis.2020/results/ind/analysis/gatk_best_practice_snp/f03_concat_bcftools` produced by Per Unneberg

## Steps for filtering the VCF
To filter the VCF file for autosome, PAR and nonPAR, the following scrips need to be run in the order numbered.

```
sbatch 1_VCF_into_A_PAR_nonPAR.sh
sbatch 2_keepSNPs_initialFilter.sh
sbatch 3_VCF_stats_pre_filter.sh
sbatch 4_VCF_stats_summary.R
sbatch 5_VCF_quality_filter.sh
sbatch 6_VCF_stats_postFilter.sh
```

Table 1. Number of SNPs after filtering
|   | SNP numbers |
| ----------- | ----------- |
| Autosomes | 6524315 |
| PAR | 301807 |
| nonPAR | 117207 |

In addition to the filtering above, several other filtering should also be performed. 

Filtering SNPs overlapping with repetitive elements

`sbatch 7_VCF_mask_repeats.sh`

Table 2. Number of SNPs after removing variant overlapping repeats
|   | SNP numbers |
| ----------- | ----------- |
| Autosomes | 6143527 |
| PAR | 285271 |
| nonPAR | 107894 |

Identifying SNPs with heterozygous genotypes in females in nonPAR as they are likely in the gametolog region

In the nonPAR, females have only 1 Z, we therefore do not expect any female individual with genotype 0/1. Some of these SNPs could be 
located in gametologous genes. If so, they do not reflect the polymorphisms segregating in the population and are substitutions between the
Z and W chromosomes.

`bash 8_VCF_remove_nonPAR_het_females.sh`

Table 2. Number of SNPs after removing variant overlapping repeats and heterzygous sites in female nonPAR
|   | SNP numbers |
| ----------- | ----------- |
| Autosomes | 6143527 |
| PAR | 285271 |
| nonPAR | 104540 |

Finally, we need to have a look at the HWE. Most sites are expected to be HWE. Interesting biological features
such as population structure can cause a deviation from HWE. However, genotyping errors can also cause 
deviations from HWE. We remove SNPs with a HWE p-value below 0.005. For nonpar, since the homozygosity in females
can cause deviations from the HWE, we detected deviations from the male samples only and removed sites from the nonpar
VCF subsequently.

`sbatch 9_VCF_remove_HWE_0.005.sh`

Table 3. Number of SNPs after removing variant overlapping repeats and heterzygous sites in female nonPAR and HWE deviation
|   | SNP numbers |
| ----------- | ----------- |
| Autosomes | 6141350 |
| PAR | 285186 |
| nonPAR | 104540 |

The filtered VCF files for autosomes, PAR and nonPAR are used for the analysis in this project:
- Linkage Disequilibrium (LD)
- Population scaled recombination rate (rho)
- Genetic diversity measures:
    - Pairwise nucleotide diversity (Pi)
    - Watterson's theta (theta)
    - Tajima's D
    - Site Frequency Spectrum (SFS)





