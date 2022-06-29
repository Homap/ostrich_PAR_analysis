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

The PAR and nonPAR VCF files are concatenated into Z VCF:

```
module load bioinfo-tools tabix/0.2.6 vcftools/0.1.16 bcftools/1.14
export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
mkdir -p ../../data/vcf/z_vcf
vcf-concat ${par_vcf} ${nonpar_vcf} | bgzip -c > ../../data/vcf/z_vcf/z_vcf.gz
bcftools sort ../../data/vcf/z_vcf/z_vcf.gz -o ../../data/vcf/z_vcf/sorted.z.vcf
bgzip ../../data/vcf/z_vcf/sorted.z.vcf
rm -f ../../data/vcf/z_vcf/z_vcf.gz
```

The filtered VCF files for autosomes, PAR and nonPAR are used for the analysis in this project:
- Linkage Disequilibrium (LD)
- Population scaled recombination rate (rho)
- Genetic diversity measures:
    - Pairwise nucleotide diversity (Pi)
    - Watterson's theta (theta)
    - Tajima's D
    - Site Frequency Spectrum (SFS)

## Conversion of scaffold to chromosome coordinates 
The Z chromosome assembly used in this study is composed of 12 scaffolds with the order as follows:

 - Order and length of Z scaffolds from Yazdi and Ellegren 2018

| number | scaffold | Start | end | orientation | segment | length |
| ------ | -------- | ----- | --- | ----------- | ------- | ------ |
| 1 | superscaffold26 | 0 | 25310599 | + | PAR | 25310599 |
| 2 | superscaffold54 | 0 | 16379243 | - | PAR | 29256470 |
| 3 | superscaffold35 | 0 | 4625539 | + | PAR | 4625539 |
| 4-a | superscaffold36 | 0 | 3516673 | - | nonPAR | 3516673 |
| 4-b | superscaffold36 | 3524263 | - | 9394175 | PAR | 5869912 |
| 5 | superscaffold62 | 0 | 2917291 | + | nonPAR | 2917291 |
| 6 | superscaffold67 | 0 | 5300260 | + | nonPAR | 5300260 |
| 7 | superscaffold69-1 | 0 | 5978518 | + | nonPAR | 5978518 |
| 8 | superscaffold93 | 0 | 4983591 | + | nonPAR | 4983591 |
| 9 | superscaffold63 | 0 | 1692925 | + | nonPAR | 1692925 |
| 10 | superscaffold88 | 0 | 624114 | + | nonPAR | 624114 |
| 11 | superscaffold83 | 0 | 782506 | + | nonPAR | 782506 |
| 12 | superscaffold92 | 0 | 2882843 | + | nonPAR | 2882843 |

Total length of Z chromosome is 80,871,604 bp.
Total length of PAR is 52,185,293 bp.
Total length of nonPAR is 28,678,721 bp.

There is a gap with Ns between coordinates 3516673 and 3524263 of 7590 bp in superscaffold36.
The gap is excluded when calculating the total length of PAR and nonPAR as it cannot be assigned to neither. 

52,185,293 + 28,678,721 + 7590 = 80,871,604 (Total Z length)

To convert to Z chromosome, superscaffold 54 must be cut at position 16379243 and inverted.
Superscaffold 36 as a whole must also gets inverted (! Do not invert PAR and nonPAR separately, the whole superscaffold must get inverted in one piece).



