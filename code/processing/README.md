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

We need to have a look at the HWE. Most sites are expected to be HWE. Interesting biological features
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

Finally, we need to remove fixed sites for alternative allele since these sites are not variable and therefore cannot be used as SNPs.

`bash 10_VCF_remove_fixed_alternative.sh`

Table 4. Number of SNPs after removing fixed sites for alternative allele
|   | SNP numbers |
| ----------- | ----------- |
| Autosomes | 5776166 |
| PAR | 268006 |
| nonPAR | 89540 |

Final VCF files are
`../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.vcf.gz`
`../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.vcf.gz`
`../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.vcf.gz`

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

## Processing of reference genome by filtering it for coverage

- Create a bed file with variants with coverage less than 5 or more than 70 <br>
```
sbatch get_background_coverage.sh ../../data/coverage/coverage_per_site.txt \
../../data/coverage/coverage_per_site_5_70.txt
```

- Background Filtering for Coverage <br>
```
module load bioinfo-tools BEDTools/2.29.2
bedtools maskfasta -fi ../../data/genome/repeatmask/Struthio_camelus.20130116.OM.fa.masked \
-bed ../../data/coverage/coverage_per_site_5_70.txt -fo ../../data/genome/black.repeat.depth.masked.fa
```

## Z chromosome assembly in this study
The Z chromosome assembly used in this study is composed of 12 scaffolds with the order as follows:

 - Order and length of Z scaffolds from Yazdi and Ellegren 2018

| number | scaffold | Start | end | orientation | segment | scaffold length | Z length |
| ------ | -------- | ----- | --- | ----------- | ------- | --------------- | -------- |
| 1 | superscaffold26 | 0 | 25310599 | + | PAR | 25310599 | 25310599 |
| 2 | superscaffold54 | 16379243 | 0 | - | PAR | 29256470 | 41689842 |
| 3 | superscaffold35 | 0 | 4625539 | + | PAR | 4625539 | 46315381 | 
| 4-a | superscaffold36 | 9394175 | 3524264 | - | PAR | 5869911 | 52185292 |
| gap | superscaffold36 | 3524264 | 3516672 | - | gap | 7592 | 52192884
| 4-b | superscaffold36 | 3516672 | 0 | - | nonPAR | 3516672 | 55709556 |
| 5 | superscaffold62 | 0 | 2917291 | + | nonPAR | 2917291 | 58626847 |
| 6 | superscaffold67 | 0 | 5300260 | + | nonPAR | 5300260 | 63927107 |
| 7 | superscaffold69-1 | 0 | 5978518 | + | nonPAR | 5978518 | 69905625 |
| 8 | superscaffold93 | 0 | 4983591 | + | nonPAR | 4983591 | 74889216 |
| 9 | superscaffold63 | 0 | 1692925 | + | nonPAR | 1692925 | 76582141 |
| 10 | superscaffold88 | 0 | 624114 | + | nonPAR | 624114 | 77206255 |
| 11 | superscaffold83 | 0 | 782506 | + | nonPAR | 782506 | 77988761 |
| 12 | superscaffold92 | 0 | 2882843 | + | nonPAR | 2882843 | 80871604 |

Total length of Z chromosome is 80,871,604 bp. <br>
Total length of PAR is 52,185,292 bp. <br>
The PAR boundary is located at the gap between 52185292 and 52192884. 
Total length of nonPAR is 28,678,720 bp. <br>

There is a gap with Ns between coordinates 3516672 and 3524264 of 7592 bp in superscaffold36.
The gap is excluded when calculating the total length of PAR and nonPAR as it cannot be assigned to neither. 

52,185,292 + 28,678,720 + 7592 = 80,871,604 (Total Z length)

## Conversion of scaffold to chromosome coordinates

Population genetics measure in this paper such as LD, rho, diversity, Tajima's D, etc are reported in windows of a certain size.
The chromosome order and orientation of ostrich Z is driven from genetic map and therefore reflects a true connection among scaffolds. 
To measure window-based statistics, we first convert the coordinates of a given statistic computed within each scaffold into chromosome level coordinates using the script `scaffold_to_chr.py`. In short, to convert to scaffold coordinates to Z chromosome, superscaffold 54 must be cut at position 16379243 and inverted. Superscaffold 36 as a whole must also gets inverted (! Do not invert PAR and nonPAR separately, the whole superscaffold must get inverted in one piece). Throughout the manuscript, we need measures of PAR and nonPAR separately. The conversion from scaffold to chromosome coordinate can therefore be done in three ways:

- Conversion of whole Z <br>
`python scaffold_to_chr.py popgene_measure_file.txt Z > popgene_measure_file.Z.coord.txt`
- Conversion of only PAR <br>
`python scaffold_to_chr.py popgene_measure_file.txt PAR > popgene_measure_file.PAR.coord.txt`
- Conversion of only nonPAR <br>
`python scaffold_to_chr.py popgene_measure_file.txt nonPAR > popgene_measure_file.nonPAR.coord.txt`

## Window-based measures

- Produce the sliding windows:

`python sliding_window.py start_length_file.txt windowsize stepsize > Z.coord.windows.bed`

- To find overlap between our computed popgene measure and a given window:

`bedtools intersect -a Z.coord.windows.bed -b popgene_measure_file.Z.coord.txt  -wao > window_statistic_overlap.txt`

Now that we have the overlap between windows and statistic of interest, we use a suitable script in each case to obtain a window-based average of our computed measure.



