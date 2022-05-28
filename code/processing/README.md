# Processing of data

This document describes the steps for reproducing the filtered data used for analyses in this project.
/proj/snic2020-16-269/private/cornwallis.2020/results/ind/analysis/gatk_best_practice_snp/f03_concat_bcftools produced by Per Unneberg

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

Filtering of SNPs with heterozygous genotypes in females in nonPAR as they are likely in the gametolog region

In the nonPAR, females have only 1 Z, we therefore do not expect any female individual with genotype 0/1. Some of these SNPs could be 
located in gametologous genes. If so, they do not reflect the polymorphisms segregating in the population and are substitutions between the
Z and W chromosomes.

```
python check_nonPAR_genotype.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz > ../data/bed/${species}.nonPAR.female_het.bed
python vcf_to_genotypes.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz ../data/bed/${species}.nonPAR.genotypes.bed
```

## Two counts files for males and females in the nonPAR
```
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_male.txt --out ../data/allele_count/${species}.nonPAR.filtered.male
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_female.txt --out ../data/allele_count/${species}.nonPAR.filtered.female
```

## Filter out heterozgous sites in females from both males and females
```
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.male.frq.count > ../data/allele_count/${species}.nonPAR.filtered.male.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.male.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count

awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.female.frq.count > ../data/allele_count/${species}.nonPAR.filtered.female.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.female.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count

awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count > ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count.bed -b ../data/bed/${species}.nonPAR.female_het.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""femalehet"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.male.repeat.femalehet.frq.count

awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count > ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count.bed -b ../data/bed/${species}.nonPAR.female_het.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""femalehet"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.female.repeat.femalehet.frq.count
```

## Obtain nonPAR allele count
The nonPAR allele count should be adjusted for females having 5 chromosomes. In VCF, females are called as diploid
in the nonPAR so 0/0 and 1/1 should be changed into 0 and 1.

```
python nonPAR_allele_count.py ../data/allele_count/${species}.nonPAR.filtered.female.repeat.femalehet.frq.count \
../data/allele_count/${species}.nonPAR.filtered.male.repeat.femalehet.frq.count > ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count
```


After filtering, SNPs in different functional categories need also to be separated into different files.
Intergenic
Intronic
CDS







# Match chicken 1 to 5 with ostrich scaffolds
python code/processing/chicken_to_ostrich_autosomes.py data/lastz/chicken_NC_chr.txt data/lastz/chicken_ostrich.lastz data/genome/Struthio_camelus.20130116.OM.fa.fai > data/lastz/gg_ostrich_macrochr.txt

# Create bed files with ostrich scaffolds matching chromosomes 1 to 5 of chicken
for chrom in {1..5}
do
echo chromosome_${chrom}
grep  chromosome_${chrom} data/lastz/gg_ostrich_macrochr.txt | awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $4"\t""0""\t"$5}' | uniq > data/bed/gg_chr${chrom}_ostrich.bed
awk '{if(NR > 1) sum+=$3} END{print sum}' data/bed/gg_chr${chrom}_ostrich.bed
done






```
export reference=/proj/snic2020-16-269/private/homap/ostrich_z/data/reference/Struthio_camelus.20130116.OM.fa
export vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf
export black=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black.txt
export black_male=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black_male.txt
export black_female=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black_female.txt
```

## Create a bed file with PAR and non-PAR regions and autosomes
```
python scaf_length_to_bed.py $reference ../data/bed/z_scaf.bed \
../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed
grep -v -f ../data/bed/Z_scaffolds.txt ../data/reference/Struthio_camelus.20130116.OM.fa.fai \
| awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $1"\t""0""\t"$2}' > ../data/bed/autosomes.bed
```








## Background Filtering for Coverage
```
sbatch get_depth.sh ../data/vcf/${species}.all.vcf.gz ../data/vcf/${species}.all
```
# Get a list of sites with coverage below or above the threshold used for SNP filtering
```
sbatch get_background_coverage.sh ../data/vcf/${species}.all.ldepth.mean ../data/vcf/${species}.all.ldepth.mean.filter.bed
```

## Obtain data for the plot of boundary coverage
```
mkdir -p ../result/boundary_coverage
module load bioinfo-tools samtools
bamaddr=../../../cornwallis.2020/data/interim/map/bwa/dedup 
samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_107.bam ${bamaddr}/P1878_108.bam ${bamaddr}/P1878_109.bam ${bamaddr}/P1878_110.bam ${bamaddr}/P1878_111.bam \
> ../result/boundary_coverage/male_boundary_coverage.txt

samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_112.bam ${bamaddr}/P1878_113.bam ${bamaddr}/P1878_114.bam ${bamaddr}/P1878_115.bam ${bamaddr}/P1878_116.bam \
> ../result/boundary_coverage/female_boundary_coverage.txt

scp homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/result/boundary_coverage/male_boundary_coverage.txt
scp homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/result/boundary_coverage/female_boundary_coverage.txt
```

Use the repeat masked genome where repeats are already N. For analysis where you need the number 
of filtered invariant sites, you can simply count non-N sites in each window and that will be 
the number of bases in window.

```
module load bioinfo-tools BEDTools
bedtools maskfasta -fi ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked \
-bed ../data/vcf/${species}.all.ldepth.mean.filter.bed -fo ../data/reference/${species}.repeat.depth.masked.fa
```

## Get allele count
```
module load bioinfo-tools vcftools BEDTools
echo "Autosome"
vcftools --gzvcf ../data/vcf/${species}.A.filtered.vcf.gz --counts --out ../data/allele_count/${species}.A
echo "PAR"
vcftools --gzvcf ../data/vcf/${species}.PAR.filtered.vcf.gz --counts --out ../data/allele_count/${species}.PAR
echo "nonPAR"
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --out ../data/allele_count/${species}.nonPAR
```






## Find overlap between windows and each of the functional categories
```
sbatch window_overlap.sh
```
