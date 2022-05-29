#!/bin/bash -l
#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J FilterVariants
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools tabix vcftools vcflib

variant_VCF=$1
variant_subset=$2
output_prefix=$3

# subset vcf
bcftools view $variant_VCF | vcfrandomsample -r 0.05 > $variant_subset
# compress vcf
bgzip $variant_subset
# index vcf
bcftools index ${variant_subset}.gz
# calculate allele frequency
vcftools --gzvcf ${variant_subset}.gz --freq2 --out $output_prefix --max-alleles 2
# mean depth of coverage per individual
vcftools --gzvcf ${variant_subset}.gz --depth --out $output_prefix
# mean depth per site
vcftools --gzvcf ${variant_subset}.gz --site-mean-depth --out $output_prefix
# site quality
vcftools --gzvcf ${variant_subset}.gz --site-quality --out $output_prefix
# proportion of missing data per individual
vcftools --gzvcf ${variant_subset}.gz --missing-indv --out $output_prefix
# proportion of missing data per site
vcftools --gzvcf ${variant_subset}.gz --missing-site --out $output_prefix
# calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf ${variant_subset}.gz --het --out $output_prefix
