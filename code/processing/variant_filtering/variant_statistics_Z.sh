#!/bin/bash -l
#SBATCH -A snic2020-5-639
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J FilterVariants
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools tabix vcftools 

variant_VCF=$1
output_prefix=$2

# calculate allele frequency
vcftools --gzvcf $variant_VCF --freq2 --out $output_prefix --max-alleles 2
# mean depth of coverage per individual
vcftools --gzvcf $variant_VCF --depth --out $output_prefix
# mean depth per site
vcftools --gzvcf $variant_VCF --site-mean-depth --out $output_prefix
# site quality
vcftools --gzvcf $variant_VCF --site-quality --out $output_prefix
# proportion of missing data per individual
vcftools --gzvcf $variant_VCF --missing-indv --out $output_prefix
# proportion of missing data per site
vcftools --gzvcf $variant_VCF --missing-site --out $output_prefix
# calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $variant_VCF --het --out $output_prefix