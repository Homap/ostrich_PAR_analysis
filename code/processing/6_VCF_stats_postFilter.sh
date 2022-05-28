#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J vcf_stats
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools/1.14 tabix/0.2.6 vcftools/0.1.16 

export a_vcf=../../data/vcf/a_vcf/a_vcf.filtered.vcf.gz
export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.vcf.gz

mkdir -p ../../data/vcf/a_vcf/stats_postFilter
mkdir -p ../../data/vcf/par_vcf/stats_postFilter
mkdir -p ../../data/vcf/nonpar_vcf/stats_postFilter

for variant_VCF in ${nonpar_vcf} ${par_vcf} ${a_vcf} 
do
output_prefix=$(echo $variant_VCF | cut -f6 -d "/" | cut -f1 -d "_")
output_path=$(echo $variant_VCF | cut -f1-5 -d "/")
echo $variant_VCF $output_prefix $output_path
# Count the number of $variant_VCF variants
bcftools view -H $variant_VCF | wc -l > ${output_path}/stats_postFilter/${output_prefix}.allCounts
# calculate allele frequency
vcftools --gzvcf $variant_VCF --freq2 --max-alleles 2 --out ${output_path}/stats_postFilter/$output_prefix 
# mean depth of coverage per individual
vcftools --gzvcf $variant_VCF --depth --out ${output_path}/stats_postFilter/$output_prefix
# mean depth per site
vcftools --gzvcf $variant_VCF --site-mean-depth --out ${output_path}/stats_postFilter/$output_prefix
# site quality
vcftools --gzvcf $variant_VCF --site-quality --out ${output_path}/stats_postFilter/$output_prefix
# proportion of missing data per individual
vcftools --gzvcf $variant_VCF --missing-indv --out ${output_path}/stats_postFilter/$output_prefix
# proportion of missing data per site
vcftools --gzvcf $variant_VCF --missing-site --out ${output_path}/stats_postFilter/$output_prefix
# calculate hwe
vcftools --gzvcf $variant_VCF --hardy --out ${output_path}/stats_postFilter/$output_prefix
done
