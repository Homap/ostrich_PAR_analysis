#!/usr/bin/bash

module load bioinfo-tools bcftools/1.14 tabix/0.2.6 vcftools/0.1.16 

export a_vcf=../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz
export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz

mkdir -p ../../data/vcf/a_vcf/stats_finalFilter
mkdir -p ../../data/vcf/par_vcf/stats_finalFilter
mkdir -p ../../data/vcf/nonpar_vcf/stats_finalFilter

for variant_VCF in ${nonpar_vcf} ${par_vcf} ${a_vcf} 
do
output_prefix=$(echo $variant_VCF | cut -f6 -d "/" | cut -f1 -d "_")
output_path=$(echo $variant_VCF | cut -f1-5 -d "/")
echo $variant_VCF $output_prefix $output_path
# Count the number of $variant_VCF variants
bcftools view -H $variant_VCF | wc -l > ${output_path}/stats_finalFilter/${output_prefix}.allCounts
# calculate allele frequency
vcftools --gzvcf $variant_VCF --freq2 --max-alleles 2 --out ${output_path}/stats_finalFilter/$output_prefix 
# mean depth of coverage per individual
vcftools --gzvcf $variant_VCF --depth --out ${output_path}/stats_finalFilter/$output_prefix
# mean depth per site
vcftools --gzvcf $variant_VCF --site-mean-depth --out ${output_path}/stats_finalFilter/$output_prefix
# site quality
vcftools --gzvcf $variant_VCF --site-quality --out ${output_path}/stats_finalFilter/$output_prefix
# proportion of missing data per individual
vcftools --gzvcf $variant_VCF --missing-indv --out ${output_path}/stats_finalFilter/$output_prefix
# proportion of missing data per site
vcftools --gzvcf $variant_VCF --missing-site --out ${output_path}/stats_finalFilter/$output_prefix
done