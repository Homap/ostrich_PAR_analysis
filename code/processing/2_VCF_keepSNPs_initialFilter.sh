#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J vcf_keep_snps
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools/0.1.16 htslib/1.14

export a_vcf=../../data/vcf/a_vcf/black.A.vcf.gz
export par_vcf=../../data/vcf/par_vcf/black.PAR.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/black.nonPAR.vcf.gz

for variant_VCF in ${nonpar_vcf} ${par_vcf} ${a_vcf} 
do
output_prefix=$(echo $variant_VCF | cut -f1-6 -d ".")
echo $output_prefix
vcftools --gzvcf $variant_VCF \
--remove-indels --min-alleles 2 --max-alleles 2 --remove-filtered-all \
--recode --recode-INFO-all \
--stdout | bgzip -c > ${output_prefix}.biallelic.snp.vcf.gz
done