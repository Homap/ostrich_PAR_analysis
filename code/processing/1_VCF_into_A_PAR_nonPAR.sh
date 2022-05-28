#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J VCF.into.A.PAR.nonPAR
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools/0.1.16 htslib/1.14

export vcf_in=../../data/vcf/original_vcf/black.all.vcf.gz

mkdir -p ../../data/vcf/a_vcf
mkdir -p ../../data/vcf/par_vcf
mkdir -p ../../data/vcf/nonpar_vcf

export a_vcf=../../data/vcf/a_vcf/black.A.vcf.gz
export par_vcf=../../data/vcf/par_vcf/black.PAR.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/black.nonPAR.vcf.gz

vcftools --gzvcf ${vcf_in} --bed ../../data/bed/autosomes.minusC.bed --recode --recode-INFO-all --stdout | bgzip -c > ${a_vcf}
vcftools --gzvcf ${vcf_in} --bed ../../data/bed/par_scaf.bed --recode --recode-INFO-all --stdout | bgzip -c > ${par_vcf}
vcftools --gzvcf ${vcf_in} --bed ../../data/bed/nonpar_scaf.bed --recode --recode-INFO-all --stdout | bgzip -c > ${nonpar_vcf}



