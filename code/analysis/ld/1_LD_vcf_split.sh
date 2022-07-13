#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J VCF.into.A.PAR.nonPAR
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools/0.1.16 htslib/1.14

# Separating autosomal scaffolds

mkdir -p ../../../data/vcf/ld_vcf/a_vcf/chr4
mkdir -p ../../../data/vcf/ld_vcf/a_vcf/chr5

export a_vcf=../../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz

# Chromosome 4 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr4_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
vcftools --gzvcf ${a_vcf} --chr $scaffold --recode --recode-INFO-all --stdout | \
bgzip -c > ../../../data/vcf/ld_vcf/a_vcf/chr4/${scaffold}.chr4.vcf.gz
done

# Chromosome 5 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr5_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
vcftools --gzvcf ${a_vcf} --chr $scaffold \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/a_vcf/chr5/${scaffold}.chr5.vcf.gz
done

# Separating Z scaffolds ---------------------------------------------------------------------------------
mkdir -p ../../../data/vcf/ld_vcf/z_vcf

export z_vcf=../../../data/vcf/z_vcf/z_vcf.gz

# chromosome Z
for scaffold in $(cat ../../../data/bed/z_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
vcftools --gzvcf ${z_vcf} --chr $scaffold \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/${scaffold}.chrz.vcf.gz
done

# PAR : This is the same as the scaffolds in Z except for superscaffold36
vcftools --gzvcf ${z_vcf} --chr superscaffold36 --from-bp 3524263 --to-bp 9394175 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.par.vcf.gz
vcftools --gzvcf ${z_vcf} --chr superscaffold36 --from-bp 1 --to-bp 3516673 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.nonpar.vcf.gz

# (100 Kb and 500 Kb)
# PAR farthest away
vcftools --gzvcf ${z_vcf} --chr superscaffold26 --from-bp 1 --to-bp 500000 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/supersaffold26.500Kb.vcf.gz

# PAR mid region
vcftools --gzvcf ${z_vcf} --chr superscaffold54 --from-bp 15879243 --to-bp 16379243 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/supersaffold54.500Kb.vcf.gz

# PAR boundary in PAR
vcftools --gzvcf ${z_vcf} --chr superscaffold36 --from-bp 3524263 --to-bp 4024263 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/supersaffold36.500Kb.vcf.gz

# PAR boundary spanning to nonPAR
vcftools --gzvcf ${z_vcf} --chr superscaffold36 --from-bp 3474263 --to-bp 3574263 \
--recode --recode-INFO-all --stdout | bgzip -c > ../../../data/vcf/ld_vcf/z_vcf/supersaffold36.100Kb.spanning.boundary.vcf.gz
