#!/bin/bash -l

#SBATCH -A snic2020-5-639
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J filter.A.DP
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07

for species in blue red #black
do
echo $species
vcftools --gzvcf ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz \
--bed ../data/bed/${species}.A.filtered.sites.coverage.bed \
--recode --recode-INFO-all --stdout | gzip -c \
> ../data/vcf/${species}.A.PASS.biallelic.nomissing.DP.vcf.gz
done
