#!/bin/bash -l

#SBATCH -A snic2020-15-128
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 20:00:00
#SBATCH -J filter.A.hwe
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
# Select only SNPs

module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07

for species in blue red
do
echo $species
vcftools --gzvcf ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz \
--bed ../data/vcf/hwe/${species}.A.hwe.filtered.bed \
--recode --recode-INFO-all --stdout | gzip -c > \
../data/vcf/${species}.A.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz
done
