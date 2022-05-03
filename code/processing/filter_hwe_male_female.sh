#!/bin/bash -l

#SBATCH -A snic2019-3-17
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 05:00:00
#SBATCH -J filter.A.hwe
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
# Select only SNPs

module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07

for species in red.female blue.female #black.male blue.male red.male black.female blue.female red.female
do
echo $species
vcftools --gzvcf ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz \
--bed ../data/vcf/hwe/${species}.A.hwe.filtered.bed \
--recode --recode-INFO-all --stdout | gzip -c > \
../data/vcf/${species}.A.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz
done
