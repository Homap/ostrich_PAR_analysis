#!/bin/bash -l

#SBATCH -A snic2019-3-17 
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J filtering
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07

for species in black.male blue.male red.male black.female blue.female red.female
do
echo $species
echo "number of autosomal SNPs before filtering"
grep -v "#" ../data/vcf/${species}.vcf | wc
vcftools --vcf ../data/vcf/${species}.vcf --bed ../data/bed/autosomes_with_SNPs.bed \
--recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.A.vcf.gz
echo "number of autosomal SNPs"
zgrep -v "#" ../data/vcf/${species}.A.vcf.gz | wc

vcftools --gzvcf ../data/vcf/${species}.A.vcf.gz --remove-filtered-all \
--recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.A.PASS.vcf.gz
echo "number of autosomal SNPs that has PASS flag"
zgrep -v "#" ../data/vcf/${species}.A.PASS.vcf.gz | wc

vcftools --gzvcf ../data/vcf/${species}.A.PASS.vcf.gz \
--min-alleles 2 --max-alleles 2 --recode --recode-INFO-all \
--stdout | gzip -c > ../data/vcf/${species}.A.PASS.biallelic.vcf.gz
echo "number of autosomal SNPs that has PASS flag and are biallelic"
zgrep -v "#" ../data/vcf/${species}.A.PASS.biallelic.vcf.gz | wc

vcftools --gzvcf ../data/vcf/${species}.A.PASS.biallelic.vcf.gz --max-missing 1.0 \
--recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz
echo "number of autosomal SNPs that has PASS flag and are biallelic and have no missing data"
zgrep -v "#" ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz | wc
done

# Calculate HWE for each site
mkdir -p ../data/vcf/hwe
for species in black.male blue.male red.male black.female blue.female red.female
do
echo $species
vcftools --gzvcf ../data/vcf/${species}.A.PASS.biallelic.nomissing.vcf.gz --hardy \
--out ../data/vcf/hwe/${species}.A
done
