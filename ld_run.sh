#!/bin/bash -l

#SBATCH -A snic2020-15-128 
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH -J LD
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools

# Since LD analysis is computationally heavy, we separate the VCF by scaffold
mkdir -p ../data/LD
mkdir -p ../data/LD/scaf.split

for species in black blue red  
do
echo $species
for scaffold in $(cat ../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed \
| grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
(vcftools --gzvcf ../data/vcf/${species}.all.snp.filtered.vcf.gz --positions \
<(grep $scaffold ../data/vcf/LD_vcf/${species}.Z.pos) --recode --out \
../data/LD/scaf.split/${species}.${scaffold})
done
done

# Estimate LD decay per scaffold
for species in black blue red
do
echo $species
for scaffold in $(cat ../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed \
| grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
./PopLDdecay/bin/PopLDdecay -InVCF \
../data/LD/scaf.split/${species}.${scaffold}.recode.vcf \
-OutStat ../data/LD/scaf.split/${species}.${scaffold}.LDdecay
perl PopLDdecay/bin/Plot_OnePop.pl -inFile \
../data/LD/scaf.split/${species}.${scaffold}.LDdecay.stat.gz \
-output ../data/LD/scaf.split/${species}.${scaffold}.LDdecay
done
done

# Calculate pairwise LD between SNPs
bash pairwise_ld_run.sh

