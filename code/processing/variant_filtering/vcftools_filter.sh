#!/bin/bash -l
#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 4:00:00
#SBATCH -J vcftools_Filter
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools tabix vcftools vcflib
# set filters
VCF_IN=$1
VCF_OUT=$2
MISS=1.0
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=$3

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN --remove-filtered-all --min-alleles 2 --max-alleles 2 \
--max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT

