#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00:00
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

vcf=$1
output=$2

./PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat $output -OutType 8
