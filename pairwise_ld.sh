#!/bin/bash -l

#SBATCH -A snic2020-15-128 
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH -J LD
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools 

vcf=$1
output=$2

./PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat $output -OutType 8
