#!/bin/bash -l

#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J depth
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL


module load bioinfo-tools vcftools

vcftools --gzvcf $1 --site-mean-depth --out $2
