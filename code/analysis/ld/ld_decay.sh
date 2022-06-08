#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

vcf=$1
output=$2
input_bin=$3
output_bin=$4

./PopLDdecay/bin/PopLDdecay -InVCF $vcf -OutStat $output

perl PopLDdecay/bin/Plot_OnePop.pl -inFile $input_bin -output $output_bin
