#!/usr/bin/bash

for species in black blue red
do
echo $species
for scaffold in $(cat ../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed \
| grep -v "^chrom" |cut -f1 | sort | uniq)
do 
echo $scaffold
sbatch pairwise_ld.sh ../data/LD/scaf.split/${species}.${scaffold}.recode.vcf \
../data/LD/scaf.split/${species}.${scaffold}.pairwise -OutType 8
done
done

