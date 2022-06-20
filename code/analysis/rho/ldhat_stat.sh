#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

rates=$1
locs=$2
out_prefix=$3

./LDhat/stat -input $rates -burn 1000 -loc $locs -prefix $out_prefix