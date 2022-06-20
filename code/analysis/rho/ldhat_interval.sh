#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 03:00:00
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

sites=$1
locs=$2
lk=$3
out_prefix=$4

./LDhat/interval -seq $sites -loc $locs -lk $lk -prefix $out_prefix -its 25000000 -bpen 5 -samp 5000
