#!/bin/bash -l
#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J ldhat
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

sites=$1
locs=$2
lk=$3
out_prefix=$4


./interval -seq $sites -loc $locs -lk $lk -prefix $out_prefix -its 10000000 -bpen 5 -samp 2000
