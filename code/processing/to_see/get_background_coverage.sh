#!/bin/bash -l

#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J filter.A.DP
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

all_depth=$1
all_depth_out=$2

python nonvariant_filtering.py $all_depth > $all_depth_out
