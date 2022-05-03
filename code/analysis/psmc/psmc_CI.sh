#!/bin/bash -l

#SBATCH -A snic2019-3-17 
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH -J psmcCI
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools psmc

outdir=$1
split_psmcfa=$2


seq 10 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" \
    -o $outdir/round_repeat-{}.psmc $split_psmcfa | sh
