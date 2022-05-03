#!/bin/bash -l

#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:30:00
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL


module load bioinfo-tools samtools

#scaffold=$1
#bam=$2
coverage_out=$1
median=$2

#samtools depth -r $scaffold $bam > $coverage_out

sort -k3n $coverage_out | awk -f median.awk > $median
