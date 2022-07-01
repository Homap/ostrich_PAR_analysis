#!/bin/bash -l

#SBATCH -A snic2022-22-149
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J depth
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL


module load bioinfo-tools samtools/1.14

samtools depth -a -H ../../data/bam/P1878_107.bam ../../data/bam/P1878_108.bam \
../../data/bam/P1878_109.bam ../../data/bam/P1878_110.bam ../../data/bam/P1878_111.bam \
../../data/bam/P1878_112.bam ../../data/bam/P1878_113.bam ../../data/bam/P1878_114.bam \
../../data/bam/P1878_115.bam ../../data/bam/P1878_116.bam -o ../../data/coverage/coverage_per_site.txt
