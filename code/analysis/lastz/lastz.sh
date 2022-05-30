#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3-00:00:00
#SBATCH -J LASTZ
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

ref=$1
que=$2
out=$3

module load bioinfo-tools lastz/1.04.00

lastz_32 $ref[multiple] $que M=254 K=4500 L=3000 Y=15000 C=2 T=2 --matchcount=10000 \
--ambiguous=iupac \
--format=general:name1,start1,end1,strand1,name2,start2+,end2+,strand2,score > $out
