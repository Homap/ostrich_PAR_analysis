#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J ld_window
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

echo "generating windowed pairwise LD"
export ld_window_z=../../../data/ld/ld_window/z
for scaffold in $(cat ../../../data/bed/z_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
zcat ${ld_window_z}/${scaffold}.pairwise.LD.gz | \
grep -v "^#" | awk '($9>=500 && $9<=50000)' > \
${ld_window_z}/${scaffold}.pairwise.LD.05-500
Rscript LD_slidingWin.R ${ld_window_z}/${scaffold}.pairwise.LD.05-500 \
../../../data/bed/z_scaf.bed 200000 50000
done

echo "Concatenate all scaffolds for each subspecies LD"
cat ${ld_window_z}/*.pairwise.LD.05-500.200kbBin.50.kbStep.out >> ${ld_window_z}/Z.scaffolds.LD.05-500.200kbBin.50.kbStep.out

echo "converting scaffold coordinates into chromosome coordinates"
python ../../processing/scaffold_to_chr.py ${ld_window_z}/Z.scaffolds.LD.05-500.200kbBin.50.kbStep.out > ${ld_window_z}/Z.LD.05-500.200kbBin.50.tab
