#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J ld_window
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

# Chromosome Z -------------------------------------------------------------------------------------------
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

echo "Concatenate all scaffolds"
awk 'NR==1' ${ld_window_z}/superscaffold26.pairwise.LD.05-500.200kbBin.50.kbStep.out > ${ld_window_z}/header
cat ${ld_window_z}/*.pairwise.LD.05-500.200kbBin.50.kbStep.out | grep -v "Dprime" | cat ${ld_window_z}/header - >> ${ld_window_z}/Z.scaffolds.LD.05-500.200kbBin.50.kbStep.out

echo "converting scaffold coordinates into chromosome coordinates"
python ../../processing/scaffold_to_chr.py ${ld_window_z}/Z.scaffolds.LD.05-500.200kbBin.50.kbStep.out LD > ${ld_window_z}/Z.LD.05-500.200kbBin.50.tab

echo "windowed pairwise for chromosome 4"
# Chromosome 4 -------------------------------------------------------------------------------------------
export ld_window_autosome=../../../data/ld/ld_window/autosome

for scaffold in $(cat ../../../data/lastz/gg_chr4_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
zcat ${ld_window_autosome}/chr4/${scaffold}.pairwise.LD.gz | \
grep -v "^#" | awk '($9>=500 && $9<=50000)' > \
${ld_window_autosome}/chr4/${scaffold}.pairwise.LD.05-500
Rscript LD_slidingWin.R ${ld_window_autosome}/chr4/${scaffold}.pairwise.LD.05-500 \
../../../data/lastz/gg_chr4_ostrich.bed 200000 50000
done

echo "windowed pairwise for chromosome 5"
# Chromosome 5 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr5_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
zcat ${ld_window_autosome}/chr5/${scaffold}.pairwise.LD.gz | \
grep -v "^#" | awk '($9>=500 && $9<=50000)' > \
${ld_window_autosome}/chr5/${scaffold}.pairwise.LD.05-500
Rscript LD_slidingWin.R ${ld_window_autosome}/chr5/${scaffold}.pairwise.LD.05-500 \
../../../data/lastz/gg_chr5_ostrich.bed 200000 50000
done

echo "concatenating pairwise LD for chromosome 4 and 5"
awk 'NR==1' ${ld_window_autosome}/chr4/scaffold1404.pairwise.LD.05-500.200kbBin.50.kbStep.out > ${ld_window_autosome}/header

cat ${ld_window_autosome}/chr4/*.pairwise.LD.05-500.200kbBin.50.kbStep.out >> ${ld_window_autosome}/chr4/chr4.pairwise.LD.05-500.200kbBin.50.kbStep.out
cat ${ld_window_autosome}/chr5/*.pairwise.LD.05-500.200kbBin.50.kbStep.out >> ${ld_window_autosome}/chr5/chr5.pairwise.LD.05-500.200kbBin.50.kbStep.out

cat ${ld_window_autosome}/chr4/chr4.pairwise.LD.05-500.200kbBin.50.kbStep.out ${ld_window_autosome}/chr5/chr5.pairwise.LD.05-500.200kbBin.50.kbStep.out >> ${ld_window_autosome}/chr4_5_pairwise_LD.tab

grep -v "Dprime" ${ld_window_autosome}/chr4_5_pairwise_LD.tab | cat ${ld_window_autosome}/header - > temp && mv temp ${ld_window_autosome}/chr4_5_pairwise_LD.tab
