#!/usr/bin/bash

echo "generating windowed pairwise rho"
export rho_z=../../../data/rho/ldhat_rho/z
for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
Rscript rho_slidingWin.R ${rho_z}/par/${scaffold}.per.site.txt \
../../../data/bed/par_scaf.bed 200000 50000 ${rho_z}
done

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
Rscript rho_slidingWin.R ${rho_z}/nonpar/${scaffold}.per.site.txt \
../../../data/bed/nonpar_scaf.bed 200000 50000 ${rho_z}
done

echo -e 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_z}/chrZ.map.length.txt
echo -e 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_z}/chrZ.200kbBin.50.kbStep.rates.out

for scaffold in 26 #54 35 36 
do
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/par/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/par/superscaffold${scaffold}.*.200kbBin.50.kbStep.rates.out | cat - >> ${rho_z}/chrZ.200kbBin.50.kbStep.rates.out
done

for scaffold in 36 62 67 69-1 93 63 88 83 92
do
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/nonpar/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/nonpar/superscaffold${scaffold}.*.200kbBin.50.kbStep.rates.out | cat - >> ${rho_z}/chrZ.200kbBin.50.kbStep.rates.out
done

echo "Add chromosome Z coordinates"

python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.map.length.txt rho > temp && mv temp ${rho_z}/chrZ.map.length.txt
python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.200kbBin.50.kbStep.rates.out rho > temp && mv temp ${rho_z}/chrZ.200kbBin.50.kbStep.rates.out
