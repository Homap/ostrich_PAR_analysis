#!/usr/bin/bash

export rho_z=../../../data/rho/ldhat_rho/z

echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_z}/chrZ.map.length.txt
echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_z}/chrZ.per.site.txt

for scaffold in 26 54 35 36 
do
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/par/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/par/superscaffold${scaffold}.per.site.txt | cat - >> ${rho_z}/chrZ.per.site.txt
done

for scaffold in 36 62 67 69-1 93 63 88 83 92
do
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/nonpar/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
    awk 'NR>1' ../../../data/rho/ldhat_rho/z/nonpar/superscaffold${scaffold}.per.site.txt | cat - >> ${rho_z}/chrZ.per.site.txt
done

sed -i '/^$/d' ${rho_z}/chrZ.per.site.txt 

echo "Add chromosome Z coordinates"

python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.map.length.txt rho > temp && mv temp ${rho_z}/chrZ.map.length.txt
python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.per.site.txt rho > temp && mv temp ${rho_z}/chrZ.per.site.txt
