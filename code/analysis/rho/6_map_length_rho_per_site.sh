#!/usr/bin/bash

# -----------------------------------------------------------------------------------------------------------
# Create stat summary output directories
# -----------------------------------------------------------------------------------------------------------

mkdir -p ../../../data/rho/ldhat_rho/chr4 \
         ../../../data/rho/ldhat_rho/chr5 \
         ../../../data/rho/ldhat_rho/z/par \
         ../../../data/rho/ldhat_rho/z/nonpar 

export rho_chr4=../../../data/rho/ldhat_rho/chr4
export rho_chr5=../../../data/rho/ldhat_rho/chr5
export rho_par=../../../data/rho/ldhat_rho/z/par
export rho_nonpar=../../../data/rho/ldhat_rho/z/nonpar

export stat_chr4=../../../data/rho/ldhat_stat/chr4
export stat_chr5=../../../data/rho/ldhat_stat/chr5
export stat_par=../../../data/rho/ldhat_stat/z/par
export stat_nonpar=../../../data/rho/ldhat_stat/z/nonpar

export chr4_dir=../../../data/rho/ldhat_input/chr4
export chr5_dir=../../../data/rho/ldhat_input/chr5
export par_dir=../../../data/rho/ldhat_input/z/par
export nonpar_dir=../../../data/rho/ldhat_input/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# - Autosome
# -----------------------------------------------------------------------------------------------------------

scaffold=superscaffold11
echo $scaffold

count=$(ls ${stat_chr4}/${scaffold}.* | wc -l)
echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_chr4}/${scaffold}.map.length.txt
echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_chr4}/${scaffold}.per.site.txt

for win_num in `seq $count`
do
    echo $win_num
    python stat_rho.py ${stat_chr4}/${scaffold}.${win_num}.res.txt ${chr4_dir}/${scaffold}.1000.200.${win_num}.pos.txt \
    ${rho_chr4}/${scaffold}.per.site.txt ${rho_chr4}/${scaffold}.map.length.txt $scaffold
done

scaffold=superscaffold8
echo $scaffold

count=$(ls ${stat_chr5}/${scaffold}.* | wc -l)
echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_chr5}/${scaffold}.map.length.txt
echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_chr5}/${scaffold}.per.site.txt

for win_num in `seq $count`
do
    echo $win_num
    python stat_rho.py ${stat_chr5}/${scaffold}.${win_num}.res.txt ${chr5_dir}/${scaffold}.1000.200.${win_num}.pos.txt \
    ${rho_chr5}/${scaffold}.per.site.txt ${rho_chr5}/${scaffold}.map.length.txt $scaffold
done

# -----------------------------------------------------------------------------------------------------------
# - PAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    count=$(ls ${stat_par}/${scaffold}.* | wc -l)
    echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_par}/${scaffold}.map.length.txt
    echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_par}/${scaffold}.per.site.txt
    for win_num in `seq $count`
    do
    echo $win_num
    python stat_rho.py ${stat_par}/${scaffold}.${win_num}.res.txt ${par_dir}/${scaffold}.1000.200.${win_num}.pos.txt \
    ${rho_par}/${scaffold}.per.site.txt ${rho_par}/${scaffold}.map.length.txt $scaffold
    done
done

# -----------------------------------------------------------------------------------------------------------
# - nonPAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    count=$(ls ${stat_nonpar}/${scaffold}.* | wc -l)
    echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_nonpar}/${scaffold}.map.length.txt
    echo -e -n 'Scaffold\tLocus_start\tLocus_end\tMean_rho\tMedian\tL95\tU95' > ${rho_nonpar}/${scaffold}.per.site.txt
    for win_num in `seq $count`
    do
    echo $win_num
    python stat_rho.py ${stat_nonpar}/${scaffold}.${win_num}.res.txt ${nonpar_dir}/${scaffold}.1000.200.${win_num}.pos.txt \
    ${rho_nonpar}/${scaffold}.per.site.txt ${rho_nonpar}/${scaffold}.map.length.txt $scaffold
    done
done


