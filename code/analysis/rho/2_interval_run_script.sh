#!/usr/bin/bash

# -----------------------------------------------------------------------------------------------------------
# Create output directories:
# -----------------------------------------------------------------------------------------------------------

mkdir -p ../../../data/rho/ldhat_output/chr4 \
         ../../../data/rho/ldhat_output/chr5 \
         ../../../data/rho/ldhat_output/z/par \
         ../../../data/rho/ldhat_output/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# Set path to ldhat output directories:
# -----------------------------------------------------------------------------------------------------------

export ldhat_output_chr4=../../../data/rho/ldhat_output/chr4
export ldhat_output_chr5=../../../data/rho/ldhat_output/chr5
export ldhat_output_par=../../../data/rho/ldhat_output/z/par
export ldhat_output_nonpar=../../../data/rho/ldhat_output/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# Set path to ldhat input directories:
# -----------------------------------------------------------------------------------------------------------

export chr4_dir=../../../data/rho/ldhat_input/chr4
export chr5_dir=../../../data/rho/ldhat_input/chr5
export par_dir=../../../data/rho/ldhat_input/z/par
export nonpar_dir=../../../data/rho/ldhat_input/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# Set path to likelihood look up tables:
# -----------------------------------------------------------------------------------------------------------

export LUT_20=../../../data/rho/ldhat_input/lk_LUT/auto_PARnew_lk.txt
export LUT_10=../../../data/rho/ldhat_input/lk_LUT/nonPARnew_lk.txt

# -----------------------------------------------------------------------------------------------------------
# - Autosome
# -----------------------------------------------------------------------------------------------------------

scaffold=superscaffold11
echo $scaffold
for mcmc in 1 2 3
do
echo $mcmc
for window in ${chr4_dir}/${scaffold}.*.sites.txt
do
    win_num=$(echo $window | cut -f8 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num}.${mcmc} --output ${ldhat_output_chr4}/${scaffold}.${win_num}.${mcmc} \
    ldhat_interval.sh ${chr4_dir}/${scaffold}.1000.200.${win_num}.sites.txt ${chr4_dir}/${scaffold}.1000.200.${win_num}.locs.txt $LUT_20 ${ldhat_output_chr4}/${scaffold}.${win_num}.${mcmc}.
done
done

scaffold=superscaffold8
echo $scaffold
for mcmc in 1 2 3
do
echo $mcmc
for window in ${chr5_dir}/${scaffold}.*.sites.txt
do
    win_num=$(echo $window | cut -f8 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num}.${mcmc} --output ${ldhat_output_chr5}/${scaffold}.${win_num}.${mcmc} \
    ldhat_interval.sh ${chr5_dir}/${scaffold}.1000.200.${win_num}.sites.txt ${chr5_dir}/${scaffold}.1000.200.${win_num}.locs.txt $LUT_20 ${ldhat_output_chr5}/${scaffold}.${win_num}.${mcmc}.
done
done

# -----------------------------------------------------------------------------------------------------------
# - PAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for mcmc in 1 2 3
    do
    echo $mcmc
        for window in ${par_dir}/${scaffold}.*.sites.txt
        do
        win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
        echo $win_num
        sbatch --job-name ${scaffold}.${win_num}.${mcmc} --output ${ldhat_output_par}/${scaffold}.${win_num}.${mcmc} \
        ldhat_interval.sh ${par_dir}/${scaffold}.1000.200.${win_num}.sites.txt \
        ${par_dir}/${scaffold}.1000.200.${win_num}.locs.txt $LUT_20 ${ldhat_output_par}/${scaffold}.${win_num}.${mcmc}.
        done
    done
done

# -----------------------------------------------------------------------------------------------------------
# - nonPAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for mcmc in 1 2 3
    do
    echo $mcmc
        for window in ${nonpar_dir}/${scaffold}.*.sites.txt
        do
        win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
        echo $win_num
        sbatch --job-name ${scaffold}.${win_num}.${mcmc} --output ${ldhat_output_nonpar}/${scaffold}.${win_num}.${mcmc} \
        ldhat_interval.sh ${nonpar_dir}/${scaffold}.1000.200.${win_num}.sites.txt \
        ${nonpar_dir}/${scaffold}.1000.200.${win_num}.locs.txt $LUT_10 ${ldhat_output_nonpar}/${scaffold}.${win_num}.${mcmc}.
        done
    done
done

# Rerunning failed jobs

sbatch --job-name superscaffold11.175.3 --output ${ldhat_output_chr4}/superscaffold11.175.3 \
    ldhat_interval.sh ${chr4_dir}/superscaffold11.1000.200.175.sites.txt ${chr4_dir}/superscaffold11.1000.200.175.locs.txt $LUT_20 ${ldhat_output_chr4}/superscaffold11.175.3.