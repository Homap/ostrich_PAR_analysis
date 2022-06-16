#!/usr/bin/bash

# -----------------------------------------------------------------------------------------------------------
# Create MCMC summary output directories
# -----------------------------------------------------------------------------------------------------------

mkdir -p ../../../data/rho/ldhat_mcmc/chr4 \
         ../../../data/rho/ldhat_mcmc/chr5 \
         ../../../data/rho/ldhat_mcmc/z/par \
         ../../../data/rho/ldhat_mcmc/z/nonpar 

# -----------------------------------------------------------------------------------------------------------
# Set path to ldhat input directories:
# -----------------------------------------------------------------------------------------------------------

export mcmc_chr4=../../../data/rho/ldhat_mcmc/chr4
export mcmc_chr5=../../../data/rho/ldhat_mcmc/chr5
export mcmc_par=../../../data/rho/ldhat_mcmc/z/par
export mcmc_nonpar=../../../data/rho/ldhat_mcmc/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# Set path to ldhat input and output directories
# -----------------------------------------------------------------------------------------------------------

export ldhat_output_chr4=../../../data/rho/ldhat_output/chr4
export ldhat_output_chr5=../../../data/rho/ldhat_output/chr5
export ldhat_output_par=../../../data/rho/ldhat_output/z/par
export ldhat_output_nonpar=../../../data/rho/ldhat_output/z/nonpar

export chr4_dir=../../../data/rho/ldhat_input/chr4
export chr5_dir=../../../data/rho/ldhat_input/chr5
export par_dir=../../../data/rho/ldhat_input/z/par
export nonpar_dir=../../../data/rho/ldhat_input/z/nonpar

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
    grep 'LK' ${ldhat_output_chr4}/${scaffold}.${win_num}.${mcmc} | \
    awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' \
    > ${mcmc_chr4}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
    grep 'LK' ${ldhat_output_chr5}/${scaffold}.${win_num}.${mcmc} | \
    awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' \
    > ${mcmc_chr5}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
        grep 'LK' ${ldhat_output_par}/${scaffold}.${win_num}.${mcmc} | \
        awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' \
        > ${mcmc_par}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
        grep 'LK' ${ldhat_output_nonpar}/${scaffold}.${win_num}.${mcmc} | \
        awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' \
        > ${mcmc_nonpar}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
        done
    done
done