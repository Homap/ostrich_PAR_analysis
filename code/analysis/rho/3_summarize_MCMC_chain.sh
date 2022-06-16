#!/usr/bin/bash

## Create the Tracer input to check for the distribution of likelihood, block number and map length
export ldhat_output_chr4=../../../data/rho/ldhat_output/chr4
export ldhat_output_chr5=../../../data/rho/ldhat_output/chr5
export ldhat_output_par=../../../data/rho/ldhat_output/z/par
export ldhat_output_nonpar=../../../data/rho/ldhat_output/z/nonpar

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
    awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' > ${ldhat_output_chr4}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
    awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' > ${ldhat_output_chr5}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
        awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' > ${ldhat_output_par}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
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
        awk 'BEGIN{print "iterations""\t""likelihood""\t""blocks""\t""map_length"} {print NR"\t"$5"\t"$10"\t"$15}' > ${ldhat_output_nonpar}/${scaffold}.${win_num}.${mcmc}.chainreport.txt
        done
    done
done