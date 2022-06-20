#!/usr/bin/python

# -----------------------------------------------------------------------------------------------------------
# Create stat summary output directories
# -----------------------------------------------------------------------------------------------------------

mkdir -p ../../../data/rho/ldhat_stat/chr4 \
         ../../../data/rho/ldhat_stat/chr5 \
         ../../../data/rho/ldhat_stat/z/par \
         ../../../data/rho/ldhat_stat/z/nonpar 

export stat_chr4=../../../data/rho/ldhat_stat/chr4
export stat_chr5=../../../data/rho/ldhat_stat/chr5
export stat_par=../../../data/rho/ldhat_stat/z/par
export stat_nonpar=../../../data/rho/ldhat_stat/z/nonpar

export ldhat_output_chr4=../../../data/rho/ldhat_output/chr4
export ldhat_output_chr5=../../../data/rho/ldhat_output/chr5
export ldhat_output_par=../../../data/rho/ldhat_output/z/par
export ldhat_output_nonpar=../../../data/rho/ldhat_output/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# - Autosome
# -----------------------------------------------------------------------------------------------------------
# We have run MCMC chain 3 times, we take the first chain.

scaffold=superscaffold11
echo $scaffold

for window in ${chr4_dir}/${scaffold}.*.sites.txt
do
    win_num=$(echo $window | cut -f8 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} ldhat_stat.sh ${ldhat_output_chr4}/${scaffold}.${win_num}.1.rates.txt ${stat_chr4}/${scaffold}.${win_num}.
done

scaffold=superscaffold8
echo $scaffold

for window in ${chr5_dir}/${scaffold}.*.sites.txt
do
    win_num=$(echo $window | cut -f8 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} ldhat_stat.sh ${ldhat_output_chr5}/${scaffold}.${win_num}.1.rates.txt ${stat_chr5}/${scaffold}.${win_num}.
done

# -----------------------------------------------------------------------------------------------------------
# - PAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for window in ${par_dir}/${scaffold}.*.sites.txt
    do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} ldhat_stat.sh ${ldhat_output_par}/${scaffold}.${win_num}.1.rates.txt ${stat_par}/${scaffold}.${win_num}.
    done
done

# -----------------------------------------------------------------------------------------------------------
# - nonPAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for window in ${nonpar_dir}/${scaffold}.*.sites.txt
    do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} ldhat_stat.sh ${ldhat_output_nonpar}/${scaffold}.${win_num}.1.rates.txt ${stat_nonpar}/${scaffold}.${win_num}.
    done
done



