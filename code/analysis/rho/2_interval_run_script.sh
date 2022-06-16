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

export ldhat_output=../../../data/rho/ldhat_output/chr4
export ldhat_output=../../../data/rho/ldhat_output/chr5
export ldhat_output=../../../data/rho/ldhat_output/z/par
export ldhat_output=../../../data/rho/ldhat_output/z/nonpar

# -----------------------------------------------------------------------------------------------------------
# Set path to ldhat input directories:
# -----------------------------------------------------------------------------------------------------------

export chr4_dir=../../../data/rho/ldhat_input/chr4
export chr5_dir=../../../data/rho/ldhat_input/chr4
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
for window in ${chr4_dir}/${scaffold}.*
do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} --output ${ldhat_output}/${scaffold}.${win_num} \
    ldhat_interval.sh ${scaffold}.2000.200.${win_num}.sites.txt ${scaffold}.2000.200.${win_num}.locs.txt $LUT_20 ${ldhat_output}/${scaffold}.${win_num}.
done

scaffold=superscaffold8
echo $scaffold
for window in ${chr4_dir}/${scaffold}.*
do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} --output ${ldhat_output}/${scaffold}.${win_num} \
    ldhat_interval.sh ${scaffold}.2000.200.${win_num}.sites.txt ${scaffold}.2000.200.${win_num}.locs.txt $LUT_20 ${ldhat_output}/${scaffold}.${win_num}.
done

# -----------------------------------------------------------------------------------------------------------
# - PAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for window in ${par_dir}/${scaffold}.*
    do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} --output ${ldhat_output}/${scaffold}.${win_num} \
    ldhat_interval.sh ${scaffold}.2000.200.${win_num}.sites.txt ${scaffold}.2000.200.${win_num}.locs.txt $LUT ${ldhat_output}/${scaffold}.${win_num}.
    done
done

# -----------------------------------------------------------------------------------------------------------
# - nonPAR
# -----------------------------------------------------------------------------------------------------------

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    for window in ${nonpar_dir}/${scaffold}.*
    do
    win_num=$(echo $window | cut -f9 -d "/" | cut -f4 -d ".")
    echo $win_num
    sbatch --job-name ${scaffold}.${win_num} --output ${ldhat_output}/${scaffold}.${win_num} \
    ldhat_interval.sh ${scaffold}.2000.200.${win_num}.sites.txt ${scaffold}.2000.200.${win_num}.locs.txt $LUT ${ldhat_output}/${scaffold}.${win_num}.
    done
done
