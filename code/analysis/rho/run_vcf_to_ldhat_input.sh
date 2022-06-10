# - autosome
# For autosomes, we are only interested to have a mean of rho to compare with the PAR. To reduce the computation time, only the longest scaffolds from each
# chromosome 4 and 5 is chosen for the rho analysis.

# Make output directories
mkdir -p ../../../data/rho/ldhat_input/chr4 \
         ../../../data/rho/ldhat_input/chr5 \
         ../../../data/rho/ldhat_input/z


# - chr4
export chr4_dir=../../../data/vcf/ld_vcf/a_vcf/chr4
python vcf_to_ldhat_out.py ${chr4_dir}/superscaffold11.chr4.vcf.gz ../../../data/lastz/gg_chr4_ostrich.bed 2000 500 superscaffold11 1000 ../../../data/rho


# - chr5
export chr5_dir=../../../data/vcf/ld_vcf/a_vcf/chr5
python vcf_to_ldhat_out.py ${chr5_dir}/superscaffold8.chr5.vcf.gz ../../../data/lastz/gg_chr5_ostrich.bed 2000 500 superscaffold8 1000


- PAR
    ```
    export z_dir=../../../data/vcf/z_vcf
    for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
    do
        echo $scaffold
        python vcf_to_ldhat_out.py vcf_rho/z_vcf.gz ../../../data/bed/par_scaf.bed 2000 500 $scaffold 1000
    done
    ```

- nonPAR
    ```
    export z_dir=../../../data/vcf/z_vcf
    for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
    do
        echo $scaffold
        python vcf_to_ldhat_out.py vcf_rho/black.PAR.hwe.filtered.vcf.gz ../../../data/bed/nonpar_scaf.bed 2000 500 $scaffold 1000
    done 
    ```