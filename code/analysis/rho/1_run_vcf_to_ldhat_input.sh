#!/usr/bin/bash

module load bcftools/1.14 htslib/1.14 vcftools/0.1.16

echo "Make output directories"
mkdir -p ../../../data/rho/ldhat_input/chr4 \
         ../../../data/rho/ldhat_input/chr5 \
         ../../../data/rho/ldhat_input/z/par \
         ../../../data/rho/ldhat_input/z/nonpar

# For autosomes, we are only interested to have a mean of rho to compare with the PAR. To reduce the computation time, only the longest scaffolds from each
# chromosome 4 and 5 is chosen for the rho analysis.

snp_size=1000
snp_overlap=200
num_win=10000

echo "----------------------chr4----------------------"
export chr4_dir=../../../data/vcf/ld_vcf/a_vcf/chr4
bcftools sort ${chr4_dir}/superscaffold11.chr4.vcf.gz | bgzip -c > ${chr4_dir}/superscaffold11.chr4.sorted.vcf.gz
python vcf_to_ldhat_out.py ${chr4_dir}/superscaffold11.chr4.sorted.vcf.gz ../../../data/lastz/gg_chr4_ostrich.bed $snp_size $snp_overlap superscaffold11 $num_win ../../../data/rho/ldhat_input/chr4

echo "----------------------chr5----------------------"
export chr5_dir=../../../data/vcf/ld_vcf/a_vcf/chr5
bcftools sort ${chr5_dir}/superscaffold8.chr5.vcf.gz | bgzip -c > ${chr5_dir}/superscaffold8.chr5.sorted.vcf.gz
python vcf_to_ldhat_out.py ${chr5_dir}/superscaffold8.chr5.sorted.vcf.gz ../../../data/lastz/gg_chr5_ostrich.bed $snp_size $snp_overlap superscaffold8 $num_win ../../../data/rho/ldhat_input/chr5

echo "----------------------chrZ--PAR-----------------"
export z_dir=../../../data/vcf/z_vcf
for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    python vcf_to_ldhat_out.py ${z_dir}/sorted.z.vcf.gz ../../../data/bed/par_scaf.bed $snp_size $snp_overlap $scaffold $num_win ../../../data/rho/ldhat_input/z/par
done
echo "----------------------chrZ--nonPAR--------------"
echo "Produce a male - VCF for nonPAR"

if [[ ! -f ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz ]]
then
    echo "Loading modules"
    export nonpar_vcf=../../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
    vcftools --gzvcf ${nonpar_vcf} --keep ../../../data/samples/black_male.txt --recode --recode-INFO-all --stdout | bgzip -c > \
    ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
fi

echo "Produce nonPAR windows for male nonPAR vcf"
export male_nonPAR=../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
bcftools sort ${male_nonPAR} | bgzip -c > ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.sorted.vcf.gz

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    python vcf_to_ldhat_out.py ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.sorted.vcf.gz \
     ../../../data/bed/nonpar_scaf.bed $snp_size $snp_overlap $scaffold $num_win ../../../data/rho/ldhat_input/z/nonpar
done


