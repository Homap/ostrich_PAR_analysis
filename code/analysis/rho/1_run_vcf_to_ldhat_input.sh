#!/usr/bin/bash

module load bioinfo-tools bcftools/1.14 htslib/1.14 vcftools/0.1.16 

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
export par_vcf=../../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz
zgrep -v "#" ${par_vcf} | wc -l
for scaffold in $(cat ../../../data/bed/par_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    python vcf_to_ldhat_out.py ${par_vcf} ../../../data/bed/par_scaf.bed $snp_size $snp_overlap $scaffold $num_win ../../../data/rho/ldhat_input/z/par
done
echo "----------------------chrZ--nonPAR--------------"
echo "Produce a male - VCF for nonPAR"

if [[ ! -f ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz ]]
then
    echo "Loading modules"
    export nonpar_vcf=../../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz
    vcftools --gzvcf ${nonpar_vcf} --keep ../../../data/samples/black_male.txt --recode --recode-INFO-all --stdout | bgzip -c > \
    ../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz
fi

echo "Produce nonPAR windows for male nonPAR vcf"
export male_nonPAR=../../../data/vcf/nonpar_vcf/male_nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz
zgrep -v "#" ${male_nonPAR} | wc -l

for scaffold in $(cat ../../../data/bed/nonpar_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq) 
do
    echo $scaffold
    python vcf_to_ldhat_out.py ${male_nonPAR} \
     ../../../data/bed/nonpar_scaf.bed $snp_size $snp_overlap $scaffold $num_win ../../../data/rho/ldhat_input/z/nonpar
done


