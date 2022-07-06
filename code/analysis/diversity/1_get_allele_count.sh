#!/usr/bin/bash

module load bioinfo-tools vcftools/0.1.16 htslib/1.14 bcftools/1.14

# make output directory
mkdir -p ../../../data/diversity/allele_count
mkdir -p ../../../data/diversity/allele_count/chr4
mkdir -p ../../../data/diversity/allele_count/chr5
mkdir -p ../../../data/diversity/allele_count/par
mkdir -p ../../../data/diversity/allele_count/nonpar
mkdir -p ../../../data/diversity/allele_count/z

export ac_out=../../../data/diversity/allele_count

# -------------------------------------------------------------------------------------------------------------
export par_vcf=../../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.vcf.gz
export nonpar_vcf=../../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
export z_vcf=../../../data/vcf/z_vcf/sorted.z.vcf.gz
# -------------------------------------------------------------------------------------------------------------
export a_vcf=../../../data/vcf/ld_vcf/a_vcf
# -------------------------------------------------------------------------------------------------------------
echo "Autosome"
# -------------------------------------------------------------------------------------------------------------
echo "Chr4"
for scaffold_vcf in ${a_vcf}/chr4/*chr4.vcf.gz
do
echo $scaffold_vcf
scaffold=$(echo $scaffold_vcf | cut -f9 -d "/" | cut -f1 -d ".")
bcftools sort $scaffold_vcf -o ${a_vcf}/chr4/${scaffold}.chr4.sorted
bgzip ${a_vcf}/chr4/${scaffold}.chr4.sorted
vcftools --gzvcf ${a_vcf}/chr4/${scaffold}.chr4.sorted.gz --counts --out ${ac_out}/chr4/${scaffold}
done
# -------------------------------------------------------------------------------------------------------------
echo "Chr5"
for scaffold_vcf in ${a_vcf}/chr5/*chr5.vcf.gz
do
echo $scaffold_vcf
scaffold=$(echo $scaffold_vcf | cut -f9 -d "/" | cut -f1 -d ".")
bcftools sort $scaffold_vcf -o ${a_vcf}/chr5/${scaffold}.chr5.sorted
bgzip ${a_vcf}/chr5/${scaffold}.chr5.sorted
vcftools --gzvcf ${a_vcf}/chr5/${scaffold}.chr5.sorted.gz --counts --out ${ac_out}/chr5/${scaffold}
done
# -------------------------------------------------------------------------------------------------------------
echo "PAR"
vcftools --gzvcf ${par_vcf} --counts --out ${ac_out}/par/PAR
echo "nonPAR"
vcftools --gzvcf ${nonpar_vcf} --counts --out ${ac_out}/nonpar/nonPAR
echo "Z"
vcftools --gzvcf ${z_vcf} --counts --out ${ac_out}/z/z