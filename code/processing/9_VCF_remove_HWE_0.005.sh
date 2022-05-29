#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J VCFHWE
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools/1.14 tabix/0.2.6 vcftools/0.1.16

export a_vcf=../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.vcf.gz
export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.vcf.gz

for variant_VCF in ${par_vcf} ${a_vcf} 
do
    output_prefix=$(echo $variant_VCF | cut -f6 -d "/" | cut -f1 -d "_")
    output_path=$(echo $variant_VCF | cut -f1-5 -d "/")
    echo $variant_VCF $output_prefix $output_path
    vcftools --gzvcf ${variant_VCF} --hwe 0.005 \
    --recode --recode-INFO-all --stdout | bgzip -c > ${output_path}/${output_prefix}_vcf.filtered.repeatmasked.hwe.vcf.gz
    bcftools view -H ${output_path}/${output_prefix}_vcf.filtered.repeatmasked.hwe.vcf.gz | wc -l \
    > ${output_path}/stats_postFilter/${output_prefix}.repeatmasked.hwe.allCounts
done

# In nonPAR, the homozygosity in females can cause HWE deviation. For this reason, we only remove 
# deviations from HWE based on male VCF.

export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.vcf.gz

# Create a new VCF for nonpar with male only
vcftools --gzvcf ${nonpar_vcf} --keep ../../data/samples/black_male.txt --recode --recode-INFO-all --stdout | bgzip -c > \
../../data/vcf/nonpar_vcf/nonpar_male_vcf.gz

# Calculate HWE for the male nonpar vcf
vcftools --gzvcf ../../data/vcf/nonpar_vcf/nonpar_male_vcf.gz --hardy --out ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male

# Create a bed file from the nonpar_male.hwe
awk '$6 < 0.005' ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male.hwe | \
awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $1"\t"$2-1"\t"$2}' \
> ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male.hwe.bed

vcftools --gzvcf ${nonpar_vcf} \
--exclude-bed ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male.hwe.bed \
--recode --recode-INFO-all --stdout | bgzip -c > \
../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
bcftools view -H ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz | wc -l \
> ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar.repeatmasked.nofemalehet.hwe.allCounts

# Remove intermediate files
rm -f ../../data/vcf/nonpar_vcf/nonpar_male_vcf.gz
rm -f ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male.hwe
rm -f ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar_male.hwe.bed
