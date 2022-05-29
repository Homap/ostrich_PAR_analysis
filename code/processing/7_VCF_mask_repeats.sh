#!/bin/bash -l

#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J VCFrepeatmask
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools/1.14 tabix/0.2.6 vcftools/0.1.16

if [[ ! -f ../../data/bed/ostrich_repeats.bed ]]
then
    grep -v "#" ../../data/genome/repeatmask/Struthio_camelus.20130116.OM.fa.out.gff | \
    awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $1"\t"$4-1"\t"$5}' > ../../data/bed/ostrich_repeats.bed
fi

export a_vcf=../../data/vcf/a_vcf/a_vcf.filtered.vcf.gz
export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.vcf.gz

for variant_VCF in ${nonpar_vcf} ${par_vcf} ${a_vcf} 
do
    output_prefix=$(echo $variant_VCF | cut -f6 -d "/" | cut -f1 -d "_")
    output_path=$(echo $variant_VCF | cut -f1-5 -d "/")
    echo $variant_VCF $output_prefix $output_path
    vcftools --gzvcf ${variant_VCF} --exclude-bed  ../../data/bed/ostrich_repeats.bed \
    --recode --recode-INFO-all --stdout | bgzip -c > ${output_path}/${output_prefix}_vcf.filtered.repeatmasked.vcf.gz
    bcftools view -H ${output_path}/${output_prefix}_vcf.filtered.repeatmasked.vcf.gz | wc -l \
    > ${output_path}/stats_postFilter/${output_prefix}.repeatmasked.allCounts
done
