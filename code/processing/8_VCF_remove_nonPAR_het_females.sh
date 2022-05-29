#!/bin/bash -l

module load bioinfo-tools bcftools/1.14 tabix/0.2.6 vcftools/0.1.16

python check_nonPAR_genotype.py ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.vcf.gz > ../../data/bed/nonPAR_female_heterozygous.bed

vcftools --gzvcf ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.vcf.gz \
--exclude-bed ../../data/bed/nonPAR_female_heterozygous.bed \
--recode --recode-INFO-all --stdout | bgzip -c > ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.vcf.gz

bcftools view -H ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.vcf.gz | wc -l \
    > ../../data/vcf/nonpar_vcf/stats_postFilter/nonpar.repeatmasked.nofemalehet.allCounts