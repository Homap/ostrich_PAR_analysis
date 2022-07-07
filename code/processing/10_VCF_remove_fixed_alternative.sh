#!/usr/bin/bash

module load bioinfo-tools vcftools/0.1.16 htslib/1.14 bcftools/1.14

export par_vcf=../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.vcf.gz
export nonpar_vcf=../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.vcf.gz
export a_vcf=../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.vcf.gz

echo "PAR"
vcftools --gzvcf ${par_vcf} --counts --out ../../data/vcf/par_vcf/PAR
python fixed_sites.py ../../data/vcf/par_vcf/PAR.frq.count > ../../data/vcf/par_vcf/PAR.fixed.alternative.txt 
vcftools --gzvcf ${par_vcf} --exclude-positions ../../data/vcf/par_vcf/PAR.fixed.alternative.txt --recode --stdout | gzip -c > ../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.vcf.gz
bcftools sort ../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.vcf.gz -o ../../data/vcf/par_vcf/par_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz

echo "nonPAR"
vcftools --gzvcf ${nonpar_vcf} --counts --out ../../data/vcf/nonpar_vcf/nonPAR  
python fixed_sites.py ../../data/vcf/nonpar_vcf/nonPAR.frq.count > ../../data/vcf/nonpar_vcf/nonPAR.fixed.alternative.txt
vcftools --gzvcf ${nonpar_vcf} --exclude-positions ../../data/vcf/nonpar_vcf/nonPAR.fixed.alternative.txt --recode --stdout | gzip -c > ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.vcf.gz
bcftools sort ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.vcf.gz -o ../../data/vcf/nonpar_vcf/nonpar_vcf.filtered.repeatmasked.nofemalehet.hwe.snps.sorted.vcf.gz

echo "autosome"
vcftools --gzvcf ${a_vcf} --counts --out ../../data/vcf/a_vcf/a
python fixed_sites.py ../../data/vcf/a_vcf/a.frq.count > ../../data/vcf/a_vcf/a.fixed.alternative.txt
vcftools --gzvcf ${a_vcf} --exclude-positions ../../data/vcf/a_vcf/a.fixed.alternative.txt --recode --stdout | gzip -c > ../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.vcf.gz
bcftools sort ../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.vcf.gz -o ../../data/vcf/a_vcf/a_vcf.filtered.repeatmasked.hwe.snps.sorted.vcf.gz