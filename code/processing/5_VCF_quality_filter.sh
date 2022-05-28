#!/bin/bash -l

#SBATCH -A snic2022-22-149
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 10:00:00
#SBATCH -J Filtering_VCF
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools GATK/4.2.0.0 vcftools/0.1.16 htslib/1.14

reference=../../data/genome/Struthio_camelus.20130116.OM.fa

# Generate index file for VCF
echo "Generating index file"
for vcf in ../../data/vcf/a_vcf/black.A.biallelic.snp.vcf.gz ../../data/vcf/par_vcf/black.PAR.biallelic.snp.vcf.gz ../../data/vcf/nonpar_vcf/black.nonPAR.biallelic.snp.vcf.gz
do
echo $vcf
if [[ ! -f $vcf ]]
then 
    prefix=$(echo $vcf | cut -f5 -d "/")
    java -Xmx12g -jar $GATK_HOME/gatk-package-4.2.0.0-local.jar IndexFeatureFile -I $vcf
fi
done

echo "Filtering with GATK"
for vcf in ../../data/vcf/a_vcf/black.A.biallelic.snp.vcf.gz ../../data/vcf/par_vcf/black.PAR.biallelic.snp.vcf.gz ../../data/vcf/nonpar_vcf/black.nonPAR.biallelic.snp.vcf.gz
do
echo $vcf
prefix=$(echo $vcf | cut -f5 -d "/")
java -Xmx12g -jar $GATK_HOME/gatk-package-4.2.0.0-local.jar VariantFiltration \
-R $reference \
-V $vcf \
-O ../../data/vcf/${prefix}/${prefix}.filterflags.vcf.gz \
--cluster-window-size 10 --cluster-size 3 \
--filter-expression "QUAL < 30.0 " --filter-name "VeryLowQual" \
--filter-expression "QD < 2.0 " --filter-name "LowQD" \
--filter-expression "FS > 60.0" --filter-name "StrandBiasFishers" \
--filter-expression "MQ < 40.0" --filter-name "Mapping Quality" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
--filter-expression "SOR > 3.0" --filter-name "SOR" 

# perform the filtering with vcftools
echo "Removing filtered sites"
vcftools --gzvcf ../../data/vcf/${prefix}/${prefix}.filterflags.vcf.gz --remove-filtered-all \
--min-meanDP 5 --max-meanDP 70 \
--minDP 5 --maxDP 70 \
--max-missing 1.0 --recode --recode-INFO-all --stdout | bgzip -c > ../../data/vcf/${prefix}/${prefix}.filtered.vcf.gz

done

