#!/bin/bash -l
#SBATCH -A snic2020-5-639
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J FilterVariants
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
 
module load bioinfo-tools GATK/4.1.4.1 

reference=/proj/snic2020-16-269/private/cornwallis.2020/data/external/ref/Struthio_camelus.20130116.OM.fa
vcf_in=$1
vcf_out=$2
 
java -Xmx50g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar VariantFiltration \
-R $reference \
-V $vcf_in \
-O $vcf_out \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
--filter-name "snp_filters"

 
