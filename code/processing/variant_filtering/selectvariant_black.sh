#!/bin/bash -l
#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J SelectVariants
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
 
module load bioinfo-tools GATK/4.1.4.1 

reference=/proj/snic2020-16-269/private/cornwallis.2020/data/external/ref/Struthio_camelus.20130116.OM.fa
vcf_in=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/black.all.vcf.gz
vcf_out=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/black.all.snp.vcf.gz

java -Xmx40g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf_in \
-sn P1878_107 \
-sn P1878_108 \
-sn P1878_109 \
-sn P1878_110 \
-sn P1878_111 \
-sn P1878_112 \
-sn P1878_113 \
-sn P1878_114 \
-sn P1878_115 \
-sn P1878_116 \
--select-type-to-include SNP \
-O $vcf_out 


