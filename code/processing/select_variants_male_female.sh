#!/bin/bash -l

#SBATCH -A snic2020-5-639
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J MFselect-filter
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
# Select only SNPs

module load bioinfo-tools
module load GATK/4.1.4.1

export reference=/proj/uppstore2017180/private/homap/ostrich_z/data/reference/Struthio_camelus.20130116.OM.fa
export vcf=/proj/uppstore2017180/private/homap/ostrich_z/data/vcf/all.bqsr.gatk.annotated.combined.snp.recal.vcf.gz
export outdir=/proj/uppstore2017180/private/homap/ostrich_z/data/vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_107 \
-sn P1878_108 \
-sn P1878_109 \
-sn P1878_110 \
-sn P1878_111 \
-O ${outdir}/black.male.vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_117 \
-sn P1878_118 \
-sn P1878_119 \
-sn P1878_120 \
-sn P1878_121 \
-O ${outdir}/blue.male.vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_127 \
-sn P1878_128 \
-sn P1878_129 \
-sn P1878_130 \
-sn P1878_131 \
-O ${outdir}/red.male.vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_112 \
-sn P1878_113 \
-sn P1878_114 \
-sn P1878_115 \
-sn P1878_116 \
-O ${outdir}/black.female.vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_123 \
-sn P1878_124 \
-sn P1878_125 \
-sn P1878_126 \
-O ${outdir}/blue.female.vcf

java -Xmx60g -jar $GATK_HOME/gatk-package-4.1.4.1-local.jar SelectVariants \
-R $reference \
-V $vcf \
-sn P1878_132 \
-sn P1878_134 \
-sn P1878_135 \
-sn P1878_136 \
-O ${outdir}/red.female.vcf
