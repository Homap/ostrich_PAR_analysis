#!/bin/bash -l

mkdir -p ../../../data/ld/ld_decay/autosome/chr4
mkdir -p ../../../data/ld/ld_decay/autosome/chr5
mkdir -p ../../../data/ld/ld_decay/z

# Chromosome 4 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr4_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD ld_decay.sh ../../../data/vcf/ld_vcf/a_vcf/chr4/${scaffold}.chr4.vcf.gz \
../../../data/ld/ld_decay/autosome/chr4/${scaffold}.LDdecay \
../../../data/ld/ld_decay/autosome/chr4/${scaffold}.LDdecay.stat.gz \
../../../data/ld/ld_decay/autosome/chr4/${scaffold}.LDdecay
done

# Chromosome 5 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr5_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD ld_decay.sh ../../../data/vcf/ld_vcf/a_vcf/chr5/${scaffold}.chr5.vcf.gz \
../../../data/ld/ld_decay/autosome/chr5/${scaffold}.LDdecay \
../../../data/ld/ld_decay/autosome/chr5/${scaffold}.LDdecay.stat.gz \
../../../data/ld/ld_decay/autosome/chr5/${scaffold}.LDdecay
done

# Separating Z scaffolds ---------------------------------------------------------------------------------
# chromosome Z
for scaffold in $(cat ../../../data/bed/z_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/${scaffold}.chrz.vcf.gz \
../../../data/ld/ld_decay/z/${scaffold}.LDdecay \
../../../data/ld/ld_decay/z/${scaffold}.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/${scaffold}.LDdecay
done

# PAR : This is the same as the scaffolds in Z except for superscaffold36
sbatch -J superscaffold36.par.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.par.vcf.gz \
../../../data/ld/ld_decay/z/superscaffold36.par.LDdecay \
../../../data/ld/ld_decay/z/superscaffold36.par.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/superscaffold36.par.LDdecay
sbatch -J superscaffold36.nonpar.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.nonpar.vcf.gz \
../../../data/ld/ld_decay/z/superscaffold36.nonpar.LDdecay \
../../../data/ld/ld_decay/z/superscaffold36.nonpar.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/superscaffold36.nonpar.LDdecay

# (100 Kb and 500 Kb)
# PAR farthest away
sbatch -J supersaffold26.500Kb.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold26.500Kb.vcf.gz \
../../../data/ld/ld_decay/z/supersaffold26.500Kb.LDdecay \
../../../data/ld/ld_decay/z/supersaffold26.500Kb.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/supersaffold26.500Kb.LDdecay 

# PAR mid region
sbatch -J supersaffold54.500Kb.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold54.500Kb.vcf.gz \
../../../data/ld/ld_decay/z/supersaffold54.500Kb.LDdecay \
../../../data/ld/ld_decay/z/supersaffold54.500Kb.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/supersaffold54.500Kb.LDdecay

# PAR boundary in PAR
sbatch -J supersaffold36.500Kb.LD ld_decay.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold36.500Kb.vcf.gz \
../../../data/ld/ld_decay/z/supersaffold36.500Kb.LDdecay \
../../../data/ld/ld_decay/z/supersaffold36.500Kb.LDdecay.stat.gz \
../../../data/ld/ld_decay/z/supersaffold36.500Kb.LDdecay


