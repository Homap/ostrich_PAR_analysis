#!/bin/bash -l

mkdir -p ../../../data/ld/ld_window/autosome/chr4
mkdir -p ../../../data/ld/ld_window/autosome/chr5
mkdir -p ../../../data/ld/ld_window/z

# Chromosome 4 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr4_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/a_vcf/chr4/${scaffold}.chr4.vcf.gz \
../../../data/ld/ld_window/autosome/chr4/${scaffold}.pairwise
done

# Chromosome 5 -------------------------------------------------------------------------------------------
for scaffold in $(cat ../../../data/lastz/gg_chr5_ostrich.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/a_vcf/chr5/${scaffold}.chr5.vcf.gz \
../../../data/ld/ld_window/autosome/chr5/${scaffold}.pairwise
done

# Separating Z scaffolds ---------------------------------------------------------------------------------
# chromosome Z
for scaffold in $(cat ../../../data/bed/z_scaf.bed | grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
sbatch -J ${scaffold}.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/${scaffold}.chrz.vcf.gz \
../../../data/ld/ld_window/z/${scaffold}.pairwise
done

# PAR : This is the same as the scaffolds in Z except for superscaffold36
sbatch -J superscaffold36.par.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.par.vcf.gz \
../../../data/ld/ld_window/z/superscaffold36.par.pairwise
sbatch -J superscaffold36.nonpar.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/superscaffold36.nonpar.vcf.gz \
../../../data/ld/ld_window/z/superscaffold36.nonpar.pairwise

# (100 Kb and 500 Kb)
# PAR farthest away
sbatch -J supersaffold26.500Kb.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold26.500Kb.vcf.gz \
../../../data/ld/ld_window/z/supersaffold26.500Kb.pairwise

# PAR mid region
sbatch -J supersaffold54.500Kb.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold54.500Kb.vcf.gz \
../../../data/ld/ld_window/z/supersaffold54.500Kb.pairwise

# PAR boundary in PAR
sbatch -J supersaffold36.500Kb.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold36.500Kb.vcf.gz \
../../../data/ld/ld_window/z/supersaffold36.500Kb.pairwise

# PAR boundary spanning to nonPAR
sbatch -J supersaffold36.100Kb.LD pairwise_ld.sh ../../../data/vcf/ld_vcf/z_vcf/supersaffold36.100Kb.spanning.boundary.vcf.gz \
../../../data/ld/ld_window/z/supersaffold36.100Kb.spanning.boundary.pairwise



