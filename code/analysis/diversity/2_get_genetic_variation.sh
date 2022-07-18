#!/usr/bin/bash

## Measure of genetic diversity (pi and theta) and folded SFS

mkdir -p ../../../data/diversity/genetic_variation/par
mkdir -p ../../../data/diversity/genetic_variation/nonpar
mkdir -p ../../../data/diversity/genetic_variation/z
mkdir -p ../../../data/diversity/genetic_variation/chr4
mkdir -p ../../../data/diversity/genetic_variation/chr5

export outdir=../../../data/diversity/genetic_variation

# - PAR

for w_size in 200000 1000000
do
echo $w_size
python SFS_measures.py ../../../data/diversity/allele_count/par/PAR.intergenic.frq.count 20 ../../../data/genomic_features/z_scaf.${w_size}.intergene.overlap.density.sorted.txt \
all > ${outdir}/par/PAR.${w_size}.sfs.txt
done

# - nonPAR

for w_size in 200000 1000000
do
echo $w_size
python SFS_measures.py ../../../data/diversity/allele_count/nonpar/nonPAR.intergenic.frq.count 15 ../../../data/genomic_features/z_scaf.${w_size}.intergene.overlap.density.sorted.txt \
all > ${outdir}/nonpar/nonPAR.${w_size}.sfs.txt
done

# - Diversity measures for whole Z and conversion to Z coordinates

for w_size in 200000 1000000
do
echo $w_size
awk 'NR>1' ${outdir}/nonpar/nonPAR.${w_size}.sfs.txt | cat ${outdir}/par/PAR.${w_size}.sfs.txt - > ${outdir}/z/Z.${w_size}.sfs.txt 
python ../../processing/scaffold_to_chr.py ${outdir}/z/Z.${w_size}.sfs.txt Z > ${outdir}/z/Z.${w_size}.chr.coordinates.sfs.txt
python processing_coordinates.py ${outdir}/z/Z.${w_size}.sfs.txt > temp && mv -f temp ${outdir}/z/Z.${w_size}.chr.coordinates.sfs.txt
done


# - chr4

cat ../../../data/diversity/allele_count/chr4/*intergenic.frq.count | grep -v "CHROM" | awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $0}' \
> ../../../data/diversity/allele_count/chr4/chr4.intergenic.frq.count
for w_size in 200000
do
echo $w_size
python SFS_measures.py ../../../data/diversity/allele_count/chr4/chr4.intergenic.frq.count 20 ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.density.sorted.txt all > ${outdir}/chr4/chr4.${w_size}.sfs.txt
done


# - chr5

cat ../../../data/diversity/allele_count/chr5/*intergenic.frq.count | grep -v "CHROM" | awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $0}' \
> ../../../data/diversity/allele_count/chr5/chr5.intergenic.frq.count
for w_size in 200000 
do
echo $w_size
python SFS_measures.py ../../../data/diversity/allele_count/chr5/chr5.intergenic.frq.count 20 ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.density.sorted.txt all > ${outdir}/chr5/chr5.${w_size}.sfs.txt
done

## Get the distribution of SFS for PAR, nonPAR, chr4 and chr5

mkdir -p ../../../data/diversity/sfs/par
mkdir -p ../../../data/diversity/sfs/nonpar
mkdir -p ../../../data/diversity/sfs/chr4
mkdir -p ../../../data/diversity/sfs/chr5

export outdir_sfs=../../../data/diversity/sfs

echo "PAR"
python get_sfs.py ../../../data/diversity/allele_count/par/PAR.intergenic.frq.count black PAR > ${outdir_sfs}/par/PAR.intergenic.foldedSFS.txt
echo "nonPAR"
python get_sfs.py ../../../data/diversity/allele_count/nonpar/nonPAR.intergenic.frq.count black nonPAR > ${outdir_sfs}/nonpar/nonPAR.intergenic.foldedSFS.txt
echo "Autosome"
awk 'NR > 1' ../../../data/diversity/allele_count/chr5/chr5.intergenic.frq.count | cat ../../../data/diversity/allele_count/chr4/chr4.intergenic.frq.count - > ../../../data/diversity/allele_count/chr4_5_intergenic.frq.count
python get_sfs.py ../../../data/diversity/allele_count/chr4/chr4.intergenic.frq.count black A > ${outdir_sfs}/chr4/chr4.foldedSFS.txt
python get_sfs.py ../../../data/diversity/allele_count/chr5/chr5.intergenic.frq.count black A > ${outdir_sfs}/chr5/chr5.foldedSFS.txt
python get_sfs.py ../../../data/diversity/allele_count/chr4_5_intergenic.frq.count black A > ${outdir_sfs}/chr4_5.foldedSFS.txt
