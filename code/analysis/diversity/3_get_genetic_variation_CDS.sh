#!/usr/bin/bash

mkdir -p ../../../data/diversity/genetic_variation_cds/par
mkdir -p ../../../data/diversity/genetic_variation_cds/nonpar
mkdir -p ../../../data/diversity/genetic_variation_cds/z
mkdir -p ../../../data/diversity/genetic_variation_cds/chr4
mkdir -p ../../../data/diversity/genetic_variation_cds/chr5

export outdir=../../../data/diversity/genetic_variation_cds

# - Annotate PAR, nonPAR, chr4 and chr5 genomic sites into zerofold and fourfold sites
cut -f1 ../../../data/bed/par_scaf.bed | grep -f - ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff > ../../../data/genome/PAR.gff
cut -f1 ../../../data/bed/nonpar_scaf.bed | grep -f - ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff > ../../../data/genome/nonPAR.gff

cut -f1  ../../../data/lastz/gg_chr4_ostrich.bed | grep -f - ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff > ../../../data/genome/chr4.gff
cut -f1  ../../../data/lastz/gg_chr5_ostrich.bed | grep -f - ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff > ../../../data/genome/chr5.gff

python annotate/annotate_sites.py ../../../data/genome/Struthio_camelus.20130116.OM.fa ../../../data/genome/PAR.gff > ${outdir}/par/par.annotated.txt
python annotate/annotate_sites.py ../../../data/genome/Struthio_camelus.20130116.OM.fa ../../../data/genome/nonPAR.gff > ${outdir}/nonpar/nonpar.annotated.txt
python annotate/annotate_sites.py ../../../data/genome/Struthio_camelus.20130116.OM.fa ../../../data/genome/chr4.gff > ${outdir}/chr4/chr4.annotated.txt
python annotate/annotate_sites.py ../../../data/genome/Struthio_camelus.20130116.OM.fa ../../../data/genome/chr5.gff > ${outdir}/chr5/chr5.annotated.txt

# - Measure diversity in 4-fold and 0-fold sites

for seg in PAR chr4 chr5
do
seg_chr=$(echo $seg | tr '[:upper:]' '[:lower:]')
echo $seg_chr
python count_total_num_syn_nonsyn.py ${outdir}/${seg_chr}/${seg_chr}.annotated.txt ../../../data/diversity/allele_count/${seg_chr}/${seg}.frq.count 20 > ${outdir}/${seg_chr}/${seg_chr}.pi.txt
python fourfold_zerofold_diversity.py ${outdir}/${seg_chr}/${seg_chr}.pi.txt > ${outdir}/${seg_chr}/${seg_chr}.pnps.txt
done

python count_total_num_syn_nonsyn.py ${outdir}/nonpar/nonpar.annotated.txt ../../../data/diversity/allele_count/nonpar/nonPAR.frq.count 15 > ${outdir}/nonpar/nonpar.pi.txt
python fourfold_zerofold_diversity.py ${outdir}/nonpar/nonpar.pi.txt > ${outdir}/nonpar/nonpar.pnps.txt


# To calculate SFS for 4fold and 0fold, all I need is to know which SNP is 0fold and 4fold and overlap it with the allele count data
for seg in PAR chr4 chr5 
do
seg_chr=$(echo $seg | tr '[:upper:]' '[:lower:]')
echo $seg_chr
python count_0fold_4fold_overlap.py ../../../data/diversity/genetic_variation_cds/${seg_chr}/${seg_chr}.annotated.txt \
../../../data/diversity/allele_count/${seg_chr}/${seg}.frq.count 20 ../../../data/diversity/genetic_variation_cds/${seg_chr}

python get_sfs.py ../../../data/diversity/genetic_variation_cds/${seg_chr}/fourfold_counts.txt black ${seg} > ../../../data/diversity/genetic_variation_cds/${seg_chr}/fourfold_${seg_chr}_SFS.txt
python get_sfs.py ../../../data/diversity/genetic_variation_cds/${seg_chr}/zerofold_counts.txt black ${seg} > ../../../data/diversity/genetic_variation_cds/${seg_chr}/zerofold_${seg_chr}_SFS.txt
done

python count_0fold_4fold_overlap.py ../../../data/diversity/genetic_variation_cds/nonpar/nonpar.annotated.txt \
../../../data/diversity/allele_count/nonpar/nonPAR.frq.count 15 ../../../data/diversity/genetic_variation_cds/nonpar

python get_sfs.py ../../../data/diversity/genetic_variation_cds/nonpar/fourfold_counts.txt black nonPAR > ../../../data/diversity/genetic_variation_cds/nonpar/fourfold_nonpar_SFS.txt
python get_sfs.py ../../../data/diversity/genetic_variation_cds/nonpar/zerofold_counts.txt black nonPAR > ../../../data/diversity/genetic_variation_cds/nonpar/zerofold_nonpar_SFS.txt