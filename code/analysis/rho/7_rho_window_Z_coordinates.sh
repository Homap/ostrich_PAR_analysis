#!/usr/bin/bash

#-------------------------------------------------------------------------------------------------
# Window based Z 200 Kb in addition to 1Mb window for par
#-------------------------------------------------------------------------------------------------
export rho_z=../../../data/rho/ldhat_rho/z

# Concatenate all per scaffold rhos for PAR and nonPAR
for scaffold in 26 54 35
do
    awk 'NR>1' ${rho_z}/par/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
done

for scaffold in 62 67 69-1 93 63 88 83 92
do
    awk 'NR>1' ${rho_z}/nonpar/superscaffold${scaffold}.map.length.txt | cat - >> ${rho_z}/chrZ.map.length.txt
done

awk 'NR>1' ${rho_z}/par/superscaffold36.map.length.txt | cat ${rho_z}/nonpar/superscaffold36.map.length.txt - | awk 'NR>1' >> ${rho_z}/chrZ.map.length.txt
# Add the header to the concatenated rho file
awk 'BEGIN{print "Scaffold""\t""Locus_start""\t""Locus_end""\t""Mean_rho""\t""Median""\t""L95""\t""U95"}{print $0}' ${rho_z}/chrZ.map.length.txt > ${rho_z}/temp && mv -f ${rho_z}/temp ${rho_z}/chrZ.map.length.txt


#-------------------------------------------------------------------------------------------------
# 200 Kb with 50 Kb overlap for the comparison of rho with LD
#-------------------------------------------------------------------------------------------------
module load bioinfo-tools BEDTools/2.29.2

# Find overlap between windows and rho estimates
awk 'NR>1' ${rho_z}/chrZ.map.length.txt | bedtools intersect -a ../../../data/sliding_window/z_scaf_200000_50000.bed -b - -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/chrZ.map.length.Z.coord.200kb50Kb.txt

# Summarise windows and their overlapping regions
python get_recombination_per_window.py ${rho_z}/chrZ.map.length.Z.coord.200kb50Kb.txt \
../../../data/sliding_window/z_scaf_200000_50000.bed > ${rho_z}/200Kb50kb_rho_Z.txt

# Conversion of scaffold to chromosome coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/200Kb50kb_rho_Z.txt Z > ${rho_z}/200Kb50Kb_rho_Z.chr.coord.txt

# Remove intermediate files
rm -f ${rho_z}/chrZ.map.length.Z.coord.200kb50Kb.txt ${rho_z}/200Kb50kb_rho_Z.txt

#-------------------------------------------------------------------------------------------------
# 200 Kb without overlap for the comparison of rho with genetic diversity and genomic features
#-------------------------------------------------------------------------------------------------
# Find overlap between windows and rho estimates
awk 'NR>1' ${rho_z}/chrZ.map.length.txt | bedtools intersect -a ../../../data/sliding_window/z_scaf_200000_200000.bed -b - -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/chrZ.map.length.Z.coord.200kb200Kb.txt

# Summarise windows and their overlapping regions
python get_recombination_per_window.py ${rho_z}/chrZ.map.length.Z.coord.200kb200Kb.txt \
../../../data/sliding_window/z_scaf_200000_200000.bed > ${rho_z}/200Kb200kb_rho_Z.txt

# Conversion of scaffold to chromosome coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/200Kb200kb_rho_Z.txt Z > ${rho_z}/200Kb200Kb_rho_Z.chr.coord.txt

# Remove intermediate files
rm -f ${rho_z}/chrZ.map.length.Z.coord.200kb200Kb.txt ${rho_z}/200Kb200kb_rho_Z.txt

#-------------------------------------------------------------------------------------------------
# 1Mb PAR window without overlap for using rho with recombination frequency from genetic map to obtain Ne and use in simulation
#-------------------------------------------------------------------------------------------------
# Find the overlap between scaffold windows and the rho
awk 'NR>1' ${rho_z}/chrZ.map.length.txt | bedtools intersect -a ../../../data/sliding_window/par_scaf_1000000_1000000.bed -b - -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/1Mb1Mb_scaf.PAR.txt

# Summarise windows and their overlapping regions
python get_recombination_per_window.py ${rho_z}/1Mb1Mb_scaf.PAR.txt \
../../../data/sliding_window/par_scaf_1000000_1000000.bed > ${rho_z}/1Mb1Mb_rho_scaf.PAR.txt

# Conversion of scaffold to chromosome coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/1Mb1Mb_rho_scaf.PAR.txt PAR > ${rho_z}/1Mb1Mb_rho_chrom.PAR.txt

# Re-calculate windows based on Z
bedtools intersect -a ../../../data/sliding_window/z_par_1000000_1000000.bed  -b ${rho_z}/1Mb1Mb_rho_chrom.PAR.txt -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""rho_per_site""\t""Mean_rho""\t""overlap"}{print $0}' > ${rho_z}/step2.PAR.1000000.txt

python get_recombination_per_window_step2.py ${rho_z}/step2.PAR.1000000.txt ../../../data/sliding_window/z_par_1000000_1000000.bed > ${rho_z}/1Mb1Mb.rho.chrom.PAR.txt

# Remove intermediate files
rm -f ${rho_z}/1Mb1Mb_scaf.PAR.txt ${rho_z}/1Mb1Mb_rho_scaf.PAR.txt ${rho_z}/1Mb1Mb_rho_chrom.PAR.txt ${rho_z}/step2.PAR.1000000.txt 

#-------------------------------------------------------------------------------------------------
# 1Mb Z window without overlap for comparing rho with genetic map
#-------------------------------------------------------------------------------------------------
# Find the overlap between scaffold windows and the rho
awk 'NR>1' ${rho_z}/chrZ.map.length.txt | bedtools intersect -a ../../../data/sliding_window/z_scaf_1000000_1000000.bed -b - -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/1Mb1Mb_scaf.Z.txt

# Summarise windows and their overlapping regions
python get_recombination_per_window.py ${rho_z}/1Mb1Mb_scaf.Z.txt \
../../../data/sliding_window/z_scaf_1000000_1000000.bed > ${rho_z}/1Mb1Mb_rho_scaf.Z.txt

# Conversion of scaffold to chromosome coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/1Mb1Mb_rho_scaf.Z.txt PAR > ${rho_z}/1Mb1Mb_rho_chrom.Z.txt

# Re-calculate windows based on Z
bedtools intersect -a ../../../data/sliding_window/z_chrom_1000000_1000000.bed  -b ${rho_z}/1Mb1Mb_rho_chrom.Z.txt -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""rho_per_site""\t""Mean_rho""\t""overlap"}{print $0}' > ${rho_z}/step2.Z.1000000.txt

python get_recombination_per_window_step2.py ${rho_z}/step2.Z.1000000.txt ../../../data/sliding_window/z_chrom_1000000_1000000.bed > ${rho_z}/1Mb1Mb.rho.chrom.Z.txt

# Remove intermediate files
rm -f ${rho_z}/1Mb1Mb_scaf.Z.txt ${rho_z}/1Mb1Mb_rho_scaf.Z.txt ${rho_z}/1Mb1Mb_rho_chrom.Z.txt ${rho_z}/step2.Z.1000000.txt 