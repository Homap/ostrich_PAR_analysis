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

# Conversion of scaffold to chromosome coordinates
# For plotting we need the whole Z coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.map.length.txt Z > ${rho_z}/chrZ.map.length.Z.coord.txt

# For simulation we need the PAR coordinates
python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.map.length.txt PAR > ${rho_z}/chrZ.map.length.PAR.coord.txt

# For statistical analyses, we need both PAR and nonPAR
python ../../processing/scaffold_to_chr.py ${rho_z}/chrZ.map.length.txt nonPAR > ${rho_z}/chrZ.map.length.nonPAR.coord.txt

# Produce sliding windows
export sw=../../../data/sliding_window
# For simulation, we use window size of 1Mb since smaller than that, the Ne calculation is too noisy
python ../../processing/sliding_window.py ../../../data/bed/z_par.bed 1000000 1000000 > ${sw}/Z.coord.1Mb.windows.PAR.txt

# For statistic, comparison between rho, genetic map and LD and genetic diversity, we use 200 Kb windows
python ../../processing/sliding_window.py ../../../data/bed/z_chrom.bed 200000 200000 > ${sw}/Z.coord.200Kb.windows.Z.txt
python ../../processing/sliding_window.py ../../../data/bed/z_par.bed 200000 200000  > ${sw}/Z.coord.200Kb.windows.PAR.txt
python ../../processing/sliding_window.py ../../../data/bed/z_nonpar.bed 200000 200000  > ${sw}/Z.coord.200Kb.windows.nonPAR.txt

# To find overlap between rho and a given window:
module load bioinfo-tools BEDTools/2.29.2
bedtools intersect -a ${sw}/Z.coord.1Mb.windows.PAR.txt -b ${rho_z}/chrZ.map.length.PAR.coord.txt -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/chrZ.map.length.PAR.coord.1Mbwindow.txt

for chr_seg in Z PAR nonPAR 
do
bedtools intersect -a ${sw}/Z.coord.200Kb.windows.${chr_seg}.txt -b ${rho_z}/chrZ.map.length.${chr_seg}.coord.txt -wao | \
awk 'BEGIN{print "CHROM""\t""win_start""\t""win_end""\t""chr""\t""start""\t""end""\t""scaffold""\t""start""\t""end""\t""Mean_rho""\t""Median""\t""L95""\t""U95""\t""overlap"}{print $0}' > ${rho_z}/chrZ.map.length.${chr_seg}.coord.200Kb.window.txt
done

# Summarise windows and their overlapping regions
python get_recombination_per_window.py ${rho_z}/chrZ.map.length.PAR.coord.1Mbwindow.txt ${sw}/Z.coord.1Mb.windows.PAR.txt > ${rho_z}/1Mb_rho_PAR.txt

python get_recombination_per_window.py ${rho_z}/chrZ.map.length.Z.coord.200Kb.window.txt ${sw}/Z.coord.200Kb.windows.Z.txt > ${rho_z}/200Kb_rho_Z.txt

python get_recombination_per_window.py ${rho_z}/chrZ.map.length.PAR.coord.200Kb.window.txt ${sw}/Z.coord.200Kb.windows.PAR.txt > ${rho_z}/200Kb_rho_PAR.txt

python get_recombination_per_window.py ${rho_z}/chrZ.map.length.nonPAR.coord.200Kb.window.txt ${sw}/Z.coord.200Kb.windows.nonPAR.txt > ${rho_z}/200Kb_rho_nonPAR.txt

# remove intermediate files
rm -f ${rho_z}/chrZ.map.length.Z.coord.200Kb.window.txt \
{rho_z}/chrZ.map.length.PAR.coord.200Kb.window.txt \
{rho_z}/chrZ.map.length.nonPAR.coord.200Kb.window.txt \
${rho_z}/chrZ.map.length.Z.coord.txt \
${rho_z}/chrZ.map.length.PAR.coord.txt \
${rho_z}/chrZ.map.length.nonPAR.coord.txt \
${rho_z}/chrZ.map.length.txt