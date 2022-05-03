#!/usr/bin/bash

#*********************************************************************************************#
# Calculating genetic diversity for males and females
#*********************************************************************************************#
module load bioinfo-tools BEDTools
# Get overlap between windows and element
#*********************************************************************************************#
# Repeats
echo "Repeats"
#*********************************************************************************************#
# merge overlapping repeat coordinates
repeat_addr='../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.out'
awk -v OFS="\t" '$1=$1' ${repeat_addr} | awk '{$6=$6-1; print $5"\t"$6"\t"$7}' | awk 'NR>2' | mergeBed -i - | grep -f ../data/bed/Z_scaffolds.txt > ../data/bed/SC.allRepeats.bed
for window in 200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/window/ostrich.Z.${window}.bed -b ../data/bed/SC.allRepeats.bed -wao > ../data/pi/ostrich.Z.${window}.AllRepeat.overlap.txt
python get_sum_overlapping.py ../data/pi/ostrich.Z.${window}.AllRepeat.overlap.txt > ../data/pi/ostrich.Z.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the outout
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file. 
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/repeat_coverage_masked_fasta/${subspecies}.repeat.cov.masked.fa
python get_density_feature.py ../data/pi/ostrich.Z.${window}.AllRepeat.overlap.sum.txt ${fasta} ../data/pi/${subspecies}.Z.${window}.AllRepeat.overlap.density.txt 
sort -k1,1 -k2n,2 ../data/pi/${subspecies}.Z.${window}.AllRepeat.overlap.density.txt | awk 'NR>1' > ../data/pi/${subspecies}.Z.${window}.AllRepeat.overlap.density.sorted.txt
done
done




# CDS coordinates

# Intron coordinates

# Intergenic coordinates

# Hard-masking background

# Get overlap of windows with genomic features

# Annotate SNPs for the genomic features

# Calculate pi per site for males

# Calculate pi pere site for females

# Resampling interegenic pi to create confidence intervals

# Prepare data for plotting in R

# Calculate FST between males and females
