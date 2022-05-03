#!/bin/bash -l
#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 14:00:00
#SBATCH -J window_overlap
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools vcftools BEDTools 

# For each window, get the total number of sites with:
# 1. All bases filtered for coverage and repeats
# 2. Intergenic bases filtered for coverage and repeats
# 3. Intronic bases filtered for coverage and repeats
#*********************************************************************************************#
# Get windows per scaffolds
#*********************************************************************************************#
echo "Create windows"
for window in 1Kb 10Kb 200Kb 500Kb 1000Kb 100Kb 
do
w_size=$(echo $window | sed 's/Kb/000/g')
#python sliding_window_ZA.py ../data/bed/autosomes.bed ${w_size} > ../data/bed/autosome.${window}.bed
python sliding_window_ZA.py ../data/bed/par_scaf.bed ${w_size} > ../data/bed/par_scaf.${window}.bed
python sliding_window_ZA.py ../data/bed/nonpar_scaf.bed ${w_size} > ../data/bed/nonpar_scaf.${window}.bed
done

for window in 0.5Kb 
do
w_size=500
python sliding_window_ZA.py ../data/bed/par_scaf.bed ${w_size} > ../data/bed/par_scaf.${window}.bed
done
#*********************************************************************************************#
#*********************************************************************************************#
# Autosomes get overlap between windows and element
#*********************************************************************************************#
# Repeats
echo "Repeats"
#*********************************************************************************************#
# merge overlapping repeat coordinates
for window in 200Kb 500Kb 1000Kb #100Kb
do
echo ${window}
bedtools intersect -a ../data/bed/autosome.${window}.bed -b ../data/bed/ostrich_repeats.bed -wao > ../data/bed/autosome.${window}.AllRepeat.overlap.txt
python get_sum_overlapping.py ../data/bed/autosome.${window}.AllRepeat.overlap.txt > ../data/bed/autosome.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the output
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file.
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/autosome.${window}.AllRepeat.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.autosome.${window}.AllRepeat.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.autosome.${window}.AllRepeat.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.autosome.${window}.AllRepeat.overlap.density.sorted.txt
# output header is: ["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
done
done
#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
for window in 500Kb 1000Kb #100Kb #200Kb 
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/autosome.${window}.bed -b - -wao > ../data/bed/autosome.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../data/bed/autosome.${window}.CDS.overlap.txt > ../data/bed/autosome.${window}.CDS.overlap.sum.txt
for subspecies in black #blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/autosome.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intron
echo "Intron"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb #200Kb 500Kb 1000Kb
do
echo ${window}
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/autosome.${window}.bed -b - -wao > ../data/bed/autosome.${window}.intronic.overlap
python get_sum_overlapping.py ../data/bed/autosome.${window}.intronic.overlap > ../data/bed/autosome.${window}.intronic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/autosome.${window}.intronic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.autosome.${window}.intronic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.autosome.${window}.intronic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.autosome.${window}.intronic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intergenic
echo "Intergenic"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb #200Kb 500Kb 1000Kb
do
bedtools intersect -a ../data/bed/autosome.${window}.bed -b ../data/bed/ostrich.intergene.coord -wao > ../data/bed/autosome.${window}.intergenic.overlap
python get_sum_overlapping.py ../data/bed/autosome.${window}.intergenic.overlap > ../data/bed/autosome.${window}.intergenic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/autosome.${window}.intergenic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.autosome.${window}.intergenic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.autosome.${window}.intergenic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.autosome.${window}.intergenic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
#PAR
#*********************************************************************************************#
# Repeats
echo "Repeats"
#*********************************************************************************************#
# merge overlapping repeat coordinates
for window in 100Kb #200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/bed/par_scaf.${window}.bed -b ../data/bed/ostrich_repeats.bed -wao > ../data/bed/par_scaf.${window}.AllRepeat.overlap.txt
python get_sum_overlapping.py ../data/bed/par_scaf.${window}.AllRepeat.overlap.txt > ../data/bed/par_scaf.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the output
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file.
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/par_scaf.${window}.AllRepeat.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.par_scaf.${window}.AllRepeat.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.par_scaf.${window}.AllRepeat.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.par_scaf.${window}.AllRepeat.overlap.density.sorted.txt
# output header is: ["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
done
done
#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
#["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
for window in  500Kb 1000Kb #100Kb #200Kb
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/par_scaf.${window}.bed -b - -wao > ../data/bed/par_scaf.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../data/bed/par_scaf.${window}.CDS.overlap.txt > ../data/bed/par_scaf.${window}.CDS.overlap.sum.txt
for subspecies in black #blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/par_scaf.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intron
echo "Intron"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb # 200Kb 500Kb 1000Kb
do
echo ${window}
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/par_scaf.${window}.bed -b - -wao > ../data/bed/par_scaf.${window}.intronic.overlap
python get_sum_overlapping.py ../data/bed/par_scaf.${window}.intronic.overlap > ../data/bed/par_scaf.${window}.intronic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/par_scaf.${window}.intronic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.par_scaf.${window}.intronic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.par_scaf.${window}.intronic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.par_scaf.${window}.intronic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intergenic
echo "Intergenic"
#*************************************************************************************************************************************************************************************************************************
for window in 0.5Kb #1Kb 10Kb #500Kb 1000Kb #100Kb 
do
bedtools intersect -a ../data/bed/par_scaf.${window}.bed -b ../data/bed/ostrich.intergene.coord -wao > ../data/bed/par_scaf.${window}.intergenic.overlap
python get_sum_overlapping.py ../data/bed/par_scaf.${window}.intergenic.overlap > ../data/bed/par_scaf.${window}.intergenic.overlap.sum.txt
for subspecies in black #blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/par_scaf.${window}.intergenic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.par_scaf.${window}.intergenic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.par_scaf.${window}.intergenic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.par_scaf.${window}.intergenic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# nonPAR
# Repeats
echo "Repeats"
#*********************************************************************************************#
# merge overlapping repeat coordinates
for window in 100Kb #200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/bed/nonpar_scaf.${window}.bed -b ../data/bed/ostrich_repeats.bed -wao > ../data/bed/nonpar_scaf.${window}.AllRepeat.overlap.txt
python get_sum_overlapping.py ../data/bed/nonpar_scaf.${window}.AllRepeat.overlap.txt > ../data/bed/nonpar_scaf.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the output
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file.
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/nonpar_scaf.${window}.AllRepeat.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.nonpar_scaf.${window}.AllRepeat.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.nonpar_scaf.${window}.AllRepeat.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.nonpar_scaf.${window}.AllRepeat.overlap.density.sorted.txt
# output header is: ["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
done
done
#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 500Kb 1000Kb
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/nonpar_scaf.${window}.bed -b - -wao > ../data/bed/nonpar_scaf.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../data/bed/nonpar_scaf.${window}.CDS.overlap.txt > ../data/bed/nonpar_scaf.${window}.CDS.overlap.sum.txt
for subspecies in black #blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/nonpar_scaf.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intron
echo "Intron"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 500Kb 1000Kb
do
echo ${window}
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/nonpar_scaf.${window}.bed -b - -wao > ../data/bed/nonpar_scaf.${window}.intronic.overlap
python get_sum_overlapping.py ../data/bed/nonpar_scaf.${window}.intronic.overlap > ../data/bed/nonpar_scaf.${window}.intronic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/nonpar_scaf.${window}.intronic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.nonpar_scaf.${window}.intronic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.nonpar_scaf.${window}.intronic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.nonpar_scaf.${window}.intronic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intergenic
echo "Intergenic"
#*************************************************************************************************************************************************************************************************************************
for window in 1Kb 10Kb #500Kb 1000Kb #100Kb 
do
bedtools intersect -a ../data/bed/nonpar_scaf.${window}.bed -b ../data/bed/ostrich.intergene.coord -wao > ../data/bed/nonpar_scaf.${window}.intergenic.overlap
python get_sum_overlapping.py ../data/bed/nonpar_scaf.${window}.intergenic.overlap > ../data/bed/nonpar_scaf.${window}.intergenic.overlap.sum.txt
for subspecies in black #blue red
do
echo ${subspecies}
fasta=../data/reference/${subspecies}.repeat.depth.masked.fa
python get_density_feature.py ../data/bed/nonpar_scaf.${window}.intergenic.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.nonpar_scaf.${window}.intergenic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.nonpar_scaf.${window}.intergenic.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.nonpar_scaf.${window}.intergenic.overlap.density.sorted.txt
done
done

