export fasta=../../data/genome/reference/black.repeat.depth.masked.fa


# Extract all regions with Z scaffolds, chr4 and chr5
cut -f1 ../../../data/bed/z_scaf.bed | awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' | grep -f - ../../../data/genomic_features/ostrich.all.CDS.coord | \
awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3"\t"a[2]"_"a[3]"\t"$4}' | sort -k1,1 -k2,2n > ../../../data/genomic_features/Z_genes_CDS.coord


#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
#["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
for window in 100000 200000 1000000
do
echo $window
bedtools intersect -a ../../../data/sliding_window/Z.${window}.bed -b ../../../data/genomic_features/Z_genes_CDS.chr.coord -wao > ../../../data/genomic_features/Z.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/Z.${window}.CDS.overlap.txt > ../data/bed/par_scaf.${window}.CDS.overlap.sum.txt


python get_density_feature.py ../data/bed/par_scaf.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.par_scaf.${window}.CDS.overlap.density.sorted.txt

done

# Convert to Z coordinates
python ../../processing/scaffold_to_chr.py ../../../data/genomic_features/Z_genes_CDS.coord Z > ../../../data/genomic_features/Z_genes_CDS.chr.coord

#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 1Mb
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../data/genomic_features/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/bed/autosome.${window}.bed -b - -wao > ../data/bed/autosome.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../data/bed/autosome.${window}.CDS.overlap.txt > ../data/bed/autosome.${window}.CDS.overlap.sum.txt

python get_density_feature.py ../data/bed/autosome.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.autosome.${window}.CDS.overlap.density.sorted.txt
done




#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 1Mb
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/bed/nonpar_scaf.${window}.bed -b - -wao > ../data/bed/nonpar_scaf.${window}.CDS.overlap.txt
python get_sum_overlapping.py ../data/bed/nonpar_scaf.${window}.CDS.overlap.txt > ../data/bed/nonpar_scaf.${window}.CDS.overlap.sum.txt

python get_density_feature.py ../data/bed/nonpar_scaf.${window}.CDS.overlap.sum.txt ${fasta} ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/bed/${subspecies}.nonpar_scaf.${window}.CDS.overlap.density.sorted.txt

done