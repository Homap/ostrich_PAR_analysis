module load bioinfo-tools BEDTools/2.29.2

export fasta=../../../data/genome/black.repeat.depth.masked.fa

#  Z
for seg in par nonpar z
do
for w_size in 200000 1000000
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/Z.intergene.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/${seg}_scaf_${w_size}_${w_size}.bed -b - -wao > ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.txt > ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.sum.txt ${fasta} ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.density.txt > ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.density.txt ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.sum.txt ../../../data/genomic_features/${seg}_scaf.${w_size}.intergene.overlap.txt
done
done

# chr4

for w_size in 200000 
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/chr4.intergene.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/chr4_${w_size}_${w_size}.bed -b - -wao > ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.txt > ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.sum.txt ${fasta} ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.density.txt > ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.density.txt ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.sum.txt ../../../data/genomic_features/chr4.scaf.${w_size}.intergene.overlap.txt
done

# Chr 5
for w_size in 200000
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/chr5.intergene.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/chr5_${w_size}_${w_size}.bed -b - -wao > ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.txt > ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.sum.txt ${fasta} ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.density.txt > ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.density.txt ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.sum.txt ../../../data/genomic_features/chr5.scaf.${w_size}.intergene.overlap.txt
done
