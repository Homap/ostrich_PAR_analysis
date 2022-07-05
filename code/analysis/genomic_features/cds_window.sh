export fasta=../../../data/genome/black.repeat.depth.masked.fa

# PAR - nonPAR - Z
for seg in PAR nonPAR Z 
do
for w_size in 100000 200000 1000000
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/Z.CDS.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/${seg}.scaf.${w_size}.bed -b - -wao > ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.txt > ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.sum.txt ${fasta} ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.density.txt > ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.density.txt ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.sum.txt ../../../data/genomic_features/${seg}.scaf.${w_size}.CDS.overlap.txt
done
done

# Chr 4
for w_size in 100000 200000 1000000
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/chr4.CDS.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/chr4.${w_size}.bed -b - -wao > ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.txt > ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.sum.txt ${fasta} ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.density.txt > ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.density.txt ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.sum.txt ../../../data/genomic_features/chr4.scaf.${w_size}.CDS.overlap.txt
done

# Chr 5
for w_size in 100000 200000 1000000
do
echo $w_size
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../../../data/genomic_features/chr5.CDS.txt | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../../../data/sliding_window/chr5.${w_size}.bed -b - -wao > ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.txt
python get_sum_overlapping.py ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.txt > ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.sum.txt

python get_density_feature.py ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.sum.txt ${fasta} ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.density.txt > ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.density.sorted.txt

rm -f ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.density.txt ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.sum.txt ../../../data/genomic_features/chr5.scaf.${w_size}.CDS.overlap.txt
done
