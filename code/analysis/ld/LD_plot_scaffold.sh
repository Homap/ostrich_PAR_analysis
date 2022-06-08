# Concatenate all scaffolds for each subspecies LD

for scaffold in 26 54 35 36 62 63 67 69-1 83 88 92 93
do
cat black.superscaffold${scaffold}.pairwise.LD.05-500.200kbBin.50.kbStep.out >> LD_chromosome_plot/black.all.scaffolds.LD.05-500.200kbBin.50.kbStep.out
done 
#*************************************************************************************************************************************************************
# Correct the number format from scientific to decimal in LD output files
for subspecies in black blue red
do
echo $subspecies
python Scientific_to_Decimal.py ../data/LD/scaf.split/LD_chromosome_plot/${subspecies}.all.scaffolds.LD.05-500.200kbBin.50.kbStep.out > ../data/LD/scaf.split/LD_chromosome_plot/${subspecies}.all.scaffolds.LD.05-500.200kbBin.50.kbStep.Int.out
done
# Convert scaffold coordinate to chromosome
for subspecies in black blue red
do
python scaffold_to_chr_vcf_LD.py ../data/LD/scaf.split/LD_chromosome_plot/${subspecies}.all.scaffolds.LD.05-500.200kbBin.50.kbStep.Int.out > ../data/LD/scaf.split/LD_chromosome_plot/${subspecies}.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.out
done
# Correct column 3
awk '{if(NR>1) print $1"\t"$2"\t"$2+200000"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ../data/LD/scaf.split/LD_chromosome_plot/black.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.out > ../data/LD/scaf.split/LD_chromosome_plot/black.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out
awk '{if(NR>1) print $1"\t"$2"\t"$2+200000"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ../data/LD/scaf.split/LD_chromosome_plot/blue.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.out > ../data/LD/scaf.split/LD_chromosome_plot/blue.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out
awk '{if(NR>1) print $1"\t"$2"\t"$2+200000"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ../data/LD/scaf.split/LD_chromosome_plot/red.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.out > ../data/LD/scaf.split/LD_chromosome_plot/red.Z.coordinate.LD.05-500.200kbBin.50.kbStep.Int.v2.out