#!/usr/bin/bash

module load bioinfo-tools BEDTools
for species in black blue red
do 

# Autosome
echo $species
echo "Autosome"
awk '{if(NR>1) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.A.repeat.frq.count | \
bedtools intersect -a - -b ../data/bed/ostrich.intergene.coord -wao | \
awk '{if($11=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""INTERGENIC"}' > \
../data/allele_count/${species}.A.repeat.intergenic.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.A.repeat.intergenic.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""INTRON"}' > \
../data/allele_count/${species}.A.repeat.intergenic.intron.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.A.repeat.intergenic.intron.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""CDS"}' \
> ../data/allele_count/${species}.A.repeat.intergenic.intron.cds.frq.count

grep 'INTERGENIC' ../data/allele_count/${species}.A.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.A.intergenic
grep 'INTRON' ../data/allele_count/${species}.A.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.A.intronic
grep 'CDS' ../data/allele_count/${species}.A.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.A.cds

# PAR
echo "PAR"
awk '{if(NR>1) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.PAR.repeat.frq.count | \
bedtools intersect -a - -b ../data/bed/ostrich.intergene.coord -wao | \
awk '{if($11=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""INTERGENIC"}' > \
../data/allele_count/${species}.PAR.repeat.intergenic.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.PAR.repeat.intergenic.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""INTRON"}' > \
../data/allele_count/${species}.PAR.repeat.intergenic.intron.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.PAR.repeat.intergenic.intron.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""CDS"}' \
> ../data/allele_count/${species}.PAR.repeat.intergenic.intron.cds.frq.count

grep 'INTERGENIC' ../data/allele_count/${species}.PAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.PAR.intergenic
grep 'INTRON' ../data/allele_count/${species}.PAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.PAR.intronic
grep 'CDS' ../data/allele_count/${species}.PAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.PAR.cds

# nonPAR
echo "nonPAR"
awk '{if(NR>1) print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count | \
bedtools intersect -a - -b ../data/bed/ostrich.intergene.coord -wao | \
awk '{if($11=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""INTERGENIC"}' > \
../data/allele_count/${species}.nonPAR.repeat.intergenic.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.nonPAR.repeat.intergenic.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t""INTRON"}' > \
../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.frq.count

awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/bed/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | \
bedtools intersect -a ../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.frq.count -b - -wao | \
awk '{if($13=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""CDS"}' \
> ../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.cds.frq.count

grep 'INTERGENIC' ../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.intergenic
grep 'INTRON' ../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.intronic
grep 'CDS' ../data/allele_count/${species}.nonPAR.repeat.intergenic.intron.cds.frq.count | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.cds
done

echo "Calculating SFS measures"
for species in black 
do
echo $species
echo "Autosome"
python SFS_measures.py ../data/allele_count/${species}.A.intergenic 20 ../data/bed/${species}.autosome.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.A.intronic 20 ../data/bed/${species}.autosome.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.A.cds 20 ../data/bed/${species}.autosome.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.cds.sfs.txt
echo "PAR"
python SFS_measures.py ../data/allele_count/${species}.PAR.intergenic 20 ../data/bed/${species}.par_scaf.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.PAR.intronic 20 ../data/bed/${species}.par_scaf.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.PAR.cds 20 ../data/bed/${species}.par_scaf.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.cds.sfs.txt
echo "nonPAR"
python SFS_measures.py ../data/allele_count/${species}.nonPAR.intergenic 15 ../data/bed/${species}.nonpar_scaf.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.nonPAR.intronic 15 ../data/bed/${species}.nonpar_scaf.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.nonPAR.cds 15 ../data/bed/${species}.nonpar_scaf.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.cds.sfs.txt
done

for species in blue red
do
echo $species
echo "Autosome"
python SFS_measures.py ../data/allele_count/${species}.A.intergenic 18 ../data/bed/${species}.autosome.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.A.intronic 18 ../data/bed/${species}.autosome.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.A.cds 18 ../data/bed/${species}.autosome.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.A.cds.sfs.txt
echo "PAR"
python SFS_measures.py ../data/allele_count/${species}.PAR.intergenic 18 ../data/bed/${species}.par_scaf.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.PAR.intronic 18 ../data/bed/${species}.par_scaf.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.PAR.cds 18 ../data/bed/${species}.par_scaf.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.PAR.cds.sfs.txt
echo "nonPAR"
python SFS_measures.py ../data/allele_count/${species}.nonPAR.intergenic 14 ../data/bed/${species}.nonpar_scaf.100Kb.intergenic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.intergenic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.nonPAR.intronic 14 ../data/bed/${species}.nonpar_scaf.100Kb.intronic.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.intronic.sfs.txt
python SFS_measures.py ../data/allele_count/${species}.nonPAR.cds 14 ../data/bed/${species}.nonpar_scaf.100Kb.CDS.overlap.density.sorted.txt \
> ../data/allele_count/${species}.nonPAR.cds.sfs.txt
done