# Linkage Disequilibrium (LD)

We would like to compute LD on the autosome and Z chromosmoe. In ostrich genome assembly, only chromosome Z has been
developed into a chromosome-level assembly by use of linkage map. We lack a chromosome level assembly for the autosomes. 
In ostrich and most other birds with studied karyotype, chromosome Z has similar physical and genetic length to chromosomes
4 and 5. Using LASTZ, we obtain the match between ostrich scaffolds and chicken chromosomes 4 and 5. Chimeric scaffolds, scaffolds
containing sequences from more than one chromosome, might exist in assemblies. We therefore filter out chimeric scaffolds and keep scaffolds
with only one hit to either chromosome 4 or 5 of chicken.

The sequence in the PAR is complex 


In the following, I use PopLDdecay, following the analysis done by
Takeshi Kawakami to calculate LD on ostrich Z chromosome
Download PopLDdecay into /proj/snic2020-16-269/private/homap/ostrich_z/bin

```
git clone https://github.com/BGI-shenzhen/PopLDdecay.git
cd PopLDdecay; chmod 755 configure; ./configure;
make;
mv PopLDdecay  bin/;
```

I use the r-squared measure of LD.

export z_vcf=../../data/vcf/z_vcf/z_vcf.gz

## Concatenate PAR and non-PAR positions
```
zgrep -v "#" ../data/vcf/LD_vcf/${species}.A.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.A.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.PAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.PAR.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.nonPAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.nonPAR.pos

cat ../data/vcf/LD_vcf/${species}.PAR.pos ../data/vcf/LD_vcf/${species}.nonPAR.pos | \
sort -k1,1 -k2n > ../data/vcf/LD_vcf/${species}.Z.pos

awk '{print $1"\t"$2-1"\t"$2}' ../data/vcf/LD_vcf/${species}.Z.pos > ../data/vcf/LD_vcf/${species}.Z.pos.bed 
python scaffold_to_chr_vcf.py ../data/vcf/LD_vcf/${species}.Z.pos.bed > ../data/vcf/LD_vcf/${species}.Z.pos.coordinates.txt
```

## Analysis of LD
```
sbatch ld_run.sh
````

## Calculate mean LD in 200 kb windows by sliding window analysis
```
for scaffold in $(cat ../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed \
| grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
zcat ../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.gz | \
grep -v "^#" | awk '($9>=500 && $9<=50000)' > \
../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.05-500
Rscript LD_slidingWin_v1.R ../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.05-500 \
../data/bed/z_scaf.bed 200000 50000
done
```

## LD at the PAR boundary

Next, we are interested to get a more clear picture of LD across the PAR boundary in each sex. 
```
module load bioinfo-tools vcftools
mkdir -p ../data/LD/sex_specific
```
## Select SNPs that 100 Kb or 50 Kb left and right of the PAR boundary
```
cat ../data/vcf/LD_vcf/black.Z.pos | grep "^superscaffold36" | \
awk '($2-3524263>-100000 && $2-3524263<100000)' > ../data/LD/sex_specific/commonSNP.B.boundary200kb

cat ../data/vcf/LD_vcf/black.Z.pos | grep "^superscaffold36" | \
awk '($2-3524263>-50000 && $2-3524263<50000)' > ../data/LD/sex_specific/commonSNP.B.boundary100kb
```

## Convert vcf to plink near PAR boundary (200 kb centered at the boundary)
```
chr=superscaffold36
vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.superscaffold36.vcf.gz --positions ../data/LD/sex_specific/commonSNP.B.boundary100kb \
--plink --out ../data/LD/sex_specific/black.superscaffold36.bothsexes.100kb.boundary

plink --file ../data/LD/sex_specific/black.superscaffold36.bothsexes.100kb.boundary --r2 square --out ../data/LD/sex_specific/black.superscaffold36.bothsexes.100Kb.boundary

awk '{print $1"\t"$2-1"\t"$2}' ../data/LD/sex_specific/commonSNP.B.boundary100kb | cat ../data/bed/par_scaf.bed - > test 

python scaffold_to_chr_vcf.py test > ../data/bed/boundary.z.coordinates.txt
```

## Join ped files per sex
```
cat ../data/LD/sex_specific/*male.ped > \
../data/LD/sex_specific/BB.superscaffold36.boundary.male.ped.join
cat ../data/LD/sex_specific/*female.ped > \
../data/LD/sex_specific/BB.superscaffold36.boundary.female.ped.join
```

Modify plink infile for Haploview. Note that Haploview takes plink format. However
.ped file needs a header. Instead of using plink format, I use Linkage format that
needs .ped file as is and .info file that requires two columns indicating SNP name
and position.

```
cat ../data/LD/sex_specific/black.superscaffold36.boundary.female.map | cut -f2,4 > \
../data/LD/sex_specific/BB.superscaffold36.boundary.map
```
## Visualize by Haploview gui in Mac
Open and upload the ped and map file for male and female in the Haploview and save the figure as svg and compressed png file.
```
java -jar Haploview.jar 
```

## Copy results to results directory from LD directory in data
```
cp *.LDdecay.bin.gz ../../../result/LD
cp *pairwise.LD.gz ../../../result/LD
```

## Recombination rate from the genetic map

Obtain per window sex averaged, male and female recombination rates


module load bioinfo-tools BEDTools

./recombination.py ../data/linkage_map/LGZ3.sex_averaged.lifted.bed > \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.txt
sed 's/superscaffold54.1/superscaffold54/g' \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.txt > \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.54.txt
awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; \
else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.54.txt > \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.54.corrected.txt

grep 'female_input' ../data/linkage_map/LGZ3.male_female.lifted.bed > ../data/linkage_map/LGZ3.female.tmp
grep -w 'male_input' ../data/linkage_map/LGZ3.male_female.lifted.bed > ../data/linkage_map/LGZ3.male.tmp
./recombination.py ../data/linkage_map/LGZ3.female.tmp > ../data/linkage_map/LGZ3.female.tmp.rec.rate.txt
./recombination.py ../data/linkage_map/LGZ3.male.tmp > ../data/linkage_map/LGZ3.male.tmp.rec.rate.txt

rm -f ../data/linkage_map/LGZ3.female.tmp
rm -f ../data/linkage_map/LGZ3.male.tmp

sed 's/superscaffold54.1/superscaffold54/g' ../data/linkage_map/LGZ3.female.tmp.rec.rate.txt > \
../data/linkage_map/LGZ3.female.tmp.rec.rate.54.txt
sed 's/superscaffold54.1/superscaffold54/g' ../data/linkage_map/LGZ3.male.tmp.rec.rate.txt > \
../data/linkage_map/LGZ3.male.tmp.rec.rate.54.txt

awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
../data/linkage_map/LGZ3.female.tmp.rec.rate.54.txt > ../data/linkage_map/LGZ3.female.tmp.rec.rate.54.corrected.txt
awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
../data/linkage_map/LGZ3.male.tmp.rec.rate.54.txt > ../data/linkage_map/LGZ3.male.tmp.rec.rate.54.corrected.txt

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/linkage_map/ostrich.Z.${window}.bed \
-b ../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.54.corrected.txt \
-wao > ../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.${window}.txt
./get_recombination_per_window.py \
../data/linkage_map/LGZ3.sex_averaged.lifted.rec.rate.${window}.txt \
../data/linkage_map/ostrich.Z.${window}.bed > ../data/linkage_map/ostrich_Z_rec_per_${window}.txt
python recombination_window_forR.py ../data/linkage_map/ostrich_Z_rec_per_${window}.txt \
> ../data/linkage_map/ostrich_Z_rec_per_${window}.forR.txt
done

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for sex in male female
do
echo ${sex}
bedtools intersect -a ../data/linkage_map/ostrich.Z.${window}.bed \
-b ../data/linkage_map/LGZ3.${sex}.tmp.rec.rate.54.corrected.txt \
-wao > ../data/linkage_map/LGZ3.${sex}.lifted.rec.rate.window_overlap.${window}.txt
./get_recombination_per_window.py \
../data/linkage_map/LGZ3.${sex}.lifted.rec.rate.window_overlap.${window}.txt \
../data/linkage_map/ostrich.Z.${window}.bed > ../data/linkage_map/ostrich_Z_rec_${sex}_per_${window}.txt
python recombination_window_forR.py ../data/linkage_map/ostrich_Z_rec_${sex}_per_${window}.txt \
> ../data/linkage_map/ostrich_Z_rec_${sex}_per_${window}.forR.txt
done
done


# Get the R plot for sex-specific recombination rate for superscaffold36
cd /proj/snic2020-16-269/private/homap/ostrich_z/data/linkage_map
Rscript ../../bin/plot_recombination_rate.R ostrich_Z_rec_male_per_200Kb.forR.txt ostrich_Z_rec_female_per_200Kb.forR.txt
