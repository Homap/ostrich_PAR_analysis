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
