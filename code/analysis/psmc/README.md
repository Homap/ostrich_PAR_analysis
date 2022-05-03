# PSMC analysis

bash psmc_run.sh

python mask_fastq.py /proj/snic2020-16-269/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
../data/psmc/P1878_110.fourline_fastq_repeat_masked.fq > ../data/psmc/P1878_110.fourline_fastq_repeat_cds_masked.fq

python mask_fastq.py /proj/snic2020-16-269/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
../data/psmc/P1878_121.fourline_fastq_repeat_masked.fq > ../data/psmc/P1878_121.fourline_fastq_repeat_cds_masked.fq

python mask_fastq.py /proj/snic2020-16-269/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
../data/psmc/P1878_128.fourline_fastq_repeat_masked.fq > ../data/psmc/P1878_128.fourline_fastq_repeat_cds_masked.fq

fq2psmcfa -q20 ../data/psmc/P1878_110.fourline_fastq_repeat_cds_masked.fq > ../data/psmc/P1878_110.psmcfa
fq2psmcfa -q20 ../data/psmc/P1878_121.fourline_fastq_repeat_cds_masked.fq > ../data/psmc/P1878_121.psmcfa
fq2psmcfa -q20 ../data/psmc/P1878_128.fourline_fastq_repeat_cds_masked.fq > ../data/psmc/P1878_128.psmcfa

python extract_seq_from_fasta.py ../data/psmc/P1878_110.psmcfa list \
../data/bed/Z_scaffolds.txt > ../data/psmc/P1878_110.autosome.psmcfa

python extract_seq_from_fasta.py ../data/psmc/P1878_121.psmcfa list \
../data/bed/Z_scaffolds.txt > ../data/psmc/P1878_121.autosome.psmcfa

python extract_seq_from_fasta.py ../data/psmc/P1878_128.psmcfa list \
../data/bed/Z_scaffolds.txt > ../data/psmc/P1878_128.autosome.psmcfa

sbatch psmc_black.sh
sbatch psmc_blue.sh
sbatch psmc_red.sh

# Obtain 1000 bootstrap of PSMC to calculate confidence interval
# The black subspecies
for i in {1..100}
do
mkdir -p ../data/psmc/black.CI.${i}.dir
sbatch psmc_CI.sh ../data/psmc/black.CI.${i}.dir ../data/psmc/P1878_110.autosome.split.psmcfa
done

cd ../data/psmc

# for i in black.CI.*
# do
# echo ${i}/round_repeat-1.psmc
# cat ${i}/round_repeat-1.psmc >> black.CI.resampling.txt
# done
# cat P1878_110.autosome.psmc.out black.CI.resampling.txt > temp && mv temp black.CI.resampling.txt

mkdir -p black_resampling_CI

j=0
for i in black.CI.*
do
echo ${i}/round_repeat-1.psmc
((j=j+1))
echo $j
cp ${i}/round_repeat-1.psmc black_resampling_CI
mv black_resampling_CI/round_repeat-1.psmc black_resampling_CI/round_repeat-${j}.psmc
done
cp black.CI.1.dir/round_repeat-1.psmc black_resampling_CI
