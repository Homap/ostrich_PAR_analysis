#!/usr/bin/bash

# export black_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/black
# export blue_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/blue
# export red_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/red

#export black_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/black/P1878_110
#export blue_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/blue/P1878_121
#export red_dir=/proj/uppstore2017180/private/cornwallis/data/interim/individual/red/P1878_128 


ref='/proj/uppstore2017180/private/homap/ostrich_Z_diversity/data/reference/Struthio_camelus.20130116.OM.fa'
outdir='/proj/uppstore2017180/private/homap/ostrich_z/data/psmc'

while read -r ind
do
echo ${ind}.merge.dup.bam 
sbatch psmc_black.sh $ref \
${outdir}/${ind}.merge.dup.bam \
${outdir}/${ind}.fq.gz \
${outdir}/${ind}.fourline.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/SC.genome.allRepeats.bed \
${outdir}/${ind}.fourline_fastq_repeat_masked.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
${outdir}/${ind}.fourline_fastq_repeat_cds_masked.fq\
${outdir}/${ind}.psmcfa \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/Z_scaffolds.txt \
${outdir}/${ind}.autosome.psmcfa \
${outdir}/${ind}.autosome.psmcfa.out \
${outdir}/${ind}.autosome.psmc.image \
${outdir}/${ind}.split_autosome_psmcfa_out
done < ../data/samples/black.psmc

while read -r ind
do
echo ${ind}.merge.dup.bam 
sbatch psmc_blue.sh $ref \
${outdir}/${ind}.merge.dup.bam \
${outdir}/${ind}.fq.gz \
${outdir}/${ind}.fourline.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/SC.genome.allRepeats.bed \
${outdir}/${ind}.fourline_fastq_repeat_masked.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
${outdir}/${ind}.fourline_fastq_repeat_cds_masked.fq\
${outdir}/${ind}.psmcfa \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/Z_scaffolds.txt \
${outdir}/${ind}.autosome.psmcfa \
${outdir}/${ind}.autosome.psmcfa.out \
${outdir}/${ind}.autosome.psmc.image \
${outdir}/${ind}.split_autosome_psmcfa_out
done < ../data/samples/blue.psmc



while read -r ind
do
echo ${ind}.merge.dup.bam 
sbatch psmc_red.sh $ref \
${outdir}/${ind}.merge.dup.bam \
${outdir}/${ind}.fq.gz \
${outdir}/${ind}.fourline.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/SC.genome.allRepeats.bed \
${outdir}/${ind}.fourline_fastq_repeat_masked.fq \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/CDS_coordinates.bed \
${outdir}/${ind}.fourline_fastq_repeat_cds_masked.fq\
${outdir}/${ind}.psmcfa \
/proj/uppstore2017180/private/homap/ostrich_z/data/bed/Z_scaffolds.txt \
${outdir}/${ind}.autosome.psmcfa \
${outdir}/${ind}.autosome.psmcfa.out \
${outdir}/${ind}.autosome.psmc.image \
${outdir}/${ind}.split_autosome_psmcfa_out
done < ../data/samples/red.psmc

#${outdir}/${ind}.coverage \
#${outdir}/${ind}.mean.coverage \
#10 \

# while read -r ind
# do
# echo ${ind}.merge.dup.bam 
# sbatch psmc.sh $ref ${black_dir}/${ind}/${ind}.merge.dup.bam \
# black/${ind}.fq.gz black/${ind}.psmcfa black/${ind}.autosome.psmcfa \
# black/${ind}.autosome.psmcfa.out black/${ind}.autosome.psmc.image
# done < ../../data/samples/black.args

# while read -r ind
# do
# echo ${ind}.merge.dup.bam 
# sbatch psmc.sh $ref ${blue_dir}/${ind}/${ind}.merge.dup.bam \
# blue/${ind}.fq.gz blue/${ind}.psmcfa blue/${ind}.autosome.psmcfa \
# blue/${ind}.autosome.psmcfa.out blue/${ind}.autosome.psmc.image
# done < ../../data/samples/blue.args

# while read -r ind
# do
# echo ${ind}.merge.dup.bam 
# sbatch psmc.sh $ref ${red_dir}/${ind}/${ind}.merge.dup.bam \
# red/${ind}.fq.gz red/${ind}.psmcfa red/${ind}.autosome.psmcfa \
# red/${ind}.autosome.psmcfa.out red/${ind}.autosome.psmc.image
# done < ../../data/samples/red.args
