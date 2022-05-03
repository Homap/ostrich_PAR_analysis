#!/bin/bash -l

#SBATCH -A snic2019-3-17 
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J psmc
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools samtools bcftools psmc seqtk perl_modules

echo "Run PSMC"
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ../data/psmc/P1878_110.autosome.psmc.out ../data/psmc/P1878_110.autosome.psmcfa

echo "Create the PSMC plot"
psmc_plot.pl -g 4 -u 4.63e-09 ../data/psmc/P1878_110.autosome.psmc.image ../data/psmc/P1878_110.autosome.psmc.out

echo "Split fasta for resampling"
splitfa ../data/psmc/P1878_110.autosome.psmcfa > ../data/psmc/P1878_110.autosome.split.psmcfa


# ref=$1
# bamfile=$2
# #coverage=$3
# #meancov=$4
# fastq=$3
# fourline_fastq=$4
# repeat_bed=$5
# fourline_fastq_repeat_masked=$6
# cds_bed=$7
# fourline_fastq_repeat_cds_masked=$8
# psmcfa=$9
# zscaffold=$10
# autosome_psmcfa=$11
# autosome_psmcfa_out=$12
# psmc_image=$13
# split_autosome_psmcfa_out=$14

#echo "Now, calculating coverage"
#samtools depth --reference $ref $bamfile > $coverage
#awk '{sum += $3} END {print sum / NR}' $coverage > $meancov
#meancov=`cat $meancov`

#echo "The mean coverage of reads is $meancov"
#d=`cat $meancov|awk '{print $0/3}'` # minimum coverage
#D=`cat $meancov|awk '{print $0*2}'` # maximum coverage

# -q, -min-MQ INT #Minimum mapping quality for an alignment to be used [0]
# -Q, --min-BQ INT #Minimum base quality for a base to be considered [13]

# echo "samtools mpileup running"
# samtools mpileup -Q20 -q20 -C50 -uf $ref $bamfile | bcftools call --threads 10 -c - | \
# vcfutils.pl vcf2fq -Q 25 -d 14 -D 82 | gzip > $fastq

# echo "Covert multi line fastq into 4-line fastq"
# seqtk seq -l0 $fastq > $fourline_fastq
# echo "Mask repeats"
# python mask_fastq.py $repeat_bed $fourline_fastq > $fourline_fastq_repeat_masked
# echo "Mask CDS"
# python mask_fastq.py $cds_bed $fourline_fastq_repeat_masked > $fourline_fastq_repeat_cds_masked

# echo "creating psmc fasta"
# fq2psmcfa -q20 $fourline_fastq_repeat_cds_masked > $psmcfa

# echo "Remove Z chromosome scaffolds"
# python extract_seq_from_fasta.py $psmcfa list $zscaffold > $autosome_psmcfa

# echo "Run PSMC"
# psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $autosome_psmcfa_out $autosome_psmcfa

# echo "Create the PSMC plot"
# psmc_plot.pl -g 4 -u 4.63e-09 $psmc_image $autosome_psmcfa_out

# echo "Split fasta for resampling"
# splitfa $autosome_psmcfa_out > $split_autosome_psmcfa_out









