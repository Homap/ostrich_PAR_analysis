#_____________________________________________________________________________________________#
# B. GFF
#_____________________________________________________________________________________________#
# 1. Get basic statistics for GFF
# 2. Check overlapping transcripts in GFF, get the longest transcript
# 3. Check protein predictions, if they correspond well to chicken when blasted.
# 4. Check if the GFF correspond well with the OP assembly before and after linkage mapping.
# 5. Get gene density per window in bed format.
# 6. Get coordinates of exons, introns and intergenic regions.
# 7. Extract transcripts per individual using the SNPs in vcf to use to create alignments.
# Get some statistics of the GFF files:
# How many line are there?
zcat gff/Struthio_camelus.OM.gene.20130116.gff.gz | wc -l
# 145691
zcat ortholog/palaeo_genomes/ostrich/strcam_assembly/GCF_000698965.1_ASM69896v1_genomic.gff.gz | wc -l
# 739459
# How many genes are there?
zcat gff/Struthio_camelus.OM.gene.20130116.gff.gz | awk '{if($3=="gene") print $0}' | wc -l
# 0: no feature called gene
zcat ortholog/palaeo_genomes/ostrich/strcam_assembly/GCF_000698965.1_ASM69896v1_genomic.gff.gz | \
awk '{if($3=="gene") print $0}' | wc -l
# 16547
# How many mRNA are there?
zcat gff/Struthio_camelus.OM.gene.20130116.gff.gz | awk '{if($3=="mRNA") print $0}' | wc -l
# 16178
zcat ortholog/palaeo_genomes/ostrich/strcam_assembly/GCF_000698965.1_ASM69896v1_genomic.gff.gz | \
awk '{if($3=="mRNA") print $0}' | wc -l
# 25425
# How many types of feature (3rd column) are there?
zcat gff/Struthio_camelus.OM.gene.20130116.gff.gz | awk '{if($0 !~ /^#/)print $3}' | sort -u
# CDS, mRNA
zcat ortholog/palaeo_genomes/ostrich/strcam_assembly/GCF_000698965.1_ASM69896v1_genomic.gff.gz | \
awk '{if($0 !~ /^#/)print $3}' | sort -u
# cDNA_match, CDS, C_gene_segment, exon, gene, lnc_RNA, mRNA, pseudogene, region, rRNA, sequence_feature, 
# transcript, tRNA, V_gene_segment

# Do a lift-over between ostrich old and new assembly 
# This is done due to errors in Struthio_camelus.OM.gene.20130116.gff.gz.
# Liftover is done using the program flo downloaded into src:
wget -c https://github.com/yeban/flo/archive/master.tar.gz -O flo.tar.gz
tar xvf flo.tar.gz
mv flo-master flo

mkdir data/gff/gff_liftover 
cd data/gff/gff_liftover 

cp ../../../src/flo/opts_example.yaml flo_opts.yaml  
../../../src/flo/scripts/install.sh

# flo can only work with transcripts and their child exons and CDS. 
# Transcripts can be annotated as: mRNA, transcript, or gene. 
# However, if you have a 'gene' annotation for each transcript, 
# you will need to remove that: 
module load bioinfo-tools ruby 

../src/flo/gff_remove_feats.rb gene \
../data/ortholog/palaeo_genomes/ostrich/strcam_assembly/GCF_000698965.1_ASM69896v1_genomic.gff \
> ../data/gff/gff_liftover/GCF_000698965.1_ASM69896v1_genomic.flo.gff

awk '{if($3=="mRNA") print $0}' ../data/gff/gff_liftover/GCF_000698965.1_ASM69896v1_genomic.flo.gff | \
wc -l
# 25425

# Now edit flo_opts.yaml to indicate:  
# Location of source and target assembly in FASTA format (required).
# Location of GFF3 file(s) containing annotations on the source assembly. 
# Number of CPU cores to use (required - not auto detected). This
# cannot be greater than the number of scaffolds in the target assembly.
# In /proj/uppstore2017180/private/homap/ostrich_Z_diversity/data/gff/gff_liftover:
sbatch liftover.sh

awk '{if($3=="mRNA") print $0}' lifted.gff3 | wc -l
# 24547
awk '{if($0 !~ /^#/)print $3}' lifted.gff3  | sort -u
# cDNA_match, CDS, C_gene_segment, exon, lnc_RNA, mRNA, pseudogene
# region, transcript, tRNA, V_gene_segment

# There are many features in the gff file but we would like to keep
# only mRNA, exon, CDS.
# Extract the CDS: Check if they are fine.

# Get overlapping mRNA and choose the longest transcript

python clean_gff.py ../data/gff/gff_liftover/lifted.gff3 > ../data/gff/CDS_coordinates.txt
python extract_CDS_from_fasta.py ../data/gff/CDS_coordinates.txt \
../data/reference/Struthio_camelus.20130116.OM.fa > ../data/gff/ostrich_CDS.fa
grep '>' ../data/gff/ostrich_CDS.fa | wc 
# 20078
# Summary of gene length data
"""
Min.   :    96.0  
1st Qu.:   843.8  
Median :  1365.0  
Mean   :  1812.3  
3rd Qu.:  2246.2  
Max.   :103571.0 
"""

python2 translate/translateone.py ../data/gff/ostrich_CDS.fa > \
../data/gff/ostrich_protein.fa 2> ../data/gff/ostrich_CDS_translation_AmbigLongORF
# Remove two sequences that were turned into '---' upon translation
grep -v '>scaffold1200_rna-XM_009673897.1_XP_009672192.1'  ostrich_protein.fa \
| grep -v '\-\-\-' | grep -v 'scaffold987_rna-XM_009681896.1_XP_009680191.1' > ostrich_protein_clean.fa
# Total number of genes
grep '>' ostrich_protein_clean.fa | wc                                                                               
# 20076
# Gene starting with Methionine
grep '>' -A1 ostrich_protein_clean.fa | grep -v '>' | grep '^M' | wc                                                 
# 18548

# Blast against chicken data
bash blast.sh


# Annotation

# Two ways:

# Using cleaned gff with snpEff

# Custom script


# Create snpEff database for ostrich Z chromosome
# Download and install SnpEff
curl -v -L http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip > snpEff_latest_core.zip
unzip snpEff_latest_core.zip

#How to create a snpEff database using a gff3 and genomic DNA fasta file... (note, the chromosome names must match in the 2 files)
#NOTE: This uses /bin/tcsh...

export DBNAME=ostrich_assembly
export GFF3=/proj/uppstore2017180/private/homap/ostrich_Z_diversity/data/gff/gff_liftover/lifted.gff3
export FASTA=/proj/uppstore2017180/private/homap/ostrich_Z_diversity/data/reference/Struthio_camelus.20130116.OM.fa

# Download in /proj/snic2020-16-269/private/homap/ostrich_z/bin/
curl -v -L http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip > snpEff_latest_core.zip
unzip snpEff_latest_core.zip

export DBNAME=ostrich_snpeff
export GFF3=/proj/snic2020-16-269/private/homap/ostrich_z/data/gff/gff_liftover/lifted.gff3
export FASTA=/proj/snic2020-16-269/private/homap/ostrich_z/data/reference/Struthio_camelus.20130116.OM.fa
#Go into the snpEff directory and create a directory for your files
cd src/snpEff
mkdir -p $DBNAME

#Copy the files into snpEff's directory structure
cp $GFF3 $DBNAME
cp $FASTA $DBNAME

#Edit snpEff.config and insert your specific database information:
echo "$DBNAME.genome : $DBNAME" >> ostrich.config

#Build the database
java -jar snpEff.jar build -gff3 -v $DBNAME -c /proj/snic2020-16-269/private/homap/ostrich_z/bin/snpEff/ostrich.config
java -jar snpEff/snpEff.jar eff chr18 freebayes.chr18.all.vcf > snpeff.chr18.vcf

# Get coordinates of exons, introns and intergenic regions

# Annotate all coding sites in terms of synonymous, nonsynonymous, etc.
grep -v '^C' genes.gff > genes.clean.gff
mv genes.c


#!/bin/bash -l

#SBATCH -A snic2021-22-447
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J snpEff
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL
java -jar snpEff.jar build -gff3 -v ostrich_snpeff -c /proj/snic2020-16-269/private/homap/ostrich_z/bin/snpEff/ostrich.config -dataDir /proj/snic2020-16-269/private/homap/ostrich_z/bin/snpEff

# In Z linked genes:
grep -E 'superscaffold26|superscaffold54|superscaffold35|superscaffold36|superscaffold62|superscaffold67|superscaffold69-1|superscaffold93|superscaffold63|superscaffold88|superscaffold83|superscaffold92' ../gff/Struthio_camelus.OM.gene.20130116.gff > Z.gff 

awk '$3=="mRNA' Z.gff
# In z_linked_genes
python ../../bin/scaffold_to_chr_vcf.py Z.genes.coordinates > Z.genes.chr.coordinates