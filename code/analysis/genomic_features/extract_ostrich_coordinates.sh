#!/usr/bin/bash

module load bioinfo-tools BEDTools/2.29.2

# Correcting GFF of ostrich genome assembly
# Obtaining coordinates of intergenic, intronic and CDS regions of GFF in BED format
#*************************************************************************************************************************************************************************************************************************
echo "Reformatting GFF"
./reformat_gff.py ../../../data/genome/Struthio_camelus.OM.gene.20130116.gff > ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.gff
#*************************************************************************************************************************************************************************************************************************
# Remove overlapping mRNA
#*************************************************************************************************************************************************************************************************************************
echo "Removing overlapping mRNA"
awk '$3=="mRNA"' ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.gff | awk '{print $1"\t"$4-1"\t"$5}' | sort -k1,1 -k2,2n | \
bedtools merge -i - -c 1 -o count | awk '{if($4==1) print $0}' | bedtools intersect -a ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.gff -b - > ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff
#*************************************************************************************************************************************************************************************************************************
# Extract CDS coordinates in bed format
#*************************************************************************************************************************************************************************************************************************
echo "Extract CDS coordinates"
./extract_from_gff.py ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff CDS all > ../../../data/genomic_features/ostrich.all.CDS.coord
# Choose Z and chr4 and chr5 genes
cut -f1 ../../../data/bed/z_scaf.bed | grep -f - ../../../data/genomic_features/ostrich.all.CDS.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/Z.CDS.txt
cut -f1 ../../../data/lastz/gg_chr4_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.all.CDS.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr4.CDS.txt
cut -f1 ../../../data/lastz/gg_chr5_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.all.CDS.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr5.CDS.txt
#*************************************************************************************************************************************************************************************************************************
# Extract intronic coordinates
#*************************************************************************************************************************************************************************************************************************
echo "Extract intronic coordinates"
./extract_from_gff.py ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff intron all > ../../../data/genomic_features/ostrich.all.intron.coord
# Choose Z and chr4 and chr5 genes
cut -f1 ../../../data/bed/z_scaf.bed | grep -f - ../../../data/genomic_features/ostrich.all.intron.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/Z.intron.txt
cut -f1 ../../../data/lastz/gg_chr4_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.all.intron.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr4.intron.txt
cut -f1 ../../../data/lastz/gg_chr5_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.all.intron.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr5.intron.txt
#*************************************************************************************************************************************************************************************************************************
# Extract intergenic coordinates
#*************************************************************************************************************************************************************************************************************************
echo "Extract intergenic coordinates"
awk '{print $1"\t"$2}' ../../../data/genome/Struthio_camelus.20130116.OM.fa.fai | grep -v "^C" | sort -k1,1 > ../../../data/genome/Struthio_camelus.20130116.OM.length.txt
awk '$3=="mRNA"' ../../../data/genome/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff | awk '{print $1"\t"$4-1"\t"$5}' | sort -k1,1 -k2,2n | grep -v "^C" | \
bedtools merge -i - -c 1 -o count | bedtools complement -i - -g ../../../data/genome/Struthio_camelus.20130116.OM.length.txt > ../../../data/genomic_features/ostrich.intergene.coord
# Choose Z and chr4 and chr5 genes
cut -f1 ../../../data/bed/z_scaf.bed | grep -f - ../../../data/genomic_features/ostrich.intergene.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/Z.intergene.txt
cut -f1 ../../../data/lastz/gg_chr4_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.intergene.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr4.intergene.txt
cut -f1 ../../../data/lastz/gg_chr5_ostrich.bed | grep -f - ../../../data/genomic_features/ostrich.intergene.coord | sort -k1,1 -k2,2n > ../../../data/genomic_features/chr5.intergene.txt