# Measures of genetic variation

Using SNPs, we calculate different measures of genetic variation:

- Average number of pairwise differences (π)
- Average number of segregating sites (θ)
- Tajima's D which reflects the relationship between π and θ
- Folded site frequency spectrum

## Getting allele counts from scaffolds of chromosomes 4, 5, PAR and nonPAR

`bash 1_get_allele_count.sh`

## Calculate measures of genetic variation

`bash 2_get_genetic_variation.sh`

# Obtain SFS of synonymous and nonsynonymous sites
"""
I want to calculate the site frequency spectrum of the zerofold and fourfold sites of the PAR, nonPAR and chromosome 4 and 5 of autosomes.
This analyses will give me three figures each containing the SFS of 4fold, 0fold and neutral equilibrium for autosome, PAR and nonPAR.
I will use the 4fold sites to infer effects of population size changes and 0fold sites to infer selection using DFE-alpha. The first step is 
to find out which sites are zero fold and which sites are 4 fold. I use a program in Python that I have written for this purpose.
"""
"""
In the first step, I annotate all coding sites in the genome to see whether they are synonymous, nonsynonymous, and among those 4fold or 0fold or nonsense
which means the change in the nucleotide has led to stop codon.
"""
# In the folder, annotate_old_2, we do the following:
cd annotate_old_2
# Annotate all genomic sites for whether they are zero fold and fourfold sites
# I added a line to print out the transcripts that have stop codons
python annotate_sites.py ../../data/reference/Struthio_camelus.20130116.OM.fa ../../data/gff/Struthio_camelus.OM.gene.20130116.gff > ../../data/gff/Struthio_camelus.OM.gene.20130116.annotated.txt


# Diversity in 4-fold and 0-fold sites
python count_total_num_syn_nonsyn.py ../data/gff/Struthio_camelus.OM.gene.20130116.annotated.txt \
../data/allele_count/black.A.repeat.frq.count 20 > ../data/allele_count/black.A.Struthio_camelus.OM.gene.20130116.annotated.pi.txt

python fourfold_zerofold_diversity.py ../data/allele_count/black.A.Struthio_camelus.OM.gene.20130116.annotated.pi.txt > ../result/sfs_measures/black.A.Struthio_camelus.OM.gene.20130116.pnps.txt

python count_total_num_syn_nonsyn.py ../data/gff/genome.annotated.txt \
../data/allele_count/black.PAR.repeat.frq.count 20 > ../data/allele_count/black.PAR.annotated.pi.txt
python fourfold_zerofold_diversity.py ../data/allele_count/black.PAR.annotated.pi.txt > ../result/sfs_measures/black.PAR.pnps.txt

python count_total_num_syn_nonsyn.py ../data/gff/genome.annotated.txt \
../data/allele_count/black.nonPAR.filtered.adjusted.frq.count 15 > ../data/allele_count/black.nonPAR.annotated.pi.txt
python fourfold_zerofold_diversity.py ../data/allele_count/black.nonPAR.annotated.pi.txt > ../result/sfs_measures/black.nonPAR.pnps.txt

# Add number of sites
python fourfold_zerofold_diversity.py ../data/allele_count/black.PAR.annotated.pi.txt > ../result/sfs_measures/black.PAR.pnps.number_syn_numer_nonsyn.txt

for species in black blue red
do 
echo $species
echo "PAR"
python get_sfs.py ../data/allele_count/${species}.PAR.repeat.frq.count ${species} PAR > ../result/sfs_measures/${species}_PAR.foldedSFS.txt
echo "nonPAR"
python get_sfs.py ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count ${species} nonPAR > ../result/sfs_measures/${species}_nonPAR.foldedSFS.txt
echo "Autosome"
python get_sfs.py ../data/allele_count/${species}.A.repeat.frq.count ${species} A > ../result/sfs_measures/${species}_A.foldedSFS.txt
done
#*******

# To calculate SFS for 4fold and 0fold, all I need is to know which SNP is 0fold and 4fold and overlap it with the allele count data
python count_0fold_4fold_overlap.py ../data/gff/Struthio_camelus.OM.gene.20130116.annotated.txt ../data/allele_count/black.A.repeat.frq.count 20

python count_0fold_4fold_overlap.py ../data/gff/Struthio_camelus.OM.gene.20130116.annotated.txt ../data/allele_count/black.PAR.repeat.frq.count 20

python count_0fold_4fold_overlap.py ../data/gff/Struthio_camelus.OM.gene.20130116.annotated.txt ../data/allele_count/black.nonPAR.repeat.frq.count 15

python get_sfs.py fourfold_counts.txt black PAR > fourfold_PAR_SFS.txt
python get_sfs.py zerofold_counts.txt black PAR > zerofold_PAR_SFS.txt

python get_sfs.py fourfold_A_counts.txt black PAR > fourfold_A_SFS.txt
python get_sfs.py zerofold_A_counts.txt black PAR > zerofold_A_SFS.txt

python get_sfs.py fourfold_A_counts.txt black nonPAR > fourfold_nonPAR_SFS.txt
python get_sfs.py zerofold_A_counts.txt black nonPAR > zerofold_nonPAR_SFS.txt

# 4fold and 0fold
1
11
3009164 4298    2812    2016    1529    1264    1119    969 825 883 542
11477529 4235    2530    1709    1261    974 783 640 563 525 343


### Get the genetic diversity for superscaffold36 in windows of 1Kb and 0.5 Kb
```
grep 'superscaffold36' ../data/allele_count/black.PAR.repeat.frq.count > ../data/allele_count/black.PAR.superscaffold36.repeat.frq.count
grep 'superscaffold36' ../data/allele_count/black.nonPAR.repeat.frq.count > ../data/allele_count/black.nonPAR.superscaffold36.repeat.frq.count
python SFS_measures.py ../data/allele_count/black.PAR.superscaffold36.repeat.frq.count 20 ../data/bed/black.par_scaf.10Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.scaf36.10Kb.sfs.txt
python SFS_measures.py ../data/allele_count/black.nonPAR.superscaffold36.repeat.frq.count 15 \
../data/bed/black.nonpar_scaf.10Kb.intergenic.overlap.density.sorted.txt all > ${result}/black.nonPAR.scaf36.10Kb.sfs.txt

python SFS_measures.py ../data/allele_count/black.PAR.superscaffold36.repeat.frq.count 20 ../data/bed/black.par_scaf.1Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.scaf36.1Kb.sfs.txt
python SFS_measures.py ../data/allele_count/black.PAR.superscaffold36.repeat.frq.count 20 ../data/bed/black.par_scaf.0.5Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.scaf36.0.5Kb.sfs.txt
```
