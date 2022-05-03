# Steps for obtaining genetic diversity (pi and theta) and folded SFS together with male-female Fst

```
export result=/proj/snic2020-16-269/private/homap/ostrich_z/result/sfs_measures
```

## Measure of genetic diversity (pi and theta) and folded SFS

### Autosomes
```
python SFS_measures_autosome.py ../data/allele_count/black.A.repeat.frq.count 20 ../data/bed/black.autosome.100Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.autosome.sfs.txt
```

### PAR and nonPAR
```
python SFS_measures.py ../data/allele_count/black.PAR.repeat.frq.count 20 ../data/bed/black.par_scaf.100Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.sfs.txt
python SFS_measures.py ../data/allele_count/black.PAR.repeat.frq.count 20 ../data/bed/black.par_scaf.500Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.500Kb.sfs.txt
python SFS_measures.py ../data/allele_count/black.PAR.repeat.frq.count 20 ../data/bed/black.par_scaf.1000Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.PAR.1000Kb.sfs.txt

python SFS_measures.py ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count 15 \
../data/bed/black.nonpar_scaf.100Kb.intergenic.overlap.density.sorted.txt all > ${result}/black.nonPAR.sfs.txt
python SFS_measures.py ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count 15 \
../data/bed/black.nonpar_scaf.500Kb.intergenic.overlap.density.sorted.txt all > ${result}/black.nonPAR.500Kb.sfs.txt
python SFS_measures.py ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count 15 \
../data/bed/black.nonpar_scaf.1000Kb.intergenic.overlap.density.sorted.txt all > ${result}/black.nonPAR.1000Kb.sfs.txt
```

### Diversity measures for whole Z
```
awk 'NR>1' ${result}/black.nonPAR.sfs.txt | cat ${result}/black.PAR.sfs.txt - > ${result}/black.Z.sfs.txt
```

### Translate scafoold to chromosome Z coordinates
```
python scaffold_to_chr_vcf.py ${result}/black.Z.sfs.txt > ${result}/black.Z.coordinates.sfs.txt
```

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

### Male-Female FST
```
vcftools --gzvcf ../data/vcf/black.PAR.filtered.vcf.gz --chr 'superscaffold36' --recode --stdout | gzip -c > ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz
vcftools --gzvcf ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz --weir-fst-pop ../data/samples/black_male.txt --weir-fst-pop ../data/samples/black_female.txt \
--fst-window-size 1000 --out black_male_female_1Kb
```

### Male-Female FST
# Exclude female positions with heterozygous SNPs
```
awk '{print $1"\t"$3}' ../data/bed/${species}.nonPAR.female_het.bed > ../data/bed/${species}.nonPAR.female_het.txt
vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.pos.vcf.gz --exclude-positions ../data/bed/${species}.nonPAR.female_het.txt --stdout | gzip -c > ../data/vcf/LD_vcf/black.Z.pos.nohetfemalenonPAR.vcf.gz

vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.pos.nohetfemalenonPAR.vcf.gz --weir-fst-pop ../data/samples/black_male.txt --weir-fst-pop ../data/samples/black_female.txt --fst-window-size 100000 --out ../data/FST/black_male_female_Z_100Kb
python scaffold_to_chr_vcf.py ../data/FST/black_male_female_Z_100Kb.windowed.weir.fst > ../data/FST/black_male_female_Z_100Kb.windowed.weir.Z.coord.fst
```

### Get the distribution of SFS for autosomes, PAR and nonPAR
```
echo "PAR"
python get_sfs.py ../data/allele_count/${species}.PAR.repeat.frq.count ${species} PAR > ../result/sfs_measures/${species}_PAR.foldedSFS.txt
echo "nonPAR"
python get_sfs.py ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count ${species} nonPAR > ../result/sfs_measures/${species}_nonPAR.foldedSFS.txt
echo "Autosome"
python get_sfs.py ../data/allele_count/${species}.A.repeat.frq.count ${species} A > ../result/sfs_measures/${species}_A.foldedSFS.txt
```
