# Male - Female Fst

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
