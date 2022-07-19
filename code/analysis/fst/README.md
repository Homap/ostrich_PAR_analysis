# Male - Female Fst

`module load bioinfo-tools vcftools/0.1.16`

```
mkdir -p ../../../data/fst/z
mkdir -p ../../../data/fst/autosome/chr4
mkdir -p ../../../data/fst/autosome/chr5
```

```
vcftools --gzvcf ../../../data/vcf/z_vcf/sorted.z.vcf.gz --weir-fst-pop ../../../data/samples/black_male.txt --weir-fst-pop ../../../data/samples/black_female.txt --fst-window-size 200000 --out ../../../data/fst/z/male_female_Z_200Kb

python ../../processing/scaffold_to_chr.py ../../../data/fst/z/male_female_Z_200Kb.windowed.weir.fst Z > ../../../data/fst/z/black_male_female_Z_200Kb.windowed.weir.Z.coord.fst
```

```
for scaffold_vcf in ${a_vcf}/chr4/*chr4.vcf.gz
do
echo $scaffold_vcf
scaffold=$(echo $scaffold_vcf | cut -f9 -d "/" | cut -f1 -d ".")
vcftools --gzvcf ../../../data/vcf/ld_vcf/a_vcf/chr4/${scaffold}.chr4.sorted.gz --weir-fst-pop ../../../data/samples/black_male.txt --weir-fst-pop ../../../data/samples/black_female.txt --fst-window-size 200000 --out ../../../data/fst/autosome/chr4/${scaffold}_chr4_200Kb
done

cat ../../../data/fst/autosome/chr4/*200Kb.windowed.weir.fst | grep -v 'CHROM' > ../../../data/fst/autosome/chr4/chr4.200Kb.windowed.weir.fst
```

```
for scaffold_vcf in ${a_vcf}/chr5/*chr5.vcf.gz
do
echo $scaffold_vcf
scaffold=$(echo $scaffold_vcf | cut -f9 -d "/" | cut -f1 -d ".")
vcftools --gzvcf ../../../data/vcf/ld_vcf/a_vcf/chr5/${scaffold}.chr5.sorted.gz --weir-fst-pop ../../../data/samples/black_male.txt --weir-fst-pop ../../../data/samples/black_female.txt --fst-window-size 200000 --out ../../../data/fst/autosome/chr5/${scaffold}_chr5_200Kb
done

cat ../../../data/fst/autosome/chr5/*200Kb.windowed.weir.fst | grep -v 'CHROM' > ../../../data/fst/autosome/chr5/chr5.200Kb.windowed.weir.fst
```

```
vcftools --gzvcf ../data/vcf/black.PAR.filtered.vcf.gz --chr 'superscaffold36' --recode --stdout | gzip -c > ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz
vcftools --gzvcf ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz --weir-fst-pop ../data/samples/black_male.txt --weir-fst-pop ../data/samples/black_female.txt \
--fst-window-size 1000 --out black_male_female_1Kb
```


