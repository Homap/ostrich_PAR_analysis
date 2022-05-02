# Evolutionary dynamics of ancient recombining sex chromosomes in ostrich  

This document contains the analyses steps to reproduce the analyses 
including criteria to filter the VCF files, calculate genetic diversity (pi and theta) and 
folded site frequency spectrum in genomic categories: intergenic, intronic and CDS (0fold and 4fold) and DFE-alpha, 
calculates linkage disequilibrium (LD), population scaled recombination rate (rho) and PSMC.

## Select SNPs from the GVCF: GVCF called black.all.vcf.gz has been copied from 
```
/proj/snic2020-16-269/private/cornwallis.2020/results/ind/analysis/gatk_best_practice_snp/f03_concat_bcftools produced by Per Unneberg

sbatch selectvariant_black.sh
```

## Variant Filteration
```posh
vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf
sbatch variantfilter.sh ${vcf}/black.all.snp.vcf.gz ${vcf}/black.all.snp.filtered.vcf.gz
```

## Load the necessary modules
`module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07`

## Export reference, vcf file and the sample information

```
export reference=/proj/snic2020-16-269/private/homap/ostrich_z/data/reference/Struthio_camelus.20130116.OM.fa
export vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf
export black=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black.txt
export black_male=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black_male.txt
export black_female=/proj/snic2020-16-269/private/homap/ostrich_z/data/samples/black_female.txt
```

## Create a bed file with PAR and non-PAR regions and autosomes
```
python scaf_length_to_bed.py $reference ../data/bed/z_scaf.bed \
../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed
grep -v -f ../data/bed/Z_scaffolds.txt ../data/reference/Struthio_camelus.20130116.OM.fa.fai \
| awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $1"\t""0""\t"$2}' > ../data/bed/autosomes.bed
```

## Calculate sequencing statistics for the data
```
species=black
# Variant files - autosomes
sbatch variant_statistics_autosomes.sh ../data/vcf/${species}.A.vcf.gz ../data/vcf/${species}.A.subset.vcf ../data/vcf/${species}.A.subset
# Variant files - PAR
sbatch variant_statistics_Z.sh ../data/vcf/${species}.PAR.vcf.gz ../data/vcf/${species}.PAR
# Variant files - nonPAR
sbatch variant_statistics_Z.sh ../data/vcf/${species}.nonPAR.vcf.gz ../data/vcf/${species}.nonPAR
```

I copy the outputs in the local computer in (/Users/homapapoli/Documents/projects/sex_chr/ostrich_z/data/vcf)
and check the plots in RStudio. Based on the outputs, decide for filtering and apply below.

```
for part in PAR nonPAR
do
echo $part
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.frq .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.het .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.idepth .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.imiss .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.ldepth.mean .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.lmiss .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.lqual .
done

scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.frq .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.het .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.idepth .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.imiss .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.ldepth.mean .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.lmiss .
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.A.subset.lqual .
```

## Perform filtering with vcftools

I set the minimum coverage for all to 10 reads per site

## Autosome mean coverage is about 32X, therefore I set the max depth to 64X
```
sbatch vcftools_filter.sh ../data/vcf/${species}.A.vcf.gz ../data/vcf/${species}.A.filtered.vcf.gz 64
```
## PAR mean coverage is about 30X, I set the max depth to 60X
```
sbatch vcftools_filter.sh ../data/vcf/${species}.PAR.vcf.gz ../data/vcf/${species}.PAR.filtered.vcf.gz 60
```
## nonPAR mean coverage is about 24X, I set the max depth to 50X
```
sbatch vcftools_filter.sh ../data/vcf/${species}.nonPAR.vcf.gz ../data/vcf/${species}.nonPAR.filtered.vcf.gz 50
````

## Number of variants after filtering
```
bcftools view -H ../data/vcf/${species}.A.filtered.vcf.gz | wc -l
bcftools view -H ../data/vcf/${species}.PAR.filtered.vcf.gz | wc -l
bcftools view -H ../data/vcf/${species}.nonPAR.filtered.vcf.gz | wc -l
```

| Category | Number of SNPs |
| :--------- | :-----: |
| Autosome   | 6620899 |
| PAR  | 306392 |
| nonPAR | 59751 |


## VCF stats depth and quality after filtering
```
cd /proj/snic2020-16-269/private/homap/ostrich_z/data/vcf run the following
vcftools --gzvcf ${species}.A.filtered.vcf.gz --site-mean-depth --out ${species}.A.filtered
vcftools --gzvcf ${species}.PAR.filtered.vcf.gz --site-mean-depth --out ${species}.PAR.filtered
vcftools --gzvcf ${species}.nonPAR.filtered.vcf.gz --site-mean-depth --out ${species}.nonPAR.filtered
```

## Copy the site mean depth into local computer under /Users/homapapoli/Documents/projects/sex_chr/ostrich_z/data/vcf
```
for part in A PAR nonPAR
do
echo $part
scp  homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/${species}.${part}.filtered.ldepth.mean .
done
```

## Background Filtering for Coverage
```
sbatch get_depth.sh ../data/vcf/${species}.all.vcf.gz ../data/vcf/${species}.all
```
# Get a list of sites with coverage below or above the threshold used for SNP filtering
```
sbatch get_background_coverage.sh ../data/vcf/${species}.all.ldepth.mean ../data/vcf/${species}.all.ldepth.mean.filter.bed
```

## Obtain data for the plot of boundary coverage
```
mkdir -p ../result/boundary_coverage
module load bioinfo-tools samtools
bamaddr=../../../cornwallis.2020/data/interim/map/bwa/dedup 
samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_107.bam ${bamaddr}/P1878_108.bam ${bamaddr}/P1878_109.bam ${bamaddr}/P1878_110.bam ${bamaddr}/P1878_111.bam \
> ../result/boundary_coverage/male_boundary_coverage.txt

samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_112.bam ${bamaddr}/P1878_113.bam ${bamaddr}/P1878_114.bam ${bamaddr}/P1878_115.bam ${bamaddr}/P1878_116.bam \
> ../result/boundary_coverage/female_boundary_coverage.txt

scp homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/result/boundary_coverage/male_boundary_coverage.txt
scp homap@rackham.uppmax.uu.se:/proj/snic2020-16-269/private/homap/ostrich_z/result/boundary_coverage/female_boundary_coverage.txt
```

Use the repeat masked genome where repeats are already N. For analysis where you need the number 
of filtered invariant sites, you can simply count non-N sites in each window and that will be 
the number of bases in window.

```
module load bioinfo-tools BEDTools
bedtools maskfasta -fi ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked \
-bed ../data/vcf/${species}.all.ldepth.mean.filter.bed -fo ../data/reference/${species}.repeat.depth.masked.fa
```

## Get allele count
```
module load bioinfo-tools vcftools BEDTools
echo "Autosome"
vcftools --gzvcf ../data/vcf/${species}.A.filtered.vcf.gz --counts --out ../data/allele_count/${species}.A
echo "PAR"
vcftools --gzvcf ../data/vcf/${species}.PAR.filtered.vcf.gz --counts --out ../data/allele_count/${species}.PAR
echo "nonPAR"
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --out ../data/allele_count/${species}.nonPAR
```

# male nonPAR
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --keep ../data/samples/black_male.txt --window-pi 100000 --out ../data/allele_count/${species}.nonPAR.male

# Filter allele counts for repeats
awk '{print $1"\t"$4-1"\t"$5}' ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.out.gff | grep -v "#" > ../data/bed/ostrich_repeats.bed
echo "Autosome"
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.A.frq.count > ../data/allele_count/${species}.A.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.A.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.A.repeat.frq.count
echo "PAR"
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.PAR.frq.count > ../data/allele_count/${species}.PAR.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.PAR.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.PAR.repeat.frq.count
echo "nonPAR"
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.frq.count > ../data/allele_count/${species}.nonPAR.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.repeat.frq.count


# SNP count after repeat masking
# black
# Autosome 6222012
# PAR 289400
# nonPAR 55736


# Check heterozygote SNPs in the nonPAR in females. In the nonPAR since we should have only 1 Z in the females, we don't
# expect any female individual with genotype 0/1. Some could be gametologous genes. It is, however, important to identify
# such SNPs.
python check_nonPAR_genotype.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz > ../data/bed/${species}.nonPAR.female_het.bed
python vcf_to_genotypes.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz ../data/bed/${species}.nonPAR.genotypes.bed
python male_female_het_proportions.py > male_female_het_residuals.txt

# Figure of male and female heterozygosity

vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_male.txt --out ../data/allele_count/${species}.nonPAR.filtered.male
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_female.txt --out ../data/allele_count/${species}.nonPAR.filtered.female
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.male.frq.count > ../data/allele_count/${species}.nonPAR.filtered.male.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.male.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.female.frq.count > ../data/allele_count/${species}.nonPAR.filtered.female.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.female.frq.count.bed -b ../data/bed/ostrich_repeats.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count > ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.male.repeat.frq.count.bed -b ../data/bed/${species}.nonPAR.female_het.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""femalehet"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.male.repeat.femalehet.frq.count
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count > ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count.bed
bedtools intersect -a ../data/allele_count/${species}.nonPAR.filtered.female.repeat.frq.count.bed -b ../data/bed/${species}.nonPAR.female_het.bed -wao | \
awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""femalehet"}' | grep "PASS" | \
awk 'BEGIN{print "CHROM""\t""POS""\t""N_ALLELES""\t""N_CHR""\t""{ALLELE:COUNT}"}{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
> ../data/allele_count/${species}.nonPAR.filtered.female.repeat.femalehet.frq.count


# nonPAR allele count
python nonPAR_allele_count.py ../data/allele_count/${species}.nonPAR.filtered.female.repeat.femalehet.frq.count \
../data/allele_count/${species}.nonPAR.filtered.male.repeat.femalehet.frq.count > ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count



#*********************************************************************************************
#*********************************************************************************************
#*********************************************************************************************
#*********************************************************************************************
# Find overlap between windows and each of the functional categories
sbatch window_overlap.sh
# Make directory for genomic features per window
# mkdir -p ../result/genomic_features
# awk '{print $8}' ../result/genomic_features/black.par_scaf.100Kb.CDS.overlap.density.sorted.txt | \
# paste ../result/genomic_features/black.par_scaf.100Kb.AllRepeat.overlap.density.sorted.txt - | \
# awk 'BEGIN{print "Scaffold"" ""Start"" ""End"" ""Base_Count"" ""GC_Count"" ""Repeat_start"" ""Repeat_end"" ""Repeat_base"" ""Repeat_GC"" ""CDS_count"} \
# {print $0}' | less -S
# Calculate pi, theta, Tajima's D across Z and macrochromosomes
export result=/proj/snic2020-16-269/private/homap/ostrich_z/result/sfs_measures

# Autosomes
python SFS_measures_autosome.py ../data/allele_count/black.A.repeat.frq.count 20 ../data/bed/black.autosome.100Kb.intergenic.overlap.density.sorted.txt \
all > ${result}/black.autosome.sfs.txt
# PAR and nonPAR
# header = ["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
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

awk 'NR>1' ${result}/black.nonPAR.sfs.txt | cat ${result}/black.PAR.sfs.txt - > ${result}/black.Z.sfs.txt

python scaffold_to_chr_vcf.py ${result}/black.Z.sfs.txt > ${result}/black.Z.coordinates.sfs.txt


# Superscaffold36
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


# Annotating SNPs for intergenic and intronic regions and measure SFS measures for each functional category
bash annotate_snps.sh

# Male-Female FST
vcftools --gzvcf ../data/vcf/black.PAR.filtered.vcf.gz --chr 'superscaffold36' --recode --stdout | gzip -c > ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz
vcftools --gzvcf ../data/vcf/black.PAR.superscaffold36.filtered.vcf.gz --weir-fst-pop ../data/samples/black_male.txt --weir-fst-pop ../data/samples/black_female.txt \
--fst-window-size 1000 --out black_male_female_1Kb

# Male-Female FST
# Exclude female positions with heterozygous SNPs
awk '{print $1"\t"$3}' ../data/bed/${species}.nonPAR.female_het.bed > ../data/bed/${species}.nonPAR.female_het.txt
vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.pos.vcf.gz --exclude-positions ../data/bed/${species}.nonPAR.female_het.txt --stdout | gzip -c > ../data/vcf/LD_vcf/black.Z.pos.nohetfemalenonPAR.vcf.gz

vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.pos.nohetfemalenonPAR.vcf.gz --weir-fst-pop ../data/samples/black_male.txt --weir-fst-pop ../data/samples/black_female.txt --fst-window-size 100000 --out ../data/FST/black_male_female_Z_100Kb
python scaffold_to_chr_vcf.py ../data/FST/black_male_female_Z_100Kb.windowed.weir.fst > ../data/FST/black_male_female_Z_100Kb.windowed.weir.Z.coord.fst

# Get the distribution of SFS for autosomes, PAR and nonPAR
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
#*********************************************************************************************#
#*********************************************************************************************#
# Calculate Linkage Disequilibrium
#*********************************************************************************************#
#*********************************************************************************************#
# In the following, I use PopLDdecay, following the analysis done by
# Takeshi Kawakami to calculate LD on ostrich Z chromosome
# Download PopLDdecay into /proj/snic2020-16-269/private/homap/ostrich_z/bin
git clone https://github.com/BGI-shenzhen/PopLDdecay.git
cd PopLDdecay; chmod 755 configure; ./configure;
make;
mv PopLDdecay  bin/;

# The download and compiling went very smoothly.

# I use the r-squared measure of LD.

echo "PAR"
vcftools --gzvcf ../data/vcf/${species}.PAR.filtered.vcf.gz --hwe 0.005 --recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.PAR.hwe.filtered.vcf.gz
echo "nonPAR"
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --hwe 0.005 --recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.nonPAR.hwe.filtered.vcf.gz
echo "Autosome"
vcftools --gzvcf ../data/vcf/${species}.A.filtered.vcf.gz --hwe 0.005 --recode --recode-INFO-all --stdout | gzip -c > ../data/vcf/${species}.A.hwe.filtered.vcf.gz 

# Move the final autosomal vcf into the directory called LD_vcf
mkdir -p ../data/vcf/LD_vcf
mv ../data/vcf/${species}.*.hwe.filtered.vcf.gz ../data/vcf/LD_vcf

# Concatenate PAR and non-PAR positions
zgrep -v "#" ../data/vcf/LD_vcf/${species}.A.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.A.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.PAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.PAR.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.nonPAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.nonPAR.pos

cat ../data/vcf/LD_vcf/${species}.PAR.pos ../data/vcf/LD_vcf/${species}.nonPAR.pos | \
sort -k1,1 -k2n > ../data/vcf/LD_vcf/${species}.Z.pos

awk '{print $1"\t"$2-1"\t"$2}' ../data/vcf/LD_vcf/${species}.Z.pos > ../data/vcf/LD_vcf/${species}.Z.pos.bed 
python scaffold_to_chr_vcf.py ../data/vcf/LD_vcf/${species}.Z.pos.bed > ../data/vcf/LD_vcf/${species}.Z.pos.coordinates.txt

# Analysis of LD
# If required to do LD analysis in autosomes, use this bed file: grep -v '^C' ../data/bed/autosomes.bed > ../data/bed/autosomes.minusC.bed
sbatch ld_run.sh

# Calculate mean LD in xx kb windows by sliding window analysis
for scaffold in $(cat ../data/bed/par_scaf.bed ../data/bed/nonpar_scaf.bed \
| grep -v "^chrom" | cut -f1 | sort | uniq)
do
echo $scaffold
zcat ../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.gz | \
grep -v "^#" | awk '($9>=500 && $9<=50000)' > \
../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.05-500
Rscript LD_slidingWin_v1.R ../data/LD/scaf.split/${species}.${scaffold}.pairwise.LD.05-500 \
../data/bed/z_scaf.bed 200000 50000
done

# Get a figure of the LD for 3 subspecies for superscaffold36. The graph is
# generated in the bin directory.
Rscript LD_three_species.R \
../data/LD/scaf.split/black.superscaffold36.pairwise.LD.05-500 \
../data/LD/scaf.split/blue.superscaffold36.pairwise.LD.05-500 \
../data/LD/scaf.split/red.superscaffold36.pairwise.LD.05-500 \
../data/bed/z_scaf.bed 200000 50000

grep -v "^#" | awk '($9>=1000 && $9<=200000)' >
# Next, we are interested to get a more clear picture of LD across the
# PAR boundary in each sex. However, if we do this analysis per
# subspecies, we will have very few samples. Since Black and Blue
# subspecies are closely related and have a similar population history
# see (PSMC graphs), we pool black and blue together for each species
# and for the shared SNPs between the two subspecies.
module load bioinfo-tools vcftools
mkdir -p ../data/LD/sex_specific
# Select SNPs that are shared among black and blue
cat ../data/vcf/LD_vcf/black.Z.pos | grep "^superscaffold36" | \
awk '($2-3524263>-100000 && $2-3524263<100000)' > ../data/LD/sex_specific/commonSNP.B.boundary200kb

cat ../data/vcf/LD_vcf/black.Z.pos | grep "^superscaffold36" | \
awk '($2-3524263>-50000 && $2-3524263<50000)' > ../data/LD/sex_specific/commonSNP.B.boundary100kb


# convert vcf to plink near PAR boundary (200 kb centered at the boundary)
chr=superscaffold36
vcftools --gzvcf ../data/vcf/LD_vcf/black.Z.superscaffold36.vcf.gz --positions ../data/LD/sex_specific/commonSNP.B.boundary100kb \
--plink --out ../data/LD/sex_specific/black.superscaffold36.bothsexes.100kb.boundary

plink --file ../data/LD/sex_specific/black.superscaffold36.bothsexes.100kb.boundary --r2 square --out ../data/LD/sex_specific/black.superscaffold36.bothsexes.100Kb.boundary

awk '{print $1"\t"$2-1"\t"$2}' ../data/LD/sex_specific/commonSNP.B.boundary100kb | cat ../data/bed/par_scaf.bed - > test 

python scaffold_to_chr_vcf.py test > ../data/bed/boundary.z.coordinates.txt
# Join ped files per sex
cat ../data/LD/sex_specific/*male.ped > \
../data/LD/sex_specific/BB.superscaffold36.boundary.male.ped.join
cat ../data/LD/sex_specific/*female.ped > \
../data/LD/sex_specific/BB.superscaffold36.boundary.female.ped.join

# modify plink infile for Haploview. Note that Haploview takes plink format. However
# .ped file needs a header. Instead of using plink format, I use Linkage format that
# needs .ped file as is and .info file that requires two columns indicating SNP name
# and position
cat ../data/LD/sex_specific/black.superscaffold36.boundary.female.map | cut -f2,4 > \
../data/LD/sex_specific/BB.superscaffold36.boundary.map

# Visualize by Haploview gui in Mac
# Open and upload the ped and map file for male and female in the Haploview
# and save the figure as svg and compressed png file.
java -jar Haploview.jar # TO DO


#*********************************************************************************************#
# Copy results to results directory from LD directory in data
cp *.LDdecay.bin.gz ../../../result/LD
cp *pairwise.LD.gz ../../../result/LD
#*********************************************************************************************#
# Many analyses such as recombination rate, genetic diversity, Fst, etc, are done in windows.
# Here, we generate bed files with windows across Z scaffolds. We will then overlap these
# windows with our feature of interest with bedtools.
#*********************************************************************************************#
# Get windows per scaffolds
#*********************************************************************************************#
mkdir -p window
echo "Create windows"
for window in 200Kb 500Kb 1000Kb
do
echo $window
w_size=$(echo $window | sed 's/Kb/000/g')
python sliding_window_ZA.py ../data/bed/z_scaf.bed ${w_size} > \
../data/linkage_map/ostrich.Z.${window}.bed
python sliding_window_ZA.py ../data/bed/z_scaf.bed ${w_size} > \
../data/window/ostrich.Z.${window}.bed
done
#*********************************************************************************************#
# Recombination rate from the genetic map
# Obtain per window sex averaged, male and female recombination rates
#*********************************************************************************************#
# Load the necessary modules
module load bioinfo-tools BEDTools
#_____________________________________________________________________________________________#
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

#*********************************************************************************************#
# PSMC
#*********************************************************************************************#
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

# The blue subspecies
for i in {1..100}
do
mkdir -p ../data/psmc/blue.CI.${i}.dir
sbatch psmc_CI.sh ../data/psmc/blue.CI.${i}.dir ../data/psmc/P1878_121.autosome.split.psmcfa
done

mkdir -p blue_resampling
j=0
for i in blue.CI.*
do
echo ${i}/round_repeat-1.psmc
((j=j+1))
echo $j
cp ${i}/round_repeat-1.psmc blue_resampling
mv blue_resampling/round_repeat-1.psmc blue_resampling/round_repeat-${j}.psmc
done
cp blue.CI.1.dir/round_repeat-1.psmc blue_resampling


for i in blue.CI.*
do
echo ${i}/round_repeat-1.psmc
cat ${i}/round_repeat-1.psmc >> blue.CI.resampling.txt
done

cat P1878_121.autosome.psmc.out blue.CI.resampling.txt > temp && mv temp blue.CI.resampling.txt

# The red subspecies
for i in {1..100}
do
mkdir -p ../data/psmc/red.CI.${i}.dir
sbatch psmc_CI.sh ../data/psmc/red.CI.${i}.dir ../data/psmc/P1878_128.autosome.split.psmcfa
done

for i in red.CI.*
do
echo ${i}/round_repeat-1.psmc
cat ${i}/round_repeat-1.psmc >> red.CI.resampling.txt
done

cat P1878_128.autosome.psmc.out red.CI.resampling.txt > temp && mv temp red.CI.resampling.txt

#*********************************************************************************************#
# Estimating population scaled recombination rate (rho)
#*********************************************************************************************#
# Run LDhat with genotype data
# Producing input for ldhat run from bin
# In /Users/homapapoli/Documents/projects/ostrich_Z/ldhat_dir
# conda activate 
# conda activate ldhat
./vcf_to_ldhat_out.py black.PAR.hwe.filtered.vcf.gz par_scaf.bed 20 19 superscaffold36 100
export vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/black.PAR.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz
./vcf_to_LDhat_genotype.py $vcf par_scaf.bed 20 19 superscaffold36 100
#*********************************************************************************************#
# In personal computer /Users/homapapoli/Documents/projects/ostrich_Z/4_Ostrich_polymorphism/manuscript_tables
# Create Z coordinates 
for species in black blue red
do 
echo $species
python scaffold_to_chr_vcf.py sfs_measures/${species}.PAR.sfs.txt sfs_measures/${species}.nonPAR.sfs.txt > sfs_measures/${species}.sfs.Z.txt
done

# Prepare LD output for plotting in R
# In /proj/snic2020-16-269/private/homap/ostrich_z/data/LD/scaf.split 
mkdir -p LD_chromosome_plot
# Run LD_plot_scaffold.sh in the /proj/snic2020-16-269/private/homap/ostrich_z/bin

# Z coordinates for GFF

#*********************************************************************************************#
# Get the PAR fasta sequence
reference/black.repeat.depth.masked.fa
# Convert fasta sequence into consesuns for each individual with SNPs
../data/vcf/${species}.PAR.filtered.vcf.gz


for segment in superscaffold36:3524263-9394175 superscaffold35:1-4625539 superscaffold54:1-16379243 superscaffold26:1-25310599
do 
echo $segment
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_107 -o consensus_fasta/P1878_107.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_108 -o consensus_fasta/P1878_108.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_109 -o consensus_fasta/P1878_109.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_110 -o consensus_fasta/P1878_110.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_111 -o consensus_fasta/P1878_111.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_112 -o consensus_fasta/P1878_112.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_113 -o consensus_fasta/P1878_113.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_114 -o consensus_fasta/P1878_114.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_115 -o consensus_fasta/P1878_115.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_116 -o consensus_fasta/P1878_116.${segment}.fa
done

# Produce 2000 SNP intervals for each scaffold
# Running in /proj/snic2020-16-269/private/homap/ostrich_z/result/ldhat/ldhat_dir/
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold36 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold35 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold54 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold26 1000

# Create the exact likelihood table
./interval -seq superscaffold36.2000.500.1.sites.txt -loc superscaffold36.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold35.2000.500.1.sites.txt -loc superscaffold35.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold35.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold54.2000.500.1.sites.txt -loc superscaffold54.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold54.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold26.2000.500.1.sites.txt -loc superscaffold26.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold26.first.run -its 10000000 -bpen 5 -samp 2000

for i in {1..24} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold36.2000.500.${i}.sites.txt superscaffold36.2000.500.${i}.locs.txt superscaffold36.first.runnew_lk.txt superscaffold36.2000.500.${i}.
done
./interval -seq test/superscaffold36.2000.500.25.sites.txt -loc test/superscaffold36.2000.500.25.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold36.2000.500.25. -its 10000000 -bpen 5 -samp 2000

for i in {1..18} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold35.2000.500.${i}.sites.txt superscaffold35.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold35.2000.500.${i}.
done
./interval -seq test/superscaffold35.2000.500.19.sites.txt -loc test/superscaffold35.2000.500.19.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold35.2000.500.19. -its 10000000 -bpen 5 -samp 2000

for i in {7..58} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold54.2000.500.${i}.sites.txt superscaffold54.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold54.2000.500.${i}.
done
./interval -seq test/superscaffold54.2000.500.59.sites.txt -loc test/superscaffold54.2000.500.59.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold54.2000.500.59. -its 10000000 -bpen 5 -samp 2000


for i in {1..100} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold26.2000.500.${i}.sites.txt superscaffold26.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.2000.500.${i}.
done
./interval -seq test/superscaffold26.2000.500.101.sites.txt -loc test/superscaffold26.2000.500.101.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold26.2000.500.101. -its 10000000 -bpen 5 -samp 2000


for i in {1..24} 
do
echo $i
./interval -seq superscaffold36.2000.500.${i}.sites.txt -loc superscaffold36.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {1..18} 
do
echo $i
./interval -seq superscaffold35.2000.500.${i}.sites.txt -loc superscaffold35.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold35.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {7..58} 
do
echo $i
./interval -seq superscaffold54.2000.500.${i}.sites.txt -loc superscaffold54.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold54.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {1..100} 
do
echo $i
./interval -seq superscaffold26.2000.500.${i}.sites.txt -loc superscaffold26.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold26.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 20 5 superscaffold36 10


for i in {1..10} 
do
echo $i
./interval -seq superscaffold36.20.5.${i}.sites.txt -loc superscaffold36.20.5.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.20.5.${i}. -its 10000000 -bpen 5 -samp 2000
done

# Plot the gene models along the Z chromosome

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

# Chromosome 4 and 5
NC_006091.5     chromosome_4
NC_006092.5     chromosome_5


git remote add origin https://github.com/Homap/ostrich_Z

…or create a new repository on the command line
echo "# ostrich_PAR_analysis" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/Homap/ostrich_PAR_analysis.git
git push -u origin main
…or push an existing repository from the command line
git remote add origin https://github.com/Homap/ostrich_PAR_analysis.git
git branch -M main
git push -u origin main

# Creating SSH key
# When creating a new github directory in Uppmax, you need to add the ssh key to the directory or you cannot push to it. Find the ssh key under ~/.ssh
git add file
git commit -m "something has changed here"
git push -u origin main

git rm <file>

git commit -m "Deleted the file from the git repository"
git push -u origin main


