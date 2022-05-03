# Steps for filtering the VCF files and obtaining the input for genetic diversity, LD and rho analyses

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

## Filter allele counts for repeats
```
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
```

| Category | Number of SNPs after repeat masking |
| :--------- | :-----: |
| Autosome   | 6222012 |
| PAR  | 289400 |
| nonPAR | 55736 |


Check heterozygote SNPs in the nonPAR in females. In the nonPAR since we should have only 1 Z in the females, we don't
expect any female individual with genotype 0/1. Some could be gametologous genes. It is, however, important to identify
such SNPs.

```
python check_nonPAR_genotype.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz > ../data/bed/${species}.nonPAR.female_het.bed
python vcf_to_genotypes.py ../data/vcf/${species}.nonPAR.filtered.vcf.gz ../data/bed/${species}.nonPAR.genotypes.bed
```

## Two counts files for males and females in the nonPAR
```
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_male.txt --out ../data/allele_count/${species}.nonPAR.filtered.male
vcftools --gzvcf ../data/vcf/${species}.nonPAR.filtered.vcf.gz --counts --keep ../data/samples/${species}_female.txt --out ../data/allele_count/${species}.nonPAR.filtered.female
```

## Filter out heterozgous sites in females from both males and females
```
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
```

## Obtain nonPAR allele count
The nonPAR allele count should be adjusted for females having 5 chromosomes. In VCF, females are called as diploid
in the nonPAR so 0/0 and 1/1 should be changed into 0 and 1.

```
python nonPAR_allele_count.py ../data/allele_count/${species}.nonPAR.filtered.female.repeat.femalehet.frq.count \
../data/allele_count/${species}.nonPAR.filtered.male.repeat.femalehet.frq.count > ../data/allele_count/${species}.nonPAR.filtered.adjusted.frq.count
```

## Find overlap between windows and each of the functional categories
```
sbatch window_overlap.sh
```
