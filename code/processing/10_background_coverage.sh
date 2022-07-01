# Background Filtering for Coverage

# Coverage for the bam files
sbatch get_depth.sh

# Create a bed file with variants with coverage less than 5 or more than 70
for ind in P1878_107 P1878_108 P1878_109 P1878_110 P1878_111 P1878_112 P1878_113 P1878_114 P1878_115 P1878_116
do
echo $ind
# Add average to the python script
sbatch get_background_coverage.sh ../../data/coverage/${ind}.bed ../../data/coverage/${ind}.5.70.tofilter.bed
done

 


Use the repeat masked genome where repeats are already N. For analysis where you need the number 
of filtered invariant sites, you can simply count non-N sites in each window and that will be 
the number of bases in window.

`module load bioinfo-tools BEDTools/2.29.2`
```
bedtools maskfasta -fi ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked \
-bed ../data/vcf/${species}.all.ldepth.mean.filter.bed -fo ../data/reference/${species}.repeat.depth.masked.fa
```