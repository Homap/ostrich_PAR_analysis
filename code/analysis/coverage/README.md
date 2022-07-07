# Coverage of BAM files and at the PAR-nonPAR boundary

- Coverage for the bam files <br>
`sbatch get_depth.sh`

- Separate the coverage output to PAR, nonPAR and autosome <br>
```
cut -f1 ../../../data/bed/par_scaf.bed | grep -f - ../../../data/coverage/coverage_per_site.txt > ../../../data/coverage/coverage_per_site.par.txt
cut -f1 ../../../data/bed/nonpar_scaf.bed | grep -f - ../../../data/coverage/coverage_per_site.txt > ../../../data/coverage/coverage_per_site.nonpar.txt
cut -f1 ../../../data/lastz/gg_chr4_ostrich.bed | grep -f - ../../../data/coverage/coverage_per_site.txt > ../../../data/coverage/coverage_per_site.chr4.txt
cut -f1 ../../../data/lastz/gg_chr5_ostrich.bed | grep -f - ../../../data/coverage/coverage_per_site.txt > ../../../data/coverage/coverage_per_site.chr5.txt
cat ../../../data/coverage/coverage_per_site.chr4.txt ../../../data/coverage/coverage_per_site.chr5.txt > ../../../data/coverage/coverage_per_site_chr4_5.txt
```

- Obtain the median coverage for each individual and category <br>
```
python median_coverage.py ../../../data/coverage/coverage_per_site.nonpar.txt 10 ../../../data/samples/black_samples.txt > ../../../data/coverage/coverage_median.nonpar.txt
python median_coverage.py ../../../data/coverage/coverage_per_site.par.txt 10 ../../../data/samples/black_samples.txt > ../../../data/coverage/coverage_median.par.txt
python median_coverage.py ../../../data/coverage/coverage_per_site_chr4_5.txt 10 ../../../data/samples/black_samples.txt > ../../../data/coverage/coverage_median.chr4_5.txt
```

`rm -i ../../../data/coverage/coverage_per_site.chr4.txt ../../../data/coverage/coverage_per_site.chr5.txt`

## Coverage at the PAR-nonPAR boundary on superscaffold36

```
module load bioinfo-tools samtools/1.14
bamaddr=../../../data/bam 
samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_107.bam ${bamaddr}/P1878_108.bam ${bamaddr}/P1878_109.bam ${bamaddr}/P1878_110.bam ${bamaddr}/P1878_111.bam \
> ../../../data/coverage/male_boundary_coverage.txt

samtools depth -r superscaffold36:3510000-3530000 ${bamaddr}/P1878_112.bam ${bamaddr}/P1878_113.bam ${bamaddr}/P1878_114.bam ${bamaddr}/P1878_115.bam ${bamaddr}/P1878_116.bam \
> ../../../data/coverage/female_boundary_coverage.txt
```

