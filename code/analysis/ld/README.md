# Linkage Disequilibrium (LD)

We use r<sup>2</sup> as the measure of LD. We are interested in two quantities. Pairwise LD and LD decay. 
Pairwise LD computes LD in pairs of sites in a given window size. LD decay gives us the average distance in which
LD reaches a value where non-random association among alleles is small (for example r<sup>2</sup>= 0.1).

We compute LD on autosome and Z chromosmoe. In ostrich genome assembly, only chromosome Z has been
developed into a chromosome-level assembly by use of linkage map. We lack a chromosome level assembly for the autosomes. 
In ostrich and most other birds with studied karyotype, chromosome Z has similar physical and genetic length to chromosomes
4 and 5. Using LASTZ, we obtain the match between ostrich scaffolds and chicken chromosomes 4 and 5. Chimeric scaffolds, scaffolds
containing sequences from more than one chromosome, might exist in assemblies. We therefore filter out chimeric scaffolds and keep scaffolds
with only one hit to either chromosome 4 or 5 of chicken.

The sequence in Z is composed of the PAR and nonPAR. We obtain separate measures of LD decay for:
- Whole PAR
    - End of PAR
    - Mid PAR
    - PAR boundary
- nonPAR

In the following, we use PopLDdecay, to calculate LD. Download and install PopLDdecay:

```
git clone https://github.com/BGI-shenzhen/PopLDdecay.git
cd PopLDdecay; chmod 755 configure; ./configure;
make;
mv PopLDdecay  bin/;
```

PopLDdecay works on VCF files, we can, therefore, split the VCF to scaffolds corresponding to regions
we are interested in. These regions include chromosome 4, 5 and Z chromosome. The Z chromosome contains
the PAR and nonPAR with superscaffold36 spanning the boundary. 

`sbatch 1_LD_vcf_split.sh`

LD decay is computed by running the following script. LD decay is calculated for scaffolds residing on chromosome
4 and 5, PAR and nonPAR scaffolds. Additionally, LD decay is calculated for a 500 Kb region at the end, 500 Kb in 
the middle and 500 Kb closest to the PAR boundary.  

`sbatch 2_ld_decay_run.sh`

Pairwise LD is calculated for all categories as above in addition for a 100 Kb region spanning the PAR-nonPAR boundary
on superscaffold36.

`sbatch 3_ld_pairwise_run.sh`


## Calculate mean LD in 200 kb windows by sliding window analysis
```
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
```










