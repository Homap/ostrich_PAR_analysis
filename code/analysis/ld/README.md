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

Decay of LD is affected by recombination rate and the number of generations of recombination. Therefore, 
investigating LD decay may reveal the population recombination history. 

LD decay is computed by running the following script. LD decay is calculated for scaffolds residing on chromosome
4 and 5, PAR and nonPAR scaffolds. Additionally, LD decay is calculated for a 500 Kb region at the end, 500 Kb in 
the middle and 500 Kb closest to the PAR boundary.  

`bash 2_ld_decay_run.sh`

Table 1. Output of LD decay - We output only r<sup>2</sup>, therefore, D' is NA.
|#Dist  | Mean_r<sup>2</sup>|  Mean_D'| Sum_r<sup>2</sup>|Sum_D' | NumberPairs |
| ----- | ------- | ------- | ------ | ----- | ----------- |
|1      | 0.4240  |NA    |  17.8087 |NA  |    42 |
|2      | 0.4062 | NA    |  6.9055 | NA    |  17 |
|3      | 0.4209 | NA     | 4.2091 | NA    |  10 |

Table 2. Output of binned LD decay used for plotting.
|#Dist  | Mean_r<sup>2</sup>       | Mean_D' |Sum_r<sup>2</sup>| Sum_D' |  NumberPairs|
| ----- | -------------- | ------- | ----- | ------ | ----------- |
|10     | 0.397282681564246  |     NA  |    71.1136 |NA    |  179|
|20     | 0.31554364640884   |     NA  |    57.1134| NA     | 181|
|30     | 0.262394318181818  |     NA  |   46.1814 |NA    |  176|

Pairwise LD is calculated for all categories as above in addition for a 100 Kb region spanning the PAR-nonPAR boundary on superscaffold36.

`bash 3_ld_pairwise_run.sh`

Table 3. Output of pairwise LD calculation.
|#chr  |  Site1  | Site2 |  D'    |  LOD    | r<sup>2</sup>    | CIlow  | CIhi    |Dist|
| ---- | ------- | ----- | ------ | ------- | ------ | ------ | ------- | --- |
|superscaffold26 |1306  |  1622   | 0.2308 | 0.0458 | 0.0110 | 0.02  |  0.83   | 316|
|superscaffold26 |1306  |  1641   | 0.2857 | 0.0850 | 0.0212 | 0.03  |  0.83   | 335|
|superscaffold26 |1306  |  1766   | 1.0000 | 1.0049 | 0.2063 | 0.16   | 0.99   | 460|

## Final LD output for statistical analysis

### LD decay

- Autosomes 4 and 5
    ```
    export ld_dir=../../../data/ld/ld_decay 
    gunzip ${ld_dir}/autosome/chr4/*.stat.gz
    gunzip ${ld_dir}/autosome/chr5/*.stat.gz

    cat ${ld_dir}/autosome/chr4/*.stat > ${ld_dir}/autosome/chr4_pairwise_stat
    cat ${ld_dir}/autosome/chr5/*.stat > ${ld_dir}/autosome/chr5_pairwise_stat

    perl PopLDdecay/bin/Plot_OnePop.pl -inFile ${ld_dir}/autosome/chr4_pairwise_stat -output ${ld_dir}/autosome/chr4
    perl PopLDdecay/bin/Plot_OnePop.pl -inFile ${ld_dir}/autosome/chr5_pairwise_stat -output ${ld_dir}/autosome/chr5

    cat ${ld_dir}/autosome/chr4_pairwise_stat ${ld_dir}/autosome/chr5_pairwise_stat > ${ld_dir}/autosome/chr4_chr5_pairwise_stat
    perl PopLDdecay/bin/Plot_OnePop.pl -inFile ${ld_dir}/autosome/chr4_chr5_pairwise_stat -output ${ld_dir}/autosome/chr4_5

    rm ${ld_dir}/autosome/*png ${ld_dir}/autosome/*pdf ${ld_dir}/autosome/*stat
    ```

- Z

    `export ld_dir=../../../data/ld/ld_decay`
    `gunzip ${ld_dir}/z/*.stat.gz`

    - whole PAR
        ```
        cat ${ld_dir}/z/superscaffold26.LDdecay.stat \
        ${ld_dir}/z/superscaffold54.LDdecay.stat \
        ${ld_dir}/z/superscaffold35.LDdecay.stat \
        ${ld_dir}/z/superscaffold36.par.LDdecay.stat > ${ld_dir}/z/PAR_Lddecay_stat
        perl PopLDdecay/bin/Plot_OnePop.pl -inFile ${ld_dir}/z/PAR_Lddecay_stat \
        -output ${ld_dir}/z/par/PAR
        rm -f ${ld_dir}/z/par/*p*
        ```
    - start PAR: 500Kb farthest from the SDR
        ```
        cp ${ld_dir}/z/supersaffold26.500Kb.LDdecay.bin.gz ${ld_dir}/z/par
        mv ${ld_dir}/z/par/supersaffold26.500Kb.LDdecay.bin.gz ${ld_dir}/z/par/start.par.500Kb.bin.gz
        ```
    - mid PAR: 500Kb in the middle of the PAR
        ```
        cp ${ld_dir}/z/supersaffold54.500Kb.LDdecay.bin.gz ${ld_dir}/z/par
        mv ${ld_dir}/z/par/supersaffold54.500Kb.LDdecay.bin.gz ${ld_dir}/z/par/mid.par.500Kb.bin.gz
        ```
    - end PAR: 500 Kb closest to the SDR
        ```
        cp ${ld_dir}/z/supersaffold36.500Kb.LDdecay.bin.gz ${ld_dir}/z/par
        mv ${ld_dir}/z/par/supersaffold36.500Kb.LDdecay.bin.gz ${ld_dir}/z/par/end.par.500Kb.bin.gz
        ```
    - nonPAR
        ```
        cat ${ld_dir}/z/superscaffold62.LDdecay.stat \
        ${ld_dir}/z/superscaffold63.LDdecay.stat \
        ${ld_dir}/z/superscaffold67.LDdecay.stat \
        ${ld_dir}/z/superscaffold69-1.LDdecay.stat \
        ${ld_dir}/z/superscaffold83.LDdecay.stat \
        ${ld_dir}/z/superscaffold88.LDdecay.stat \
        ${ld_dir}/z/superscaffold92.LDdecay.stat \
        ${ld_dir}/z/superscaffold93.LDdecay.stat \
        ${ld_dir}/z/superscaffold36.nonpar.LDdecay.stat > ${ld_dir}/z/nonPAR_Lddecay_stat
        perl PopLDdecay/bin/Plot_OnePop.pl -inFile ${ld_dir}/z/nonPAR_Lddecay_stat \
        -output ${ld_dir}/z/nonpar/nonPAR
        rm -f ${ld_dir}/z/nonpar/*p*
        ```

## Calculate mean LD in 200 kb windows by sliding window analysis

For this, we are looking at each scaffold separately and we are missing the connection between the scaffolds. 
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













