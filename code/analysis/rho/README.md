# Population scaled recombination rate (ρ)

We estimate the population scaled recombination rate using the method LDhat (Auton and McVean Genome Research, 2007). We specifically use the *interval* program in the LDhat package which
"uses a Bayesian reversible-jump Markov chain Monte Carlo (rjMCMC) scheme to fit a piecewise-constant model of recombination rate variation. However,
rather than calculating the full coalescent likelihood, a composite-likelihood is employed (Hudson 2001)." We lacked phased data, we therefore run *interval* on genotypes. We keep each 
file to 2000 SNPs since it is recommended by the manual of LDhat.

## Generate input files for *interval* of LDhat
The *interval* program in LDhat requires three types of input. 

- Sites

      10 2000 2
      >P1878_107
      00022000022022
      >P1878_108
      00020200220022

- Locus

    ```
    2000 574366 L
    1
    70
    102
    798
    ```

- Likelihood Lookup Table

    ```
    20 3194
    1 0.00100
    100 101.000000
    276 #   0   0   0   0   0   0   0   0   0   2   2   0   0   2   0   4  :   -37.60  -37.50
    277 #   0   0   0   0   0   0   0   0   0   2   2   0   0   2   0   4  :   -37.60  -37.50
    ```

A likelihood lookup table is also required to run *interval*. We use the program *lkgen* to generate lookup tables from the precomputed
likelihood table available by the LDhat program since minor differences in theta do not appear to strongly influence 
the results. The average pairwise nucleotide diversity for the PAR and autosome in ostrich is 0.002 and for SDR is 0.0007.
We therefore use the likelihood table for 50 sequences and theta = 0.001 to obtain a likelihood table for 20 sequences
and use that in the *interval* program.

## Likelihood look up table
There are a number of pre-computed likelihood tables provided with the LDhat package. The one we use is lk_n50_t0.001.gz which 
is the likelihood for 50 chromosomes and a theta of 0.001. We use the *lkgen* to generate the lookup table. 

```
mkdir -p ../../../data/rho/ldhat_input/lk_LUT
export lk_LUT_path=../../../data/rho/ldhat_input/lk_LUT
```

- For autosomes and PAR: 10 individuals (20 chromosomes)

`./LDhat/lkgen -lk LDhat/lk_files/lk_n50_t0.001 -nseq 20 -prefix ${lk_LUT_path}/auto_PAR`

- For the nonPAR: 5 males (10 chromosomes)

`./LDhat/lkgen -lk LDhat/lk_files/lk_n50_t0.001 -nseq 10 -prefix ${lk_LUT_path}/nonPAR`


./interval -seq $sites -loc $locs -lk $lk -prefix $out_prefix -its 10000000 -bpen 5 -samp 2000

For the SDR, calculate Rho only in males and then for the sex-averaged recombination rate, do 2/3*(male recombination
rate).

`bash run_vcf_to_ldhat_input.sh`

The python script, `vcf_to_ldhat_out.py` produces a file ending with `scaffold.pos.txt`that stores
the original positions of sites on the scaffold, later used for plotting:

```
2000 574366 L
161
230
262
958
```

When running interval, check the likelihood curve to see if the maximum is reached. Remember if the region is too large, there might be
no LD and therefore, rho computation is wrong.


python vcf_to_ldhat_out.py ../../../data/vcf/z_vcf/z_vcf.gz ../../../data/bed/z_scaf.bed 2000 500 superscaffold26 2 ./
python vcf_to_ldhat_out.py ../../../data/vcf/z_vcf/z_vcf.gz ../../../data/bed/z_scaf.bed 1000 250 superscaffold26 2 ./
python vcf_to_ldhat_out.py ../../../data/vcf/z_vcf/z_vcf.gz ../../../data/bed/z_scaf.bed 500 100 superscaffold26 2 ./

./LDhat/interval -seq superscaffold26.2000.500.1.sites.txt -loc superscaffold26.2000.500.1.locs.txt -lk ostrich_complete_theta0.002.txt -prefix superscaffold26.2000.500.1.complete0.002 -its 10000000 -bpen 5 -samp 2000

for i in {1..2}
do
sbatch ldhat_interval.sh superscaffold26.2000.500.${i}.sites.txt superscaffold26.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.2000.${i}.
sbatch ldhat_interval.sh superscaffold26.1000.250.${i}.sites.txtsuperscaffold26.1000.250.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.1000.250.${i}.
sbatch ldhat_interval.sh superscaffold26.500.100.${i}.sites.txt superscaffold26.500.100.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.500.100.${i}.
done

python vcf_to_ldhat_out.py ../../../data/vcf/z_vcf/z_vcf.gz ../../../data/bed/z_scaf.bed 100 10 superscaffold26 2 ./
sbatch ldhat_interval.sh superscaffold26.100.10.${i}.sites.txt superscaffold26.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.2000.${i}.

grep 'LK' slurm-27664288.out | awk '{print $5}' > lk2000
grep 'LK' slurm-27664290.out | awk '{print $5}' > lk500
grep 'LK' lk100.10.bpen5 | awk '{print $5}' > lk100
grep 'LK' lk100.10.bpen20 | awk '{print $5}' > lk100.20

grep 'blocks' slurm-27664288.out | awk '{print $10}' > blocks2000
grep 'blocks' slurm-27664290.out | awk '{print $10}' > blocks500
grep 'blocks' lk100.10.bpen5 | awk '{print $10}' > blocks100
grep 'blocks' lk100.10.bpen20 | awk '{print $10}' > blocks100.20

grep 'Map' slurm-27664288.out | awk '{print $15}' > Map2000
grep 'Map' slurm-27664290.out | awk '{print $15}' > Map500
grep 'Map' lk100.10.bpen5 | awk '{print $15}' > map100
grep 'Map' lk100.10.bpen20 | awk '{print $15}' > map100.20

grep 'LK' example_lk | awk '{print $5}' > examplelk
grep 'blocks' example_lk | awk '{print $10}' > exampleblock
grep 'Map' example_lk | awk '{print $15}' > examplemap

After changing the distances to kb
grep 'LK' slurm-27675485.out | awk '{print $5}' > lk2000.kb
grep 'blocks' slurm-27675485.out | awk '{print $10}' > blocks2000.kb
grep 'Map' slurm-27675485.out | awk '{print $15}' > map2000.kb

With new likelihood table
grep 'LK' slurm-27681299.out | awk '{print $5}' > lk2000.newlk.25kb
grep 'blocks' slurm-27681299.out | awk '{print $10}' > block2000.newlk.25kb
grep 'Map' slurm-27681299.out | awk '{print $15}' > map2000.newlk.25kb

grep 'LK' slurm-27681300.out | awk '{print $5}' > lk2000.25kb
grep 'blocks' slurm-27681300.out | awk '{print $10}' > block2000.25kb
grep 'Map' slurm-27681300.out | awk '{print $15}' > map2000.25kb

Double check rho_Ne_reveresed.txt

# Run it for more chains

for i in {0..50}
do
echo $i
sbatch -J bpen.${i} ldhat_interval.sh superscaffold26.2000.500.1.sites.txt superscaffold26.2000.500.1.locs.txt ostrich_genotypenew_lk.txt superscaffold26.2000.500.1.complete.bp.${i}.25mil $i >> slurm.submission
done


# In /Users/homapapoli/Documents/projects/ostrich_Z/ldhat_dir
# conda activate 
# conda activate ldhat
./vcf_to_ldhat_out.py black.PAR.hwe.filtered.vcf.gz par_scaf.bed 20 19 superscaffold36 100
export vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/black.PAR.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz
./vcf_to_LDhat_genotype.py $vcf par_scaf.bed 20 19 superscaffold36 100
#*********************************************************************************************#
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

# Just add the non-PAR segments
The nonPAR needs to be calculated only using the male SNPs