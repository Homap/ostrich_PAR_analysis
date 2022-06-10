# Population scaled recombination rate (ρ)

We estimate the population scaled recombination rate using the method LDhat (Auton and McVean Genome Research, 2007). We specifically use the *interval* program in the LDhat package which
"uses a Bayesian reversible-jump Markov chain Monte Carlo (rjMCMC) scheme to fit a piecewise-constant model of recombination rate variation. However,
rather than calculating the full coalescent likelihood, a composite-likelihood is employed (Hudson 2001)." We lacked phased data, we therefore run *interval* on genotypes. We keep each 
file to 2000 SNPs since it is recommended by the manual of LDhat.

## Generate input files for *interval* of LDhat
The *interval* program in LDhat requires three types of input. 
    - Sites: 
      ```     
      10 2000 2
      >P1878_107
      00022000022022
      >P1878_108
      00020200220022
      ```
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
./interval -seq $sites -loc $locs -lk $lk -prefix $out_prefix -its 10000000 -bpen 5 -samp 2000

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