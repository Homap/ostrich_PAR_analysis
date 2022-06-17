# Population scaled recombination rate (ρ)

We estimate the population scaled recombination rate using the method LDhat (Auton and McVean Genome Research, 2007). We specifically use the *interval* program in the LDhat package which
"uses a Bayesian reversible-jump Markov chain Monte Carlo (rjMCMC) scheme to fit a piecewise-constant model of recombination rate variation. However,
rather than calculating the full coalescent likelihood, a composite-likelihood is employed (Hudson 2001)." We lacked phased data, we therefore run *interval* on genotypes. We keep each 
file to 2000 SNPs since it is recommended by the manual of LDhat.

## Generate input files for *interval* program of LDhat
The *interval* program in LDhat requires three types of input. 

- Likelihood Lookup Table

    ```
    20 3194
    1 0.00100
    100 101.000000
    276 #   0   0   0   0   0   0   0   0   0   2   2   0   0   2   0   4  :   -37.60  -37.50
    277 #   0   0   0   0   0   0   0   0   0   2   2   0   0   2   0   4  :   -37.60  -37.50
    ```

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

    - `./LDhat/lkgen -lk LDhat/lk_files/lk_n50_t0.001 -nseq 20 -prefix ${lk_LUT_path}/auto_PAR`

- For the nonPAR: 5 males (10 chromosomes)

    - `./LDhat/lkgen -lk LDhat/lk_files/lk_n50_t0.001 -nseq 10 -prefix ${lk_LUT_path}/nonPAR`

## Sites and locus files 

We would like to calculate the rho for the whole Z chromosome however, if running LDhat on long 
scaffolds or chromosomes, the algorithm may not coverge. Instead, we split the data into chunks 
of 2000 SNPs with 200 SNPs overlap between windows. 

`bash 1_run_vcf_to_ldhat_input.sh`

The python script in the script above called `vcf_to_ldhat_out.py` produces an additional file ending with `pos.txt`that stores
the original positions of sites on the scaffold, later used for plotting:

- Pos

```
2000 574366 L
161
230
262
958
```

## Running *interval*

We run 3 chain, each of 25,000,000 iterations and sample every 5000 interations.

`bash 2_interval_run_script.sh`

## Output likelihood, block number and map length of MCMC chain
`bash 3_summarize_MCMC_chain.sh`

## Read the MCMC summaries into R, output summary statistics

In addition to calculating summary statistics, we use the Gelman diagnostic to check for chain convergence.
The script outputs boxplots of the Gelman point estimates for likelihood, block number and map length.

`bash 4_MCMC_convergence_check.R`


For the SDR, calculate Rho only in males and then for the sex-averaged recombination rate, do 2/3*(male recombination
rate).









