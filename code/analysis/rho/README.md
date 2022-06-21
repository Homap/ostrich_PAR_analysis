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

## Check for convergence of MCMC 

In addition to calculating summary statistics, we use the Gelman diagnostic to check for chain convergence.
The script outputs boxplots of the Gelman point estimates for likelihood, block number and map length.

`Rscript 4_MCMC_convergence_check.R`

![Chr4 Gelman diagnostic](../../../data/rho/ldhat_mcmc/chr4_gelman.pdf) <br>
![Chr5 Gelman diagnostic](../../../data/rho/ldhat_mcmc/chr5_gelman.pdf) <br>
![PAR Gelman diagnostic](../../../data/rho/ldhat_mcmc/par_gelman.pdf) <br>
![nonPAR Gelman diagnostic](../../../data/rho/ldhat_mcmc/nonpar_gelman.pdf) <br>

Gelman diagnostic value of about 1 indicates convergence of MCMC chains. We observe that for our MCMC chains.

## Extract rho per window and sites from the *rates* output of *interval*

Using the *stat* program in LDhat to summarise the *rates* output

`bash 5_stat_LDhat.sh`

The per site rate from the *interval* as reported in the *rates.txt* file is per kb (kilo base)
since the *locs* input was in kb. To convert that per site, the rate for each pair of site must 
be divided by 1000. This gives per bp rate of recombination. 

The python script `stat_rho.py` takes the output of LDhat *stat* program with the list of 
positions of each window. It outputs two files: One file containing the map length for each window
and one file containing the per site rho for every pair of positions.

`bash 6_map_length_rho_per_site.sh`

Two outputs for each scaffold are as follow:

- Map length (total rho for each window)

| Scaffold | Locus_start | Locus_end | Mean_rho | Median | L95 | U95 |
| -------- | ----------- | --------- | -------- | ------ | --- | --- |
| superscaffold26 | 1306 | 190500 | 199.28876 | 199.15247 | 189.58553 | 210.24434 |

- Per site rho (per kbp rho for each pair of SNPs)

| Scaffold | Locus_start | Locus_end | Mean_rho | Median | L95 | U95 |
| -------- | ----------- | --------- | -------- | ------ | --- | --- |
| superscaffold26 | 1306 | 1622 | 8.09272 | 7.67439 | 5.81891 | 13.38059 |

To convert to bp divide Mean_rho by 1000. To obtain the map length from the per site rho, in *R* you can do:

`map_length <- sum((persite$Locus_end - persite$Locus_start)*(persite$Mean_rho/1000))`

## Map length and per site rho with Z chromosome coordinates
We concatenate the map and rates files and create one for each with the Z coordinates.

Order of Z scaffolds

| Number | Scaffold |
| ------ | -------- |
| 1 | superscaffold26 |
| 2 | superscaffold54 |
| 3 | superscaffold35 |
| 4 | superscaffold36 |
| 5 | superscaffold62 |
| 6 | superscaffold67 |
| 7 | superscaffold69-1 |
| 8 | superscaffold93 |
| 9 | superscaffold63 |
| 10| superscaffold88 |
| 11| superscaffold83 |
| 12| superscaffold92 |

`bash 7_rho_Z_coordinates.sh &`

## Use rho and recombination rate to obtain Ne for the PAR



Also in 200 Kb windows with 50 Kb overlap



For the SDR, calculate Rho only in males and then for the sex-averaged recombination rate, do 2/3*(male recombination
rate).






