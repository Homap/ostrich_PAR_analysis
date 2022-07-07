# Scripts used for processing, analysing and visualising genomic data

- **processing** <br>
Codes for filtering the VCF files, generating sliding windows and conversion of scaffold coordinates to chromosome coordinates.

- **analysis** <br>
Codes for obtaining various population genetics measure including LD, rho, genetic diversity. 
    * **coverage** : Median coverage per individual, PAR-nonPAR boundary coverage
    * **diversity** : π, θ, Tajima's D, folded SFS
    * **fst** : Male - Female Fst
    * **geneticmap** : Recombination rate from genetic map in windows
    * **genomic_features** : Count of CDS and intergenic regions per windows
    * **lastz** : Alignment of ostrich scaffolds to chicken genome assembly 
    * **ld** : LD decay and pairwise LD per window
    * **rho** : Population scaled recombination rate (ρ) 

- **statistical_analyses** <br>
Code for statistical analyses including GLM using the outputs from **analysis**.

- **simulation** <br>
Code for coalescent simulation presented in the study.

- **figures** <br>
Code to generate the figures presented in the manuscript.


