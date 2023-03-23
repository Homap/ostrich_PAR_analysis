#!/usr/bin/bash

mkdir -p ../../figures

Rscript FigureS1.R

Rscript FigureS2.R 

Rscript FigureS3.R

Rscript Figure2.R

Rscript Figure3.R

Rscript Figure4.R

# intro_figure.R contains the recombination frequency plot and the plot showing genes along the chromosome. 
# To produce the final figure, after plotting each in R, I assembled them in inkscape.
# The LD matrix in Figure S4 was produced by th origram [LDBlockShow](https://github.com/BGI-shenzhen/LDBlockShow)
