# Coalscent simulation for the whole PAR and at the PAR boundary

## Create conda environment for msprime from the yml file
conda activate par_simulation

## Install msprime
conda install pandas=1.4.3
conda install matplotlib=3.5.2
conda install -c conda-forge msprime=1.2.0

## Run the simulation for the whole PAR
python makeData_fullPAR_plot.py
