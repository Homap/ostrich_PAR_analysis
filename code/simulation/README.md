# Coalscent simulation for the whole PAR and at the PAR boundary

We evaluate whether the observed patterns of genetic variation and male-female differentiation across the PAR and close to the SDR boundary are consistent with predictions from neutral genetic theory using coalescent simulations based on (Kirkpatrick, Guerrero, & Scarpino, 2010).

## Create conda environment for msprime from the yml file
The simulation is done using the `msprime` program. To run the simulation scripts, the easiest way is to create a conda and activate a conda environment as follows:

`conda create --name par_simulation_env python=3.10` <br>
`conda activate par_simulation_env` <br>

## Install dependecies and msprime
`conda install pandas=1.4.3` <br>
`conda install matplotlib=3.5.2` <br> 
`conda install -c conda-forge msprime=1.2.0` <br>

## Run the simulation for the whole PAR
python makeData_fullPAR_plot.py
