#!/bin/bash -l
#SBATCH -A snic2022-22-149 
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J ldhat
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

split_n=$1
element_n=$1

./LDhat/complete -n 20 -rhomax 100 -n_pts 101 -theta 0.002 -split $split_n -element $element_n

