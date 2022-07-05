#!/usr/bin/bash

for window in 100Kb 200Kb 1000Kb
do
echo $window
w_size=$(echo $window | sed 's/Kb/000/g')
python ../../processing/sliding_window.py ../../../data/bed/z_chrom.bed ${w_size} > ../../../data/sliding_window/Z.${w_size}.bed
python ../../processing/sliding_window.py ../../../data/bed/z_par.bed ${w_size} > ../../../data/sliding_window/PAR.${w_size}.bed
python ../../processing/sliding_window.py ../../../data/bed/z_nonpar.bed ${w_size} > ../../../data/sliding_window/nonPAR.${w_size}.bed

python ../../processing/sliding_window.py ../../../data/bed/z_scaf.bed ${w_size} > ../../../data/sliding_window/Z.scaf.${w_size}.bed
python ../../processing/sliding_window.py ../../../data/bed/par_scaf.bed ${w_size} > ../../../data/sliding_window/PAR.scaf.${w_size}.bed
python ../../processing/sliding_window.py ../../../data/bed/nonpar_scaf.bed ${w_size} > ../../../data/sliding_window/nonPAR.scaf.${w_size}.bed

python ../../processing/sliding_window.py ../../../data/lastz/gg_chr4_ostrich.bed ${w_size} > ../../../data/sliding_window/chr4.${w_size}.bed
python ../../processing/sliding_window.py ../../../data/lastz/gg_chr5_ostrich.bed ${w_size} > ../../../data/sliding_window/chr5.${w_size}.bed
done

