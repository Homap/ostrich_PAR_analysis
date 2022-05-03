#!/usr/bin/bash

while read -r ind
do
mkdir -p black/${ind}
echo black/${ind}
cd black/${ind}
sbatch ../../psmc_CI.sh ../${ind}.autosome.psmcfa ${ind}.split.autosome.psmcfa \
../${ind}.autosome.psmcfa.out ${ind}.autosome.combined.psmc \
${ind}.autosome.combined.graph.psmc
cd ../..
done < ../../data/samples/black.args

while read -r ind
do
mkdir -p blue/${ind}
echo blue/${ind}
cd blue/${ind}
sbatch ../../psmc_CI.sh ../${ind}.autosome.psmcfa ${ind}.split.autosome.psmcfa \
../${ind}.autosome.psmcfa.out ${ind}.autosome.combined.psmc \
${ind}.autosome.combined.graph.psmc
cd ../..
done < ../../data/samples/blue.args

while read -r ind
do
mkdir -p red/${ind}
echo red/${ind}
cd red/${ind}
sbatch ../../psmc_CI.sh ../${ind}.autosome.psmcfa ${ind}.split.autosome.psmcfa \
../${ind}.autosome.psmcfa.out ${ind}.autosome.combined.psmc \
${ind}.autosome.combined.graph.psmc
cd ../..
done < ../../data/samples/red.args


