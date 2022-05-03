#!/usr/bin/python
import sys

genotypes = open("../data/bed/black.nonPAR.genotypes.norepeats.bed")

print("scaf"+"\t"+"pos"+"\t"+"male_het"+"\t"+"female_het"+"\t"+"residual")
for line in genotypes:
    line = line.rstrip().split()
    gen = line[3:13]
    genotype = []
    for element in gen:
        genotype.append(element.split(":")[0])
        male_genotype = genotype[0:5]
        female_genotype = genotype[5:10]
    #print(male_genotype)
    #print(female_genotype)
    male_het_prop = (male_genotype.count("0/1") + male_genotype.count("0|1"))/5
    female_het_prop = (female_genotype.count("0/1") + female_genotype.count("0|1"))/5
    residual = male_het_prop - female_het_prop
    #print(male_het_prop, female_het_prop)
    print(line[0]+"\t"+line[2]+"\t"+str(male_het_prop)+"\t"+str(female_het_prop)+"\t"+str(residual))

        
