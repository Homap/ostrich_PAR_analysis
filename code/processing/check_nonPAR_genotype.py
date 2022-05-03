#!/usr/bin/python
import sys
import gzip

vcf_f = gzip.open(sys.argv[1], "rt")

for line in vcf_f:
    if not line.startswith("#"):
        line = line.rstrip().split()
        gen = line[9:]
        genotype = []
        for element in gen:
            genotype.append(element.split(":")[0])
        if '0/1' in genotype[5:] or '0|1' in genotype[5:]:
            print(line[0]+"\t"+str(int(line[1])-1)+"\t"+line[1]+"\t"+"\t".join(genotype))

        
