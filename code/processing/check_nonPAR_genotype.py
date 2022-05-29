#!/usr/bin/python
import sys
import gzip

vcf_f = gzip.open(sys.argv[1], "rt")

print("chrom""\t""chromStart""\t""chromEnd""\t""107""\t""108""\t""109""\t""110""\t""111""\t""112""\t""113""\t""114""\t""115""\t""116")
for line in vcf_f:
    if not line.startswith("#"):
        line = line.rstrip().split()
        gen = line[9:]
        genotype = []
        for element in gen:
            genotype.append(element.split(":")[0])
        if '0/1' in genotype[5:] or '0|1' in genotype[5:]:
            print(line[0]+"\t"+str(int(line[1])-1)+"\t"+line[1]+"\t"+"\t".join(genotype))


