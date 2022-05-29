#!/usr/bin/python
import sys
import gzip

# Homa Papoli

# Extract genotypes from VCF file and if heterozygous (0/1 or 0|1) found in females, the SNP is recorded in the output.

vcf_f = gzip.open(sys.argv[1], "rt")

for line in vcf_f:
    if not line.startswith("#"):
        line = line.rstrip().split()
        gen = line[9:]
        genotype = []
        for element in gen:
            genotype.append(element.split(":")[0])
        print(line[0]+"\t"+str(int(line[1])-1)+"\t"+line[1]+"\t"+"\t".join(genotype))

        
