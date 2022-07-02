#!/usr/bin/python

import sys
import numpy as np


# Get all the scaffolds in the VCF with length

f = open(sys.argv[1])

#print("CHROM"+"\t"+"START"+"\t"+"END"+"\t"+"DEPTH")

for line in f:
    if not line.startswith("#CHROM"):
        line = line.rstrip().split()
        scaf = line[0]
        pos = int(line[1])        
        depth = float(np.mean([int(i) for i in line[2:]]))
        if depth < 5 or depth > 70:
            print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))

f.close()



