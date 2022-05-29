#!/usr/bin/python

import sys
import gzip
import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd

# Get all the scaffolds in the VCF with length

# Autosome VCF background
# Subsample the all sites VCF
# bcftools view ../data/vcf/black.all.vcf.gz | vcfrandomsample -r 0.012 > black.subset.vcf
# Get a list of sites with scaffold and position to filter out for depth
# Get a list of sites with repeats to filter out for repeats
# Create a masked fasta with all above as N
# Create a window file with number of filtered background sites

# Using the above filtered sites:
# Create a window of a given number of intergenic sites, for example 100K
# Create a window of a given number of intronic sites, for example 100K
# Create a window of a given number of 4-fold sites, for example 100K
# Create a window of a given number of 0-fold sites, for example 100K

# PAR background
# Get a list of sites with scaffold and position to filter out for depth
# Get a list of sites with repeats to filter out for repeats
# Create a masked fasta with all above as N
# Create a window file with number of filtered background sites 

# Using the above filtered sites:
# Create a window of a given number of intergenic sites, for example 100K
# Create a window of a given number of intronic sites, for example 100K
# Create a window of a given number of 4-fold sites, for example 100K
# Create a window of a given number of 0-fold sites, for example 100K

# nonPAR background
# Get a list of sites with scaffold and position to filter out for depth
# Get a list of sites with repeats to filter out for repeats
# Create a masked fasta with all above as N
# Create a window file with number of filtered background sites

# Using the above filtered sites:
# Create a window of a given number of intergenic sites, for example 100K
# Create a window of a given number of intronic sites, for example 100K
# Create a window of a given number of 4-fold sites, for example 100K
# Create a window of a given number of 0-fold sites, for example 100K

# Autosome PAR and nonPAR
# 

# Open autosome, PAR and non-PAR bed files

# Open the zipped VCF

# chrom   chromStart      chromEnd
# superscaffold26 0       25310599
# superscaffold54 0       16379243
# superscaffold35 0       4625539
# superscaffold36 3524263 9394175

# chrom   chromStart      chromEnd
# superscaffold36 0       3516673
# superscaffold62 0       2917291
# superscaffold67 0       5300260
# superscaffold69-1       0       5978518
# superscaffold93 0       4983591
# superscaffold63 0       1692925
# superscaffold88 0       624114
# superscaffold83 0       782506
# superscaffold92 0       2882843

# Autosome mean coverage is about 32X, therefore I set the max depth to 64X
# PAR mean coverage is about 30X, max depth to 60X
# nonPAR mean coverage is 24X, max depth to 50X

PAR_c = {"superscaffold26": [0, 25310599], "superscaffold54": [0,16379243], "superscaffold35": [0, 4625539], "superscaffold36": [3524263, 9394175]}
nonPAR = {"superscaffold36": [0, 3516673], "superscaffold62": [0, 2917291], "superscaffold67": [0, 5300260], "superscaffold69-1": [0, 5978518], 
            "superscaffold93": [0, 4983591], "superscaffold63": [0, 1692925], "superscaffold88": [0, 624114], "superscaffold83": [0, 782506], "superscaffold92": [0, 2882843]}

f = open(sys.argv[1])

#print("CHROM"+"\t"+"START"+"\t"+"END"+"\t"+"DEPTH")

for line in f:
    if not line.startswith("CHROM"):
        line = line.rstrip().split()
        scaf = line[0]
        pos = int(line[1])
        depth = float(line[2])
        if scaf in PAR_c.keys(): # If scaffold in in the PAR
            if scaf == "superscaffold36": # If the PAR is superscaffold36
                if pos >= PAR_c[scaf][0]+1 and pos <= PAR_c[scaf][1]: 
                    if depth < 10 or depth > 60:
                       print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))
                else:
                    if depth < 10 or depth > 50:
                       print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))
            elif scaf == "superscaffold54":
                if pos >= PAR_c[scaf][0]+1 and pos <= PAR_c[scaf][1]:
                    if depth < 10 or depth > 60:
                       print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))     
                else:
                    if depth < 10 or depth > 64:
                       print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth)) 
        elif scaf in nonPAR.keys():
            if depth < 10 or depth > 50:
                print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))                                                      
        else:
            if depth < 10 or depth > 64:
                print(scaf+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+str(depth))

f.close()



