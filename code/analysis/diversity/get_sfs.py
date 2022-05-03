#!/usr/bin/python

import sys
import numpy as np
from collections import Counter
import math
from operator import itemgetter
from fasta import readfasta
from popgene.basecount import seq_length
from popgene.tajimad import TajimaD
from popgene.theta import theta_w
from popgene.pi import ncr
from popgene.pi import pi_fun
from popgene.foldedSFS import folded_sfs

#********************************************************************
def alleleCountdict(snp_file_name):
    """
    Read the allele count file obtained from vcftools into a dictionary
    """
    snp_file = open(snp_file_name, "r")
    snp_dict = {}
    for line in snp_file:
        if not line.startswith("CHROM"):
            line = line.rstrip().split()
            scaf = line[0]
            info = line[1:]
            if scaf in snp_dict.keys():
                snp_dict[scaf].append(info)
            else:
                snp_dict[scaf] = [info]
    snp_file.close()
    return(snp_dict)
#********************************************************************
def genotypeArray(snp_dict):
    """
    Return a list of lists containing genotypes for a given window 
    in the form [[1,19],[18,2]] where each element is a SNP count.
    """
    genotype_array = [] # Genotype array for window
    for scaffold in snp_dict.keys():
        for snp in snp_dict[scaffold]:
            a1 = int(snp[3].split(":")[1]) # count of allele 1
            a2 = int(snp[4].split(":")[1]) # count of allele 2
            if not 0 in [a1, a2]:
                genotype_array.append([a1, a2])
    return(genotype_array)
#********************************************************************

snp_dict = alleleCountdict(sys.argv[1])
geno_array = genotypeArray(snp_dict)
sfs_dict = folded_sfs(geno_array)
sfs_list = []
for key, value in sorted(sfs_dict.items(), key=lambda item: int(item[0])):
    sfs_list.append(value)  

species = sys.argv[2]
chr = sys.argv[3]

print("count"+"\t"+species+"_"+chr+"\t"+"chr_norm")
for index, item in enumerate(sfs_list):
    print(str(index+1)+"\t"+str(item)+"\t"+str(item/sum(sfs_list)))



