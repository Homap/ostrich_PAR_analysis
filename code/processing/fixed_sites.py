#!/usr/bin/python

import sys


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

snp_dict = alleleCountdict(sys.argv[1])

for key in snp_dict:
    for snp in snp_dict[key]:
        allele1 = snp[3].split(":")[1]
        allele2 = snp[4].split(":")[1]
        if allele1 == "0" and allele2 == "20":
            print(key+"\t"+snp[0])
        elif allele1 == "20" and allele2 == "0":
            raise Exception("Error! fixed for reference cannot be in VCF")
        else:
            pass
        # print(allele1, allele2)
#print(snp_dict)