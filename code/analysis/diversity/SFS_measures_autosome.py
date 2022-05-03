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
from popgene.harmonic import harmonic_number
from popgene.slidingwindow import slidingWindow

# python SFS_measures.py ../data/allele_count/black.A.repeat.frq.count ../data/reference/black.repeat.depth.masked.fa Autosome
# python SFS_measures.py qR
# python SFS_measures.py ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count ../data/reference/black.repeat.depth.masked.fa nonPAR

# python SFS_measures.py ../data/allele_count/black.A.repeat.frq.count ../data/bed/autosomes.bed ../data/reference/black.repeat.depth.masked.fa 20 200000 200000
# python SFS_measures.py ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count ../data/bed/nonpar_scaf.bed ../data/reference/black.repeat.depth.masked.fa 15 100000 100000
#********************************************************************
# Reading scaffold length of autosomes into a dictionary
def scaf_len_dict(fasta_file_bed):
    len_f = open(fasta_file_bed, "r")
    fasta_index = {}
    for line in len_f:
        if not line.startswith("chrom"):
            line = line.rstrip().split()
            scaffold = line[0]
            len_scaffold = line[2]
            fasta_index[scaffold] = int(len_scaffold)
    return(fasta_index)
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
# #********************************************************************
def genotypeArray(win_start, win_end, snp_dict, scaffold):
    """
    Return a list of lists containing genotypes for a given window 
    in the form [[1,19],[18,2]] where each element is a SNP count.
    """
    genotype_array = [] # Genotype array for window
    het_l = []
    for snp in snp_dict[scaffold]:
        a1 = int(snp[3].split(":")[1]) # count of allele 1
        a2 = int(snp[4].split(":")[1]) # count of allele 2
        pos = int(snp[0])
        if not 0 in [a1, a2] and pos >= win_start and pos < win_end:
            het = pi_fun((a1+a2), a1, a2)[2]
            het_l.append(het)
            genotype_array.append([a1, a2])
    return(genotype_array, het_l)
# #********************************************************************
def sfs_measures(win_start, win_end, genotype_array, het_l, num_chr):
    tajimad_l = {}
    midpoint = (win_start + win_end)/2
    sfs_dict = folded_sfs(genotype_array)
    sfs_list = []
    for key, value in sorted(sfs_dict.items(), key=lambda item: int(item[0])):
        sfs_list.append(value)  
    tajimad_l[midpoint] = TajimaD(sfs_list, num_chr, sum(het_l))
    return(sum(het_l), theta_w(len(het_l), num_chr), midpoint, tajimad_l[midpoint][0], tajimad_l[midpoint][1], tajimad_l[midpoint][2],tajimad_l[midpoint][1]-sum(het_l), len(sfs_list))
# #********************************************************************
def resampling_fun(het_array, len_array, reps):
    resampling_list = [sum(het_array)]
    for i in range(0, reps):
        d = list(np.random.randint(low=0, high=len_array, size=len_array))
        overlap = itemgetter(*d)(het_array)
        resampling_list.append(sum(overlap))
    return(np.mean(resampling_list), np.percentile(resampling_list, 2.5), np.percentile(resampling_list, 97.5))
# #********************************************************************
snp_dict = alleleCountdict(sys.argv[1])
# len_dict = scaf_len_dict(sys.argv[2])
# ref = open(sys.argv[3], "r")
# fasta_dict = readfasta(ref)
num_chr = int(sys.argv[2])
# winsize = int(sys.argv[5])
# step = int(sys.argv[6])
window_file = open(sys.argv[3], "r")
region = sys.argv[4]

if region == "all":
    intervals = {}
    for line in window_file:
        line = line.rstrip().split()
        if line[0] in intervals.keys():
            intervals[line[0]].append([line[1], line[2], line[3]])
        else:
            intervals[line[0]] = [[line[1], line[2], line[3]]]
else:
    intervals = {}
    for line in window_file:
        line = line.rstrip().split()
        if line[0] in intervals.keys():
            intervals[line[0]].append([line[1], line[2], line[7]])
        else:
            intervals[line[0]] = [[line[1], line[2], line[7]]]    

# print(intervals)

# print(snp_dict)
# print(len_dict)
header = ["Scaffold", "Start", "End", "midpont", "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resampl_CI_low", "pi_resamp_CI_up"]
print("\t".join(header))
for scaffold in snp_dict.keys():
    #intervals = slidingWindow(len_dict[scaffold], winsize, step)
    for window in intervals[scaffold]:
    #for window in intervals:
        window_start = int(window[0])
        window_end = int(window[1])
        window_length = int(window[2])
        midpoint = (window_start + window_end)/2
        geno_array = genotypeArray(window_start, window_end, snp_dict, scaffold)[0]
        # print(geno_array)
        if len(geno_array) != 0:
            het_list = genotypeArray(window_start, window_end, snp_dict, scaffold)[1]
            # print(len(het_list))
            # conf_interval = resampling_fun(het_array = het_list, len_array = len(het_list), reps = 99)
            sfs = sfs_measures(window_start, window_end, geno_array, het_list, num_chr)
            pi = sfs[0]
            theta = sfs[1]
            midpoint = sfs[2]
            td = sfs[3]
            pi_td = sfs[4]
            theta_td = sfs[5]
            diff = sfs[6]
            out = [scaffold, str(window_start), str(window_end), str(midpoint), str(pi), str(theta), str(td), str(window_length)]
            print("\t".join(out))
        else:
            # midpoint = (window_start + window_end)/2
            # w_l = seq_length(fasta_dict[scaffold][window_start:window_end])
            out = [scaffold, str(window_start), str(window_end), str(midpoint), 'NA', 'NA', 'NA', str(window_length)]
            print("\t".join(out))

