#!/usr/bin/python
import sys
from fasta import readfasta

# Calculates the length of each scaffold and write out into a bed file
# for the autosomes.
# Written by Homa Papoli, July 2020
# Run:
# python scaf_length_to_bed.py fasta.fa autosome_list


f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")
f1_fasta = readfasta(f1)


scaf_list = []
for line in f2:
	line = line.strip("\n")
	scaf_list.append(line)

print("chrom"+"\t"+"chromStart"+"\t"+"chromEnd")
for scaffold in scaf_list:
	print(scaffold+"\t"+"0"+"\t"+str(len(f1_fasta[scaffold])))
