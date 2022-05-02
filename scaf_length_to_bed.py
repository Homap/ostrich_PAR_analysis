#!/usr/bin/python
import sys
from fasta import readfasta

#********************************************************************************
# Calculates the length of each scaffold and write out into a bed file
# for the Z, PAR and non-PAR.
#********************************************************************************
# Written by Homa Papoli, July 2020
# Run:
# python scaf_length_to_bed.py fasta.fa z_scaf.bed par_scaf.bed nonpar_scaf.bed
#********************************************************************************

f1 = open(sys.argv[1], "r")
f1_fasta = readfasta(f1)

# List of Z scaffolds
scaffold_list = ["superscaffold26","superscaffold54","superscaffold35",
"superscaffold36","superscaffold62","superscaffold67","superscaffold69-1",
"superscaffold93","superscaffold63","superscaffold88","superscaffold83",
"superscaffold92"]

par_list = ["superscaffold26","superscaffold54","superscaffold35",
"superscaffold36"]

nonpar_list = ["superscaffold36","superscaffold62","superscaffold67","superscaffold69-1",
"superscaffold93","superscaffold63","superscaffold88","superscaffold83",
"superscaffold92"]

z_bed = open(sys.argv[2], "w") # Open an output file to write the bed file into it.
z_bed.write("chrom"+"\t"+"chromStart"+"\t"+"chromEnd"+"\n")
for scaffold in scaffold_list:
	if scaffold == "superscaffold54":
		z_bed.write(scaffold+"\t"+"0"+"\t"+"16379243"+"\n")
	else:
		z_bed.write(scaffold+"\t"+"0"+"\t"+str(len(f1_fasta[scaffold]))+"\n")

par_bed = open(sys.argv[3], "w")
par_bed.write("chrom"+"\t"+"chromStart"+"\t"+"chromEnd"+"\n")
for par_scaf in par_list:
	if par_scaf == "superscaffold54":
		par_bed.write(par_scaf+"\t"+"0"+"\t"+"16379243"+"\n")
	elif par_scaf == "superscaffold36":
		par_bed.write(par_scaf+"\t"+"3524263"+"\t"+str(len(f1_fasta[par_scaf]))+"\n")
	else:
		par_bed.write(par_scaf+"\t"+"0"+"\t"+str(len(f1_fasta[par_scaf]))+"\n")

nonpar_bed = open(sys.argv[4], "w")
nonpar_bed.write("chrom"+"\t"+"chromStart"+"\t"+"chromEnd"+"\n")
for nonpar_scaf in nonpar_list:
	if nonpar_scaf == "superscaffold36":
		nonpar_bed.write(nonpar_scaf+"\t"+"0"+"\t"+"3516673"+"\n")
	else:
		nonpar_bed.write(nonpar_scaf+"\t"+"0"+"\t"+str(len(f1_fasta[nonpar_scaf]))+"\n")

