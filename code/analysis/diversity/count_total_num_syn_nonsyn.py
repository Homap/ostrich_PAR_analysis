#!/usr/bin/python
from __future__ import division
import sys
import math

"""
Written by Homa Papoli
11 June 2020
./count_total_num_syn_nonsyn.py ../src/snpEff/scaf26.annotate.txt 
../src/snpEff/scaf26.snp.txt > ../src/snpEff/scaf26.het.txt
"""

# To run:
# python count_total_num_syn_nonsyn.py ../data/gff/genome.annotated.txt \
# ../data/allele_count/black.nonPAR.filtered.adjusted.frq.count 15 > ../data/allele_count/black.nonPAR.annotated.pi.txt
# superscaffold54 848     849     2       10      A:4     G:6
# superscaffold54 2811    2812    2       10      G:2     A:8
# **********************************
# Function denfinitions
# **********************************
# Function for N choose r
# **********************************
def ncr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
# **********************************
# Function for Pairwise differences
# **********************************
def pi_fun(nchr, allele1, allele2):
	P_diff = allele1*allele2
	NCR = ncr(nchr, 2)
	pi = round(P_diff/NCR, 4)
	return (P_diff, NCR, pi) 

annt_list = open(sys.argv[1], "r")
snp_list = open(sys.argv[2], "r")
num_chr =int(sys.argv[3])

#print("snp dict")
snp_dict = {}
for line in snp_list:
	if not line.startswith("CHROM"):
		line = line.strip("\n").split()
		scaf_pos = line[0] + "_" + line[1]
		allele1 = line[4].split(":")[1]
		allele2 = line[5].split(":")[1]
		# print(allele1)
		if not '0' in [allele1, allele2]:
			snp_dict[scaf_pos] = allele1
#print(snp_dict)

print("scaffold"+"\t"+"position"+"\t"+"gene"+"\t"+"synonymous"+"\t"+"missense"+"\t"+"nonsense"+"\t"+"count_allele1"+"\t"+"count_allele2"+"\t"+"total_pairs_diff"+"\t"+"total_pairs_num"+"\t"+"pi")
gene_annt_count_syn = {}
gene_annt_count_nonsyn = {}
genelist = []
for line in annt_list:
	line = line.strip("\n").split()
	#print(line)
	if not line[1] == "mRNA_id":
		gene = line[1]
		syn = float(line[5])
		miss = float(line[6])
		nons = float(line[7])
		scaffold = line[0]
		pos = int(line[4]) + 1
		scaf_pos = scaffold + "_" + str(pos)
		if scaf_pos in snp_dict.keys():
			p_diff = str(pi_fun(num_chr, int(snp_dict[scaf_pos]), num_chr-int(snp_dict[scaf_pos]))[0])
			ncr1 = str(pi_fun(num_chr, int(snp_dict[scaf_pos]), num_chr-int(snp_dict[scaf_pos]))[1])
			pi = str(pi_fun(num_chr, int(snp_dict[scaf_pos]), num_chr-int(snp_dict[scaf_pos]))[2])
			print(scaf_pos.split("_")[0]+"\t"+scaf_pos.split("_")[1]+"\t"+gene+"\t"+str(syn)+"\t"+str(miss)+"\t"+str(nons)+"\t"+snp_dict[scaf_pos]+"\t"+str(num_chr-int(snp_dict[scaf_pos]))+"\t"+p_diff+"\t"+ncr1+"\t"+pi)
		else:
			print(scaf_pos.split("_")[0]+"\t"+scaf_pos.split("_")[1]+"\t"+gene+"\t"+str(syn)+"\t"+str(miss)+"\t"+str(nons)+"\t"+"0"+"\t"+"10"+"\t"+"0"+"\t"+"NA"+"\t"+"0")			
		nonsyn = miss + nons
		genelist.append(gene)
		if not gene in gene_annt_count_syn.keys():
			gene_annt_count_syn[gene] = [syn]
		else:
			gene_annt_count_syn[gene].append(syn)
		if not gene in gene_annt_count_nonsyn.keys():
			gene_annt_count_nonsyn[gene] = [nonsyn]
		else:
			gene_annt_count_nonsyn[gene].append(nonsyn)

#for geneid in list(set(genelist)):
#	print(geneid.split("-")[1]+"\t"+str(sum(gene_annt_count_syn[geneid]))+"\t"+str(sum(gene_annt_count_nonsyn[geneid])))

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone.frq.count | awk '$11!=0' > annotate_old_2/crotone_annotated_pi.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone.frq.count  > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone_annotated_pi.txt

# for geneid in list(set(genelist)):
# 	scaf_pos = geneid.split(":")[1]
# 	if scaf_pos in snp_dict.keys():
# 		print(scaf_pos.split("_")[0]+"\t"+scaf_pos.split("_")[1]+"\t"+geneid.split("-")[1]+"\t"+str(sum(gene_annt_count_syn[geneid]))+"\t"+str(sum(gene_annt_count_nonsyn[geneid]))+"\t"+snp_dict[scaf_pos]+"\t"+str(10-int(snp_dict[scaf_pos])))
# 	else:
# 		print(scaf_pos.split("_")[0]+"\t"+scaf_pos.split("_")[1]+"\t"+geneid.split("-")[1]+"\t"+str(sum(gene_annt_count_syn[geneid]))+"\t"+str(sum(gene_annt_count_nonsyn[geneid]))+"\t"+"10"+"\t"+"0")		
# #print(sum(gene_annt_count_syn["rna-XM_009667647.1"]))
# #print(sum(gene_annt_count_nonsyn["rna-XM_009667647.1"]))

