#!/usr/bin/python
from __future__ import division
import sys
import gzip
import numpy as np


def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

# Read a vcf file
vcf_file = gzip.open(sys.argv[1], "rb")
# Parse the genotype field
sample_names = []
geno_position = []
allele_bases = []
# kb = 1000
for line in vcf_file:
	line = line.rstrip().split()
	if line[0].startswith("#CHROM"):
		sample_names = line[9:]
	if not line[0].startswith("#"):
		if line[0] == "superscaffold36":
			geno_position.append(str(int(line[1])))
			allele_bases.append("".join([line[3], line[4]]))

# print(geno_position)
# print(allele_bases)

# Read the fastPHASE output
fastphase_out = open(sys.argv[2], "r")

fastphase_matrix = []
for line in fastphase_out:
	ind_haps = []
	if line.startswith("ID"):
		ind_haps.append(next(fastphase_out).rstrip().split())
		ind_haps.append(next(fastphase_out).rstrip().split())
#		print(ind_haps)
		fastphase_matrix.append(ind_haps)

# print(fastphase_matrix)
#print(len(fastphase_matrix))

fastphase_matrix_array = np.array([np.array(xi) for xi in fastphase_matrix])

number_inds = fastphase_matrix_array.shape[0]
number_snps = fastphase_matrix_array.shape[2]
# print(str(number_inds)+" "+str(number_snps)+" "+"1")
print(str(number_inds)+" "+str(len(geno_position[0:1000]))+" "+"1")
for ind_index, ind_item in enumerate(fastphase_matrix_array):
	# print(ind_index, ind_item)
	seq_hap1 = []
	seq_hap2 = []
	genotypes = zip(ind_item[0], ind_item[1])
	# print(len(genotypes[0:1000]))
	for geno_index, geno in enumerate(genotypes[0:1000]):
		allele1 = list(geno)[0]
		allele2 = list(geno)[1]
		if allele1 == allele2 == '0':
			#print(allele1, allele2, allele_bases[geno_index][0], allele_bases[geno_index][0])
			seq_hap1.append(allele_bases[geno_index][0])
			seq_hap2.append(allele_bases[geno_index][0])
		elif allele1 == allele2 == '1':
			#print(allele1, allele2, allele_bases[geno_index][1], allele_bases[geno_index][1])
			seq_hap1.append(allele_bases[geno_index][1])
			seq_hap2.append(allele_bases[geno_index][1])				
		elif allele1 == '0' and allele2 == '1':
			#print(allele1, allele2, allele_bases[geno_index][0], allele_bases[geno_index][1])
			seq_hap1.append(allele_bases[geno_index][0])
			seq_hap2.append(allele_bases[geno_index][1])
		elif allele1 == '1' and allele2 == '0':
			#print(allele1, allele2, allele_bases[geno_index][1], allele_bases[geno_index][0])
			seq_hap1.append(allele_bases[geno_index][1])
			seq_hap2.append(allele_bases[geno_index][0])				
		else:
			raise Exception('Allele encoding is incorrect!')
	print(">"+sample_names[ind_index]+"_1")
	# print(len("".join(seq_hap1)))
	for seq in chunks("".join(seq_hap1), 1000):
		print(seq)
	print(">"+sample_names[ind_index]+"_2")
	for seq in chunks("".join(seq_hap2), 1000):
		print(seq)

# print(len(geno_position))

relative_geno = [str(int(pos) - int(geno_position[0]) + 1) for pos in geno_position[0:1000]]
# print(relative_geno)

L = int(geno_position[1000]) - int(geno_position[0]) + 1

locs = open("sites.PARphased.locs", "w")
# locs.write(str(number_snps)+" "+str(L)+" "+"L"+"\n")
locs.write(str(len(relative_geno))+" "+str(L)+" "+"L"+"\n")
locs.write(" ".join(relative_geno)+"\n")


vcf_file.close()
fastphase_out.close()
locs.close()