#!/usr/bin/env python
import gzip
import sys
import io
import os
import numpy as np
import textwrap
import argparse
from argparse import ArgumentParser, HelpFormatter
import subprocess
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import math

class RawFormatter(HelpFormatter):
	def _fill_text(self, text, width, indent):
		return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

program_descripton = f'''
	Title: vcf_to_genotype_fasta.py
	Version: v2.0
	Date: 10 June 2022
	Author: Homa Papoli Yazdi
	
	To run the script, create a conda environment using the ldhat_env.yml file as below:
	conda env create --file ldhat_env.yml

	Activate the conda environment:
	conda activate ldhat

	You can run the script as follows for 20 SNPs, 19 overlapping for superscaffold36 and 10 windows:
	./vcf_to_ldhat_out.py black.PAR.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz par_scaf.bed 20 19 superscaffold36 10

	Description: Take a gzipped VCF file as an input, the region and the number of SNPs 
	to extract from VCF and outputs the LDhat inputs, the sites and locs files.

	**sites file**
	4 10 2
	>GenotypeA
	122110?000
	>GenotypeB
	1111201100
	>GenotypeC
	011111?112
	>GenotypeD
	2112210100

	**locs file**
	10 1200 L
	1 57 180 187 223 250 438 509 878 1034

	For genotype/unphased data, the convention used is 0 and 1 for the two homozygotes,
	and 2 for heterozygotes. ? is used for missing data. The first line of the sites
	file is as follows: 4: number of sequences, 10: number of SNPs, 2: unphased 
	The first line of the locs file is as follows: 10: number of SNPs, 1200: total
	length of the sequence, L: using cross-over model (the other model is the gene
	coversion model)

	List of outputs of the Python script:
	Name examples below are for superscaffold36, ran for 20 SNPs, 10 overlaps and 1 window:
	1. sites input of LDhat pairwise called : superscaffold36.20.19.1.sites.txt
	2. locs input of LDhat pairwise called : superscaffold36.20.19.1.locs.txt
	3. pairwise rho called : superscaffold36.20.19.1.pairwise.rmin.txt
	4. rho ouput for all windows called : superscaffold36.20.19.1.rho.out"
	2. Original SNP positions called : superscaffold36.20.19.1.pos.txt

	List of outputs of the LDhat pairwise script:
	# Their descript can be found in LDhat manual.
	superscaffold36.20.19.1.ldhat.fit.txt
	superscaffold36.20.19.1.ldhat.freqs.txt
	superscaffold36.20.19.1.ldhat.new_lk.txt
	superscaffold36.20.19.1.ldhat.outfile.txt
	superscaffold36.20.19.1.ldhat.rmin.txt
	superscaffold36.20.19.1.ldhat.type_table.txt
	'''

# Set script's arguments

parser = ArgumentParser(description=program_descripton, formatter_class=RawFormatter)
parser.add_argument('vcf', help='gzipped VCF')
parser.add_argument('bed', help='scaffold length file in bed format')
parser.add_argument('Nsnps', help='number of SNPs', type=int)
parser.add_argument('Noverlap', help='number of overlapping SNPs, Nsnps=50 and Noverlap=49, then there are 49 overlapping SNPs', type=int)
parser.add_argument('chr', help='chromosome or scaffold', type=str)
parser.add_argument('Nwin', help='number of windows', type=int)

args = parser.parse_args()

#**********************************************************************************
# Define the sliding window function

def slidingWindow(sequence_l, winSize, step):
	""" Returns a generator that will iterate through
	the defined chunks of input sequence. Input
	sequence must be iterable."""

	# Verify the inputs
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step cannot be larger than winSize.")
	if winSize > sequence_l:
		pass
	# Pre-compute number of chunks to emit
	numOfChunks = ((int(sequence_l-winSize)/step))+1
	numOfChunks = int(numOfChunks) 

	for i in range(0, numOfChunks*step, step):
		yield i,i+winSize
	if sequence_l > numOfChunks*step:
		yield i+winSize, sequence_l

# Function to control output character length per line. Ldhat does not accept
# more than 2000 characters per line.
def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]
#**********************************************************************************
# Read scaffold length file into a dictionary
def main():
	scaf_len = {}
	with open(args.bed, 'r') as lenfile:
		next(lenfile)
		for line in lenfile:
			line = line.rstrip().split()
			scaf_len[line[0]] = int(line[2]) - int(line[1])
	#**********************************************************************************
	# Parse the genotype field
	# print("Parsing the VCF file")
	genotype_matrix = []
	sample_names = []
	geno_position = []
	with io.TextIOWrapper(gzip.open(args.vcf, 'r')) as vcf_file:
		for line in vcf_file:
			line = line.rstrip().split()
			if line[0].startswith("#CHROM"):
				sample_names = line[9:]
			if not line[0].startswith("#"):
				if line[0] == args.chr:
					geno_position.append(line[1])
					genotype_fields = line[9:]
					genotype_matrix.append(genotype_fields)

	# Create list of windows
	windows = slidingWindow(len(genotype_matrix), args.Nsnps, args.Nsnps-args.Noverlap)
	#**********************************************************************************
	# # rho_out_name: outputs the SNP position and rho for each window
	# rho_out_name = args.chr + "." + str(args.Nsnps) + "." + str(args.Noverlap) + "." + str(args.Nwin) + ".rho.out"
	# rho_out = open(rho_out_name, "w")
	# print("chr"+"\t"+"SNP1_pos"+"\t"+"SNP2_pos"+"\t"+"length"+"\t"+"rho"+"\t"+"rho_per_site"+"\t"+"Lk")
	# rho_out.write("chr"+"\t"+"SNP1_pos"+"\t"+"SNP2_pos"+"\t"+"length"+"\t"+"rho"+"\t"+"rho_per_site"+"\t"+"Lk"+"\n")

	# # print("Running pairwise program of LDhat")
	interval_counter = 0
	for interval_index, interval in enumerate(windows):
		#print(interval) # 0, 2000
		#print(interval_index)
		if interval_index < args.Nwin:
			start = list(interval)[0]
			end = list(interval)[1]
			print(start, end)
			for index, item in enumerate(genotype_matrix):
				#print(index, start)
				index = index + start
				#print("INDEX", index)#, start, end, item)
				if start <= index <= end and index != len(genotype_matrix): 
					# print(index, start, end)
					for index_ind, item_ind in enumerate(item):
						genotype = genotype_matrix[index][index_ind].split(":")[0]
						if genotype == "0/1" or genotype == "0|1" or genotype == "2":
							genotype_matrix[index][index_ind] = '2'
						elif genotype == "1/0" or genotype == "1|0" or genotype == "2":
							genotype_matrix[index][index_ind] = '2'
						elif genotype == "0/0" or genotype == "0|0" or genotype == "0":
							genotype_matrix[index][index_ind] = '0'
						elif genotype == "1/1" or genotype == "1|1" or genotype == "1":
							genotype_matrix[index][index_ind] = '1'
						elif genotype == "./." or genotype == ".|.":
							genotype_matrix[index][index_ind] = '?'
						else:
							raise Exception("Genotype field contains incorrect notation!")
				else:
					break

			genotype_np_array = np.array([np.array(xi) for xi in genotype_matrix[start:end]])
			# print("1", genotype_np_array)

			gen_array = genotype_np_array.transpose()
			# print("2", gen_array)

			gen_array_seq = np.apply_along_axis(lambda row: row.astype('|S1').tobytes().decode('utf-8'), axis=1,arr=gen_array)
			# print("3", gen_array_seq)

			interval_counter += 1
			sites = args.chr + "." + str(args.Nsnps) + "." + str(args.Noverlap) + "." + str(interval_counter) + ".sites.txt"
			# print(outname)
			with open(sites, 'w') as sites_out:
				sites_out.write(str(gen_array.shape[0])+" "+str(gen_array.shape[1])+" "+"2"+"\n")
				for index, ind in enumerate(gen_array_seq):
					sites_out.write(">"+sample_names[index]+"\n")
					for seq in chunks(gen_array_seq[index], args.Nsnps):
						sites_out.write(seq+"\n")

			L = int(geno_position[end-1]) - int(geno_position[start]) + 1
			# print(L)
			locs = args.chr + "." + str(args.Nsnps) + "." + str(args.Noverlap) + "." + str(interval_counter) + ".locs.txt"
			with open(locs, 'w') as locs_out:
				new_coord = [str(int(coord) - int(geno_position[start]) + 1) for coord in geno_position[start:end]]
				locs_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+"\n".join(new_coord)+"\n")

			original_pos = args.chr + "." + str(args.Nsnps) + "." + str(args.Noverlap) + "." + str(interval_counter) + ".pos.txt"
			with open(original_pos, 'w') as pos_out:
				pos_out.write(str(gen_array.shape[1])+" "+str(L)+" "+ "L"+"\n"+"\n".join(geno_position[start:end])+"\n")
	return 0

if __name__ == "__main__":
	main()