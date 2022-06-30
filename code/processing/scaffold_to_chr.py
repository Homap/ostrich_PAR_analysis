#!/usr/bin/python3
import sys
from collections import defaultdict

#************************************************************************************
# Written by Homa Papoli Yazdi
#************************************************************************************

# The input file should contain at least three columns: Chromosome, start, end

def usage():
	msg = '''converts scaffold coordinates of ostrich Z chromosome into chromosome coordinates based on linkage map Yazdi and Ellegren 2018
	Usage:
	scaffold_to_chr.py <tab-separated-coordinate-file> > <output-file-name>
		'''
	print(msg)

def main():
	try:
		infile = open(sys.argv[1], "r")
		chr_seg = sys.argv[2]
	except IndexError:
		usage()
		sys.exit()

	if len(sys.argv) < 2 or infile in ('-h', '--help'):
		usage()
		sys.exit()

	scaffold_order = ["superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

	# superscaffold54 length pre LM: 29256470
	Z_scaffold_length = {'superscaffold26': 25310599, 'superscaffold54':16379243 , 'superscaffold35': 4625539, 'superscaffold36': 9394175, 'superscaffold62': 2917291, 'superscaffold67': 5300260, 'superscaffold69-1': 5978518, 'superscaffold93': 4983591, 'superscaffold63': 1692925, 'superscaffold88': 624114, 'superscaffold83': 782506, 'superscaffold92': 2882843}
	PAR_scaffold_length = {'superscaffold26': 25310599, 'superscaffold54':16379243 , 'superscaffold35': 4625539, 'superscaffold36': 5869911}
	nonPAR_scaffold_length = {'superscaffold36': 3516672, 'superscaffold62': 2917291, 'superscaffold67': 5300260, 'superscaffold69-1': 5978518, 'superscaffold93': 4983591, 'superscaffold63': 1692925, 'superscaffold88': 624114, 'superscaffold83': 782506, 'superscaffold92': 2882843}

	PAR_start = [9394175, 3524264]
	nonPAR_start = [3516672, 0]
	# Read file into a dictionary
	vcf_dict = {}
	next(infile) # Remove the header from dictionary
	for line in infile:
		line = line.strip("\n").split()
		key, value = line[0], line[1:]
		if key in vcf_dict.keys():
			vcf_dict[key].append(value)
		else:
			vcf_dict[key] = [value]

	# Change the coordinates
	if chr_seg == "Z":
		for scaffold in scaffold_order:
			if scaffold == "superscaffold26":
				for l in vcf_dict[scaffold]:
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + "\t".join(l[0:2]) + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold54":
				new_list = []
				scaf_list = []
				for l in vcf_dict[scaffold]:
					if int(l[0]) <= Z_scaffold_length[scaffold]:
						scaf_list.append(l[0:])
						l[0] = str(Z_scaffold_length["superscaffold26"] + (Z_scaffold_length[scaffold] - int(l[0]) + 1))
						l[1] = str(Z_scaffold_length["superscaffold26"] + (Z_scaffold_length[scaffold] - int(l[1]) + 1))
						new_list.append(l)
				for i in reversed(new_list):
					ind = new_list.index(i)
					print ("ChrZ" + "\t" + str(int(i[1])-1) + "\t" + str(int(i[0])+1) + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind])
		# ******************************************************************************			
			elif scaffold == "superscaffold35":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold36":
				new_list = []
				scaf_list = []
				for l in vcf_dict[scaffold]:
					scaf_list.append(l[0:])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[0]) + 1))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[1]) + 1))
					new_list.append(l)
				for i in reversed(new_list):
					ind = new_list.index(i)
					print ("ChrZ" + "\t" + str(int(i[1])-1) + "\t" + str(int(i[0])+1) + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind])
		# ******************************************************************************			
			elif scaffold == "superscaffold62":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold67":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold69-1":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold93":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold63":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold88":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold83":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + Z_scaffold_length["superscaffold88"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + Z_scaffold_length["superscaffold88"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold92":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + Z_scaffold_length["superscaffold88"] + Z_scaffold_length["superscaffold83"] + int(l[0]))
					l[1] = str(Z_scaffold_length["superscaffold26"] + Z_scaffold_length["superscaffold54"] + Z_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + Z_scaffold_length["superscaffold62"] + Z_scaffold_length["superscaffold67"] + Z_scaffold_length["superscaffold69-1"] + Z_scaffold_length["superscaffold93"] + Z_scaffold_length["superscaffold63"] + Z_scaffold_length["superscaffold88"] + Z_scaffold_length["superscaffold83"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************				
	elif chr_seg == "PAR":
		for scaffold in scaffold_order:
			if scaffold == "superscaffold26":
				for l in vcf_dict[scaffold]:
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + "\t".join(l[0:2]) + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold54":
				new_list = []
				scaf_list = []
				for l in vcf_dict[scaffold]:
					if int(l[0]) <= PAR_scaffold_length[scaffold]: 
						scaf_list.append(l[0:2])
						l[0] = str(PAR_scaffold_length["superscaffold26"] + (PAR_scaffold_length[scaffold] - int(l[0]) + 1))
						l[1] = str(PAR_scaffold_length["superscaffold26"] + (PAR_scaffold_length[scaffold] - int(l[1]) + 1))
						new_list.append(l)
				for i in reversed(new_list):
					ind = new_list.index(i)
					print ("ChrZ" + "\t" + str(int(i[1])-1) + "\t" + str(int(i[0])+1) + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind])
		# ******************************************************************************			
			elif scaffold == "superscaffold35":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold36":
				new_list = []
				scaf_list = []
				for l in vcf_dict[scaffold]:
					if int(l[0]) >= PAR_start[1]:
						scaf_list.append(l[0:2])
						l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[0]) + 1))
						l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[1]) + 1))
						new_list.append(l)
				for i in reversed(new_list):
					ind = new_list.index(i)
					print ("ChrZ" + "\t" + str(int(i[1])-1) + "\t" + str(int(i[0])+1) + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind])		
		# ******************************************************************************
	elif chr_seg == "nonPAR":
		for scaffold in scaffold_order:
			if scaffold == "superscaffold36":
				new_list = []
				scaf_list = []
				for l in vcf_dict[scaffold]:
					if int(l[0]) <= nonPAR_start[0]:
						scaf_list.append(l[0:2])
						l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[0]) + 1))
						l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + (Z_scaffold_length[scaffold] - int(l[1]) + 1))
						new_list.append(l)
				for i in reversed(new_list):
					ind = new_list.index(i)
					print ("ChrZ" + "\t" + str(int(i[1])-1) + "\t" + str(int(i[0])+1) + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind])
		# ******************************************************************************			
			elif scaffold == "superscaffold62":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold67":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold69-1":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold93":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold63":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold88":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + int(l[1]))
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold83":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + nonPAR_scaffold_length["superscaffold88"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + nonPAR_scaffold_length["superscaffold88"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************			
			elif scaffold == "superscaffold92":
				scaf_list = []
				for l in vcf_dict[scaffold]:
					ind = vcf_dict[scaffold].index(l)
					scaf_list.append(l[0:2])
					l[0] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + nonPAR_scaffold_length["superscaffold88"] + nonPAR_scaffold_length["superscaffold83"] + int(l[0]))
					l[1] = str(PAR_scaffold_length["superscaffold26"] + PAR_scaffold_length["superscaffold54"] + PAR_scaffold_length["superscaffold35"] + Z_scaffold_length["superscaffold36"] + nonPAR_scaffold_length["superscaffold62"] + nonPAR_scaffold_length["superscaffold67"] + nonPAR_scaffold_length["superscaffold69-1"] + nonPAR_scaffold_length["superscaffold93"] + nonPAR_scaffold_length["superscaffold63"] + nonPAR_scaffold_length["superscaffold88"] + nonPAR_scaffold_length["superscaffold83"] + int(l[1]))			
					print ("ChrZ" + "\t" + l[0] + "\t" + l[1] + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l[2:]))
		# ******************************************************************************

	return 0

if __name__ == "__main__":
	main()
