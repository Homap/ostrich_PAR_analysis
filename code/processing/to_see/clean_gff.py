#!/usr/bin/python
import sys

gff = open(sys.argv[1], "r")
gff_l = gff.readlines()

# We want to extract only mRNA with their corresponding 
# exon and CDS. 
#*****************************************************************************
# Read the mRNA and CDS coordinates into dictionaries
#*****************************************************************************
mRNA_coordinates = {}
for line in gff_l:
	if not line.startswith("#"):
		line = line.strip('\n').split('\t')
		if line[2] == 'mRNA':
			ID = line[8].split(";")[0].split("=")[1]
#			print(line)
#			ID = line[8].split(";")[1].split(",")[0].split(":")[1]
#			print(ID)
			key, value = line[0]+"_"+line[8].split(";")[0].split("=")[1], line[3:5] 
			if key in mRNA_coordinates.keys(): 
				mRNA_coordinates[key].append(value) 
			else:
				mRNA_coordinates[key] = value

#print(mRNA_coordinates)
# {'scaffold1490_rna-XM_009676741.1': ['351', '858'], 

# There are several mRNAs with different IDs but corresponding to exactly
# the same DNA coordinates. They are essentially duplicates and only one
# should be kept. 
rev_mRNA_dict = {}
for key, value in mRNA_coordinates.items():
	#print(key, value)
	scaf = key.split("_")[0]
	value = scaf+'_'+value[0]+'_'+value[1]
	rev_mRNA_dict.setdefault(value, set()).add(key)
#print(rev_mRNA_dict)

new_mRNA_dict = {} 
mRNA_IDs = []
for key, values in rev_mRNA_dict.items():
	if len(values) > 1:
		#print(key, values)
		#print(list(rev_mRNA_dict[key])[0])
		new_key = list(rev_mRNA_dict[key])[0]
		start = key.split("_")[1]
		end = key.split("_")[2]
		coord = [start, end]
		#print(coord)
		new_mRNA_dict[new_key] = coord
		#print(new_mRNA_dict)
		if not new_key in mRNA_IDs:
				mRNA_IDs.append(new_key.split("_")[1]+"_"+new_key.split("_")[2])
	else:
		new_mRNA_dict[list(values)[0]] = [key.split("_")[1], key.split("_")[2]]
		if not list(values)[0] in mRNA_IDs:
				mRNA_IDs.append(list(values)[0].split("_")[1]+"_"+list(values)[0].split("_")[2])
		#print(new_mRNA_dict)

#print(len(set(new_mRNA_dict.keys())))
#print(len(set(mRNA_IDs)))
#print(mRNA_IDs)

# Remove these plus take the longest transcript in overlapping mRNAs.
# Get the CDS fasta sequence of clean mRNA set.
# Translate them and then blast them to be sure they are correct.
# This will be the final set of CDS. 
# Use these CDS to calculate dNdS per gene and per window.
# To do per window, you simply need to know the coordinate of mRNA of that
# gene and count the number of dN and dS in the respective window.
# With the vcf file, you do the same per each gene and also per window for
# pN and pS.
# and there you'll have your pN, pS, dN and dS estimates.
# You can also annotate each coding site in the genome using your script and also
# snpEff with the clean gff as database. In this way, you can be double sure.
# Finish this this weekend and that would be a dream!

#print(mRNA_coordinates.keys())
#print(len(mRNA_coordinates.keys()))
#print(len(set(mRNA_coordinates.keys())))

#print(mRNA_IDs)
#print(len(set(mRNA_IDs))) # 24547
# superscaffold20_rna-XM_009689810.1 ['9963557', '9974409']
# superscaffold20_rna-XM_009689894.1 ['9963557', '9974409']
# superscaffold20_rna-XM_009689731.1 ['9963557', '9974409']

#for key in mRNA_coordinates:
# 	for coordinates in mRNA_coordinates[key]:
# 		print(key, coordinates)

# 'scaffold1490_rna-XM_009676741.1'

# scaffold1030_rna-XM_009670518.1 ['27932', '38640']
# scaffold1030_rna-XM_009672115.1 ['27932', '38640']
# scaffold1030_rna-XM_009672903.1 ['27932', '38640']
# scaffold1030_rna-XM_009673656.1 ['27932', '38640']
# scaffold1030_rna-XM_009671305.1 ['27932', '38640']
# CDS dictionary with CDS that belong to the mRNA, not any other feature.
# To do that, the Parent id of a CDS must mactch with the ID of mRNA.

CDS_coordinates = {}
for line in gff_l:
	if not line.startswith("#"):
		line = line.strip('\n').split('\t')
		#print(line)
		scaf = line[0]
		if line[2] == 'CDS':
			mRNA_id = line[8].split(";")[1].split("=")[1]
			#print(mRNA_id)
			if mRNA_id in mRNA_IDs:
				protein_id = line[8].split(";")[2].split(",")[1].split(":")[1]
				#print(protein_id)
				key = scaf+'_'+mRNA_id+'_'+protein_id
				value = [line[3], line[4], line[6], line[7]]
				if key in CDS_coordinates.keys():
					CDS_coordinates[key].append(value)
				else:
					CDS_coordinates[key] = [value]
#print(CDS_coordinates)

for key in CDS_coordinates:
	for coord in CDS_coordinates[key]:
		print(key+"\t"+"\t".join(coord))

#print(CDS_coordinates)

gff.close()
