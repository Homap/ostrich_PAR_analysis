#!/usr/bin/python
from __future__ import division
import sys

#****************************************************************************
# Written by Homa Papoli - Jan 2019
#****************************************************************************
# Script calculates recombination rate per window in cM per base pair (cM/bp)
# Logic:
# Given recombination rate per site between Marker 1 and Marker 2:
# _________|_______________|__________
#       Marker 1        Marker 2
# Case 1:
# _________|_______________|__________
#       Marker 1        Marker 2
#            |__________|
#	 Window_start    Window_end
# Real example:
# 'superscaffold62:1000000:2000000': [['superscaffold62', '207526', '2434793', '0.0', '55917482', '58144749', '1000000']]
# Calculation: 1000000*0.0/1000000
# Case 2:
#				r1		  r2
# _________|_________|_________|_
#       Marker 1   Marker 2  Marker 3
#            |__________|
#	 Window_start    Window_end
# Calculation: r1*(Overlap of window with the distance between Marker 1 and Marker 2) + r2*(Overlap of window with the distance between Marker 2 and Marker 3)
# Real example:
# 'superscaffold26:3000000:4000000': [['superscaffold26', '1113307', '3461665', '2.19131835947e-06', '1113307', '3461665', '461665'], 
#                                     ['superscaffold26', '3461665', '4790876', '4.25891750821e-06', '3461665', '4790876', '538335']]
# Calculation: (461665*2.19131835947e-06) + (538335*4.25891750821e-06) = 3.304379347206948
# Recombination rate: 3.304379347206948/(461665+538335) = 3.304379347206948e-06
#*****************************************************************************

f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

# Read the window file into a list as such: ['ChrZ:1:1000000', 'ChrZ:1000001:2000000', 'ChrZ:2000001:3000000', 'ChrZ:3000001:4000000'...]
f2_l = []
for line in f2:
	if not line.startswith("CHROM"):
		if not line.startswith("Scaffold"):
			line = line.strip("\n").split("\t")
			f2_l.append(line[0]+":"+line[1]+":"+line[2])

# print(f2_l)

# Next, read the window-popmeasure overlap file into a dictionary: {'ChrZ:1:1000000': [['ChrZ', '1306', '190500', 'superscaffold26', '1306', '190500', '199.28876', '199...]...}
# {'superscaffold26:0:1000000': [['superscaffold26', '1306', '201908', ' 219.53464', ' 219.42760', ' 209.30981', ' 230.18001', '200602']...}
f1_dict = {}
for line in f1:
	if not line.startswith("CHROM"):
		line = line.strip("\n").split("\t")
		key = line[0]+":"+line[1]+":"+line[2]
		value = line[3:]
		if key in f1_dict.keys():
			f1_dict[key].append(value)
		else:
			f1_dict[key] = [value]

# print(f1_dict)

print("window"+"\t"+"Window_start"+"\t"+"Window_end"+"\t"+"rho_per_bp"+"\t"+"rho_per_window")
for key in f2_l: # for each window in the order of the sliding window list
	if len(f1_dict[key]) == 1: # if there is only 1 occurrence of that window
		# print(key)
		if f1_dict[key][0][0] == ".": # if there is no overlap
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"NA"+"\t"+"NA")
		elif f1_dict[key][0][6] == "NA": 
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"NA"+"\t"+"NA")
		else: # if there is an overlap
			rho_per_site = float(f1_dict[key][0][6])
			window_size = int(key.split(":")[2]) - int(key.split(":")[1])
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+str(round(rho_per_site, 6))+"\t"+str(float(f1_dict[key][0][7])))
	else: # if there is more than one overlap
		# print(key)
		addup_pair_rho = 0
		total_pair = 0
		for element in f1_dict[key]:
			# print(element)
			if not element[6] == "NA":
				rho_per_site = float(element[6])
				addup_pair_rho = addup_pair_rho + (int(element[8])*rho_per_site)
				total_pair = total_pair + int(element[8])
		new_pair_rho = addup_pair_rho/total_pair
		window_size = int(key.split(":")[2]) - int(key.split(":")[1])
		#print(new_pair_rho*window_size)
		print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+str(round(new_pair_rho,6))+"\t"+str(round(new_pair_rho*window_size, 6)))






