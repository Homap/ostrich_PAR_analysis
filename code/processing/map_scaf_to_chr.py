#!/usr/bin/python
import sys
import pandas as pd
import matplotlib.pyplot as plt

lastz_chr = open("../data/lastz/chicken_NC_chr.txt")
lastz = open("../data/lastz/chicken_ostrich.lastz")
# pi = open("../data/pi/black_aut_100Kb.windowed.pi")
# pi_par = open("../data/pi/black_par_100Kb.windowed.pi")
td = open("black_A.100Kb.allsites.new.td.Tajima.D")
# td_par = open("../data/TajimaD/black_par.100Kb.allsites.Tajima.D")
# fst = open("../data/FST/black_male_female_aut_100Kb.windowed.weir.fst")
# fst_par = open("../data/FST/black_male_female_par_100Kb.windowed.weir.fst")


lastz_chr_dict = {}
for line in lastz_chr:
	line = line.rstrip().split()
	lastz_chr_dict[line[0]] = line[1]

lastz_dict = {}
for line in lastz:
	if not line.startswith("#name1"):
		line = line.rstrip().split()
		gg_scaf = line[0]
		sc_scaf = line[4]
		# if sc_scaf in lastz_dict:
		# 	lastz_dict[sc_scaf].append(gg_scaf)
		# else:
		lastz_dict[sc_scaf] = gg_scaf



# print(lastz_chr_dict)
# print(lastz_dict)

# outfile = open("pi.txt", "w")
# header = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI", "gg_chr"]
# outfile.write("\t".join(header)+"\n")
# for line in pi:
# 	if not line.startswith("CHROM"):
# 		line = line.rstrip().split()
# 		if int(line[3]) > 10:
# 			scaf = line[0]
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], line[4], lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile.write("\t".join(l)+"\n")

outfile_td = open("td.txt", "w")
header = ["CHROM", "BIN_START", "N_SNPS", "TajimaD", "gg_chr"]
outfile_td.write("\t".join(header)+"\n")
for line in td:
	if not line.startswith("CHROM"):
		line = line.rstrip().split()
		scaf = line[0]
		if scaf in lastz_dict:
			l = [scaf, line[1], line[2], line[3], lastz_chr_dict[lastz_dict[scaf]]]
			outfile_td.write("\t".join(l)+"\n")			

# outfile_fst = open("fst.txt", "w")
# header = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "gg_chr"]
# outfile_fst.write("\t".join(header)+"\n")
# for line in fst:
# 	if not line.startswith("CHROM"):
# 		line = line.rstrip().split()
# 		scaf = line[0]
# 		if float(line[4]) < 0:
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], "0", lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile_fst.write("\t".join(l)+"\n")
# 		else:
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], line[4], lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile_fst.write("\t".join(l)+"\n")	

# outfile_par = open("pi_par.txt", "w")
# header = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI", "gg_chr"]
# outfile_par.write("\t".join(header)+"\n")
# for line in pi_par:
# 	if not line.startswith("CHROM"):
# 		line = line.rstrip().split()
# 		if int(line[3]) > 10:
# 			scaf = line[0]
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], line[4], lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile_par.write("\t".join(l)+"\n")	

# outfile_td_par = open("td_par.txt", "w")
# header = ["CHROM", "BIN_START", "N_SNPS", "TajimaD", "gg_chr"]
# outfile_td_par.write("\t".join(header)+"\n")
# for line in td_par:
# 	if not line.startswith("CHROM"):
# 		line = line.rstrip().split()
# 		scaf = line[0]
# 		if scaf in lastz_dict:
# 			l = [scaf, line[1], line[2], line[3], lastz_chr_dict[lastz_dict[scaf]]]
# 			outfile_td_par.write("\t".join(l)+"\n")		

# outfile_fst_par = open("fst_par.txt", "w")
# header = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "gg_chr"]
# outfile_fst_par.write("\t".join(header)+"\n")
# for line in fst_par:
# 	if not line.startswith("CHROM"):
# 		line = line.rstrip().split()
# 		scaf = line[0]
# 		if float(line[4]) < 0:
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], "0", lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile_fst_par.write("\t".join(l)+"\n")
# 		else:
# 			if scaf in lastz_dict:
# 				l = [scaf, line[1], line[2], line[3], line[4], lastz_chr_dict[lastz_dict[scaf]]]
# 				outfile_fst_par.write("\t".join(l)+"\n")		

# pi_d = pd.read_table("pi.txt", sep='\s+')
# # print(pi_d)

# pi_d_par = pd.read_table("pi_par.txt", sep='\s+')



# chr1 = pi_d[pi_d["gg_chr"] == "chromosome_1"]
# chr2 = pi_d[pi_d["gg_chr"] == "chromosome_2"]
# chr3 = pi_d[pi_d["gg_chr"] == "chromosome_3"]
# chr4 = pi_d[pi_d["gg_chr"] == "chromosome_4"]
# chr5 = pi_d[pi_d["gg_chr"] == "chromosome_5"]
# chr6 = pi_d[pi_d["gg_chr"] == "chromosome_6"]
# chr7 = pi_d[pi_d["gg_chr"] == "chromosome_7"]
# chr8 = pi_d[pi_d["gg_chr"] == "chromosome_8"]
# chr9 = pi_d[pi_d["gg_chr"] == "chromosome_9"]
# chr10 = pi_d[pi_d["gg_chr"] == "chromosome_10"]
# chr11 = pi_d[pi_d["gg_chr"] == "chromosome_11"]
# chr12 = pi_d[pi_d["gg_chr"] == "chromosome_12"]
# chr13 = pi_d[pi_d["gg_chr"] == "chromosome_13"]
# chr14 = pi_d[pi_d["gg_chr"] == "chromosome_14"]
# chr15 = pi_d[pi_d["gg_chr"] == "chromosome_15"]
# chr17 = pi_d[pi_d["gg_chr"] == "chromosome_17"]
# chr18 = pi_d[pi_d["gg_chr"] == "chromosome_18"]
# chr19 = pi_d[pi_d["gg_chr"] == "chromosome_19"]
# chr20 = pi_d[pi_d["gg_chr"] == "chromosome_20"]
# chr21 = pi_d[pi_d["gg_chr"] == "chromosome_21"]
# chr22 = pi_d[pi_d["gg_chr"] == "chromosome_22"]
# chr23 = pi_d[pi_d["gg_chr"] == "chromosome_23"]
# chr24 = pi_d[pi_d["gg_chr"] == "chromosome_24"]
# chr25 = pi_d[pi_d["gg_chr"] == "chromosome_25"]
# chr26 = pi_d[pi_d["gg_chr"] == "chromosome_26"]
# chr27 = pi_d[pi_d["gg_chr"] == "chromosome_27"]
# #chrZ = pi_d[pi_d["gg_chr"] == "chromosome_Z"]
# chrZ = pi_d_par[pi_d_par["gg_chr"] == "chromosome_Z"]

# chrs = [
# chr1["PI"],
# chr2["PI"],
# chr3["PI"],
# chr4["PI"],
# chr5["PI"],
# chr6["PI"],
# chr7["PI"],
# chr8["PI"],
# chr9["PI"],
# chr10["PI"],
# chr11["PI"],
# chr12["PI"],
# chr13["PI"],
# chr14["PI"],
# chr15["PI"],
# chr17["PI"],
# chr18["PI"],
# chr19["PI"],
# chr20["PI"],
# chr21["PI"],
# chr22["PI"],
# chr23["PI"],
# chr24["PI"],
# chr25["PI"],
# chr26["PI"],
# chr27["PI"],
# chrZ["PI"]]

chrs_names = [
"chr1",
"chr2",
"chr3",
"chr4",
"chr5",
"chr6",
"chr7",
"chr8",
"chr9",
"chr10",
"chr11",
"chr12",
"chr13",
"chr14",
"chr15",
"chr17",
"chr18",
"chr19",
"chr20",
"chr21",
"chr22",
"chr23",
"chr24",
"chr25",
"chr26",
"chr27"]
# "PAR"]

# # fig7, ax7 = plt.subplots()
# fig = plt.figure(figsize=(10,5))
# # fig.set_title('Chromosome diversity')
# plt.boxplot(chrs)
# plt.xticks([i for i in range(1, 28)], chrs_names, rotation=30)
# plt.savefig("chr_diversity.png") 

td_d = pd.read_table("td.txt", sep='\s+')
td_d.dropna(subset = ["TajimaD"], inplace=True)
# td_par_d = pd.read_table("td_par.txt", sep='\s+')

chr1 = td_d[td_d["gg_chr"] == "chromosome_1"]
chr2 = td_d[td_d["gg_chr"] == "chromosome_2"]
chr3 = td_d[td_d["gg_chr"] == "chromosome_3"]
chr4 = td_d[td_d["gg_chr"] == "chromosome_4"]
chr5 = td_d[td_d["gg_chr"] == "chromosome_5"]
chr6 = td_d[td_d["gg_chr"] == "chromosome_6"]
chr7 = td_d[td_d["gg_chr"] == "chromosome_7"]
chr8 = td_d[td_d["gg_chr"] == "chromosome_8"]
chr9 = td_d[td_d["gg_chr"] == "chromosome_9"]
chr10 = td_d[td_d["gg_chr"] == "chromosome_10"]
chr11 = td_d[td_d["gg_chr"] == "chromosome_11"]
chr12 = td_d[td_d["gg_chr"] == "chromosome_12"]
chr13 = td_d[td_d["gg_chr"] == "chromosome_13"]
chr14 = td_d[td_d["gg_chr"] == "chromosome_14"]
chr15 = td_d[td_d["gg_chr"] == "chromosome_15"]
chr17 = td_d[td_d["gg_chr"] == "chromosome_17"]
chr18 = td_d[td_d["gg_chr"] == "chromosome_18"]
chr19 = td_d[td_d["gg_chr"] == "chromosome_19"]
chr20 = td_d[td_d["gg_chr"] == "chromosome_20"]
chr21 = td_d[td_d["gg_chr"] == "chromosome_21"]
chr22 = td_d[td_d["gg_chr"] == "chromosome_22"]
chr23 = td_d[td_d["gg_chr"] == "chromosome_23"]
chr24 = td_d[td_d["gg_chr"] == "chromosome_24"]
chr25 = td_d[td_d["gg_chr"] == "chromosome_25"]
chr26 = td_d[td_d["gg_chr"] == "chromosome_26"]
chr27 = td_d[td_d["gg_chr"] == "chromosome_27"]
#chrZ = pi_d[pi_d["gg_chr"] == "chromosome_Z"]
# chrZ = td_par_d[td_par_d["gg_chr"] == "chromosome_Z"]

chrs_td = [
chr1["TajimaD"],
chr2["TajimaD"],
chr3["TajimaD"],
chr4["TajimaD"],
chr5["TajimaD"],
chr6["TajimaD"],
chr7["TajimaD"],
chr8["TajimaD"],
chr9["TajimaD"],
chr10["TajimaD"],
chr11["TajimaD"],
chr12["TajimaD"],
chr13["TajimaD"],
chr14["TajimaD"],
chr15["TajimaD"],
chr17["TajimaD"],
chr18["TajimaD"],
chr19["TajimaD"],
chr20["TajimaD"],
chr21["TajimaD"],
chr22["TajimaD"],
chr23["TajimaD"],
chr24["TajimaD"],
chr25["TajimaD"],
chr26["TajimaD"],
chr27["TajimaD"]]
# chrZ["TajimaD"]]

plt.figure(figsize=(10,5))
# fig.set_title('Chromosome diversity')
plt.boxplot(chrs_td)
# plt.boxplot(chrZ["TajimaD"])
plt.xticks([i for i in range(1, 28)], chrs_names, rotation=30)
plt.savefig("chr_td_new.png") 

# fst_d = pd.read_table("fst.txt", sep='\s+')
# fst_d.dropna(subset = ["WEIGHTED_FST"], inplace=True)
# fst_par_d = pd.read_table("fst_par.txt", sep='\s+')

# chr1 = fst_d[fst_d["gg_chr"] == "chromosome_1"]
# chr2 = fst_d[fst_d["gg_chr"] == "chromosome_2"]
# chr3 = fst_d[fst_d["gg_chr"] == "chromosome_3"]
# chr4 = fst_d[fst_d["gg_chr"] == "chromosome_4"]
# chr5 = fst_d[fst_d["gg_chr"] == "chromosome_5"]
# chr6 = fst_d[fst_d["gg_chr"] == "chromosome_6"]
# chr7 = fst_d[fst_d["gg_chr"] == "chromosome_7"]
# chr8 = fst_d[fst_d["gg_chr"] == "chromosome_8"]
# chr9 = fst_d[fst_d["gg_chr"] == "chromosome_9"]
# chr10 = fst_d[fst_d["gg_chr"] == "chromosome_10"]
# chr11 = fst_d[fst_d["gg_chr"] == "chromosome_11"]
# chr12 = fst_d[fst_d["gg_chr"] == "chromosome_12"]
# chr13 = fst_d[fst_d["gg_chr"] == "chromosome_13"]
# chr14 = fst_d[fst_d["gg_chr"] == "chromosome_14"]
# chr15 = fst_d[fst_d["gg_chr"] == "chromosome_15"]
# chr17 = fst_d[fst_d["gg_chr"] == "chromosome_17"]
# chr18 = fst_d[fst_d["gg_chr"] == "chromosome_18"]
# chr19 = fst_d[fst_d["gg_chr"] == "chromosome_19"]
# chr20 = fst_d[fst_d["gg_chr"] == "chromosome_20"]
# chr21 = fst_d[fst_d["gg_chr"] == "chromosome_21"]
# chr22 = fst_d[fst_d["gg_chr"] == "chromosome_22"]
# chr23 = fst_d[fst_d["gg_chr"] == "chromosome_23"]
# chr24 = fst_d[fst_d["gg_chr"] == "chromosome_24"]
# chr25 = fst_d[fst_d["gg_chr"] == "chromosome_25"]
# chr26 = fst_d[fst_d["gg_chr"] == "chromosome_26"]
# chr27 = fst_d[fst_d["gg_chr"] == "chromosome_27"]
# #chrZ = pi_d[pi_d["gg_chr"] == "chromosome_Z"]
# chrZ = fst_par_d[fst_par_d["gg_chr"] == "chromosome_Z"]

# chrs_td = [
# ç,
# chr2["WEIGHTED_FST"],
# chr3["WEIGHTED_FST"],
# chr4["WEIGHTED_FST"],
# chr5["WEIGHTED_FST"],
# chr6["WEIGHTED_FST"],
# chr7["WEIGHTED_FST"],
# chr8["WEIGHTED_FST"],
# chr9["WEIGHTED_FST"],
# chr10["WEIGHTED_FST"],
# chr11["WEIGHTED_FST"],
# chr12["WEIGHTED_FST"],
# chr13["WEIGHTED_FST"],
# chr14["WEIGHTED_FST"],
# chr15["WEIGHTED_FST"],
# chr17["WEIGHTED_FST"],
# chr18["WEIGHTED_FST"],
# chr19["WEIGHTED_FST"],
# chr20["WEIGHTED_FST"],
# chr21["WEIGHTED_FST"],
# chr22["WEIGHTED_FST"],
# chr23["WEIGHTED_FST"],
# chr24["WEIGHTED_FST"],
# chr25["WEIGHTED_FST"],
# chr26["WEIGHTED_FST"],
# chr27["WEIGHTED_FST"],
# chrZ["WEIGHTED_FST"]]

# plt.figure(figsize=(10,5))
# # fig.set_title('Chromosome diversity')
# plt.boxplot(chrs_td)
# # plt.boxplot(chrZ["TajimaD"])
# plt.xticks([i for i in range(1, 28)], chrs_names, rotation=30)
# plt.savefig("chr_fst.png") 

# s54 = fst_par_d[fst_par_d["CHROM"] == "superscaffold54"]
# s26 = fst_par_d[fst_par_d["CHROM"] == "superscaffold26"]
# s35 = fst_par_d[fst_par_d["CHROM"] == "superscaffold35"]
# s36 = fst_par_d[fst_par_d["CHROM"] == "superscaffold36"]

# par_fst = [s26["WEIGHTED_FST"], s54["WEIGHTED_FST"], s35["WEIGHTED_FST"],s36["WEIGHTED_FST"]]
 
# par_names = ['s26', 's54', 's35', 's36']

# plt.figure(figsize=(10,5))
# # fig.set_title('Chromosome diversity')
# plt.boxplot(par_fst)
# # plt.boxplot(chrZ["TajimaD"])
# plt.xticks([i for i in range(1, 5)], par_names, rotation=30)
# plt.savefig("par_fst.png") 

# s54 = pi_d_par[pi_d_par["CHROM"] == "superscaffold54"]
# s26 = pi_d_par[pi_d_par["CHROM"] == "superscaffold26"]
# s35 = pi_d_par[pi_d_par["CHROM"] == "superscaffold35"]
# s36 = pi_d_par[pi_d_par["CHROM"] == "superscaffold36"]

# par_pi = [s26["PI"], s54["PI"], s35["PI"],s36["PI"]]
 
# par_names = ['s26', 's54', 's35', 's36']

# plt.figure(figsize=(10,5))
# # fig.set_title('Chromosome diversity')
# plt.boxplot(par_pi)
# # plt.boxplot(chrZ["TajimaD"])
# plt.xticks([i for i in range(1, 5)], par_names, rotation=30)
# plt.savefig("par_pi.png")

# s54 = td_par_d[td_par_d["CHROM"] == "superscaffold54"]
# s26 = td_par_d[td_par_d["CHROM"] == "superscaffold26"]
# s35 = td_par_d[td_par_d["CHROM"] == "superscaffold35"]
# s36 = td_par_d[td_par_d["CHROM"] == "superscaffold36"]

# par_td = [s26["TajimaD"], s54["TajimaD"], s35["TajimaD"],s36["TajimaD"]]
 
# par_names = ['s26', 's54', 's35', 's36']

# plt.figure(figsize=(10,5))
# # fig.set_title('Chromosome diversity')
# plt.boxplot(par_td)
# # plt.boxplot(chrZ["TajimaD"])
# plt.xticks([i for i in range(1, 5)], par_names, rotation=30)
# plt.savefig("par_td.png")
