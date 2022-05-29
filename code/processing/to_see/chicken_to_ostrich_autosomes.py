#!/usr/bin/python
import sys
import numpy as np
#-----------------------------------------------------------------
# Homa Papoli
# The script takes the lastz output between chicken and ostrich
# and the translation of chicken scaffolds to chromosome and 
# lists ostrich scaffolds that corresponds to chicken chromosomes.
#-----------------------------------------------------------------

# Open the files
# scaffold:chromosome in chicken
gg_scaf_chr = open(sys.argv[1], "r")
# chicken ostrich lastz
lastz = open(sys.argv[2], "r")
# ostrich scaffold length
scaf_length = open(sys.argv[3], "r")
# 

macro_chrs = ["chromosome_"+str(i) for i in range(1, 6)]

# Dictionary of chicken scaffold
chicken_dict = {}
for line in gg_scaf_chr:
    line = line.strip('\n').split()
    scaf, chromosome = line[0], line[1]
    if not scaf in chicken_dict.keys():
        chicken_dict[scaf] = chromosome
    else:
        raise Exception("Sorry! 1 gg chromosome must correspond to 1 scaffold only!")

# Dictionary of chicken and ostrich lastz coordinates
lastz_dict = {}
lastz_sc_key_gg_value = {}
for line in lastz:
    line = line.strip('\n').split()
    if not line[0].startswith("#"):
        gg_scaf, gg_start, gg_end = line[0], line[1], line[2]
        sc_scaf, sc_start, sc_end = line[4], line[5], line[6]
        key = gg_scaf+":"+gg_start+":"+gg_end
        value = sc_scaf+":"+sc_start+":"+sc_end
        if not sc_scaf in lastz_sc_key_gg_value.keys():
            lastz_sc_key_gg_value[sc_scaf] = [gg_scaf]
        else:
            lastz_sc_key_gg_value[sc_scaf].append(gg_scaf)
        if not key in lastz_dict.keys():
            lastz_dict[key] = value
        else:
            raise Exception("Sorry! Same scaffold and coordinate are aligning to two different regions!")

print(lastz_dict)

for key in lastz_sc_key_gg_value.keys():
    lastz_sc_key_gg_value[key] = list(set(lastz_sc_key_gg_value[key]))
print(lastz_sc_key_gg_value)
# Dictionary of scaf and scaf length
sc_len = {}
for line in scaf_length:
    line = line.rstrip().split()
    if not line[0].startswith("C"):
        scaf, length = line[0], line[1]
        sc_len[scaf] = length

# For the comparison of the autosomes with the Z chromosomes, I am only interested in autosomes 1 to 10 since
# these are the macrochromosomes in aves with similar properties such as recombination rate to the Z chromosome.
# In particular, chromosome 4 and 5 are similar in size and hence similar in recombination rate to the Z chromosome.

header = ["gg_chr", "gg_start", "gg_end", "sc_scaffold", "sc_scaffold_length", "sc_start", "sc_end"]
print("\t".join(header))
sc_scaf_dict = {}
sc_scaf_gg_match = {}
for chrom in macro_chrs:
    for key in lastz_dict.keys():
        gg_scaf = key.split(":")[0]
        if chicken_dict.get(gg_scaf) == chrom:
            if not chrom in sc_scaf_dict.keys():
                sc_scaf_dict[chrom] = [lastz_dict[key].split(":")[0]]
            else:
                sc_scaf_dict[chrom].append(lastz_dict[key].split(":")[0])
            if not lastz_dict[key].split(":")[0] in sc_scaf_gg_match.keys():
                sc_scaf_gg_match[lastz_dict[key].split(":")[0]] = [abs(int(lastz_dict[key].split(":")[2])-int(lastz_dict[key].split(":")[1]))]
            else:
                sc_scaf_gg_match[lastz_dict[key].split(":")[0]].append(abs(int(lastz_dict[key].split(":")[2])-int(lastz_dict[key].split(":")[1])))
            if len(lastz_sc_key_gg_value[lastz_dict[key].split(":")[0]]) == 1: 
                output_list = [chicken_dict[gg_scaf], key.split(":")[1], key.split(":")[2], lastz_dict[key].split(":")[0], sc_len[lastz_dict[key].split(":")[0]], \
                    lastz_dict[key].split(":")[1], lastz_dict[key].split(":")[2]]
                print("\t".join(output_list))




# for element in sc_scaf_dict.keys():
#     sc_scaf_dict[element] = list(set(sc_scaf_dict[element]))
    
# print(sc_scaf_dict)
# print(sc_len)
# print(sc_scaf_gg_match)
# for element in sc_scaf_gg_match.keys():
#     output_l = [element, np.sum(sc_scaf_gg_match[element]), sc_len[element], round(np.sum(sc_scaf_gg_match[element])/int(sc_len[element]),2)*100]
#     print("\t".join([str(i) for i in output_l]))


# close the files
gg_scaf_chr.close()
lastz.close()
