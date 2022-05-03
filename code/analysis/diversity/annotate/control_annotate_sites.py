#!/usr/bin/python

# annotation = open("house_sparrow_annotation.txt", "r")

syn = open("synonymoussites.txt", "r")
nonsyn = open("nonsynonymoussites.txt", "r")

print("populating annotation list")
annotation_list = {}
fourfold_list = []
zerofold_list = []

with open("house_sparrow_annotation.txt", "r") as annotation:
    next(annotation)
    for line in annotation:
        if not line.startswith("scaffold"):
            line = line.strip("\n").split()
            chr = line[0]
            pos = str(int(line[4])+1)
            info = [line[3], line[5], line[6], line[7]]
            element = chr + ":" + pos
            # annotation_list[element] = info
            if line[5] == "1.0":
                fourfold_list.append(element)
            else:
                pass
            if line[6] == "1.0":
                zerofold_list.append(element)
            # annotation_list.append(element)
print("DONE!")
#print(annotation_list[10000:12000])
print(len(fourfold_list))
print(len(zerofold_list))

print("syn overlap")
syn_over = open("syn_overlap.txt", "w")
for line in syn:
    line = line.strip("\n")
    if line in fourfold_list:
        syn_over.write(line + "\n")
        
print("nonsyn overlap")
nonsyn_over = open("nonsyn_overlap.txt", "w")
for line in nonsyn:
    line = line.strip("\n")
    if line in zerofold_list:
        nonsyn_over.write(line + "\n")

# /proj/snic2020-6-222/Projects/Pitaliae/working/Julio/Data/nonsynonymous_missing.txt
# /proj/snic2020-6-222/Projects/Pitaliae/working/Julio/Data/synonymous_missing.txt