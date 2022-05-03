#!/usr/bin/python

import sys

annotation_file = open(sys.argv[1], "r")

count_4fold = {}
count_0fold = {}

gene_dict = {}
for line in annotation_file:
    line = line.rstrip().split()
    # if not line[0].startswith("chr"):
    scaf_gene = line[0]+":"+line[2]
    if scaf_gene in gene_dict:
        gene_dict[scaf_gene].append(line[3:])
    else:
        gene_dict[scaf_gene] = [line[3:]]


# print(gene_dict)

gene_dict_annot = {}
for element in gene_dict.keys():
    fourfold = 0
    zerofold = 0
    for position in gene_dict[element]:
        syn = position[0]
        nonsyn = position[1]
        if syn == "1.0":
            fourfold += 1
        else:
            pass
        if nonsyn == "1.0":
            zerofold += 1
        else:
            pass
    four_zero = [fourfold, zerofold]
    # print(four_zero)
    if not element in gene_dict_annot:
        gene_dict_annot[element] = four_zero
    else:
        raise Exception("Duplicated gene ID, something is not right!")

variant_dict_annot = {}
for element in gene_dict.keys():
    fourfold_snp = 0
    zerofold_snp = 0
    fourfold_snp_count = 0
    zerofold_snp_count = 0
    for position in gene_dict[element]:
        syn = position[0]
        nonsyn = position[1]
        if syn == "1.0" and position[7] != "0":
            fourfold_snp += float(position[7])
            fourfold_snp_count += 1
        else:
            pass
        if nonsyn == "1.0" and position[7] != "0":
            zerofold_snp += float(position[7])
            zerofold_snp_count += 1
        else:
            pass
    four_zero_snp = [fourfold_snp, zerofold_snp, fourfold_snp_count, zerofold_snp_count]
    # print(four_zero)
    if not element in variant_dict_annot:
        variant_dict_annot[element] = four_zero_snp
    else:
        raise Exception("Duplicated gene ID, something is not right!")

#print(variant_dict_annot)

print("chr"+"\t"+"gene"+"\t"+"ps"+"\t"+"pn"+"\t"+"pnps"+"\t"+"NS"+"\t"+"NN"+"\t"+"S4fold"+"\t"+"S0fold"+"\t"+"whole_syn"+"\t"+"whole_nonsyn")
for gene in gene_dict_annot:
    if not 0 in variant_dict_annot[gene]:
        ps = variant_dict_annot[gene][0]/gene_dict_annot[gene][0]
        pn = variant_dict_annot[gene][1]/gene_dict_annot[gene][1]
        # print(gene, gene_dict_annot[gene], variant_dict_annot[gene], ps, pn, pn/ps)
        print(gene.split(":")[0]+"\t"+gene.split(":")[1].split("-")[0]+"\t"+str(ps)+"\t"+str(pn)+"\t"+str(pn/ps)+"\t"+str(variant_dict_annot[gene][0])+"\t"+str(variant_dict_annot[gene][1])+"\t"+str(variant_dict_annot[gene][2])+"\t"+str(variant_dict_annot[gene][3])+"\t"+str(gene_dict_annot[gene][0])+"\t"+str(gene_dict_annot[gene][1]))
