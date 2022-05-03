#!/usr/bin/python

import sys

annotation_file = open(sys.argv[1], "r")

# annotation_file = open("/proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone_annotated_pi.txt")

# annotation_file = open("/proj/snic2020-16-269/private/homap/ostrich_z/bin/annotate_old_2/house_sparrow_annotation.txt")

# annotation_file = open("/proj/snic2020-16-269/private/homap/ostrich_z/bin/annotate_old_2/genome_annotated.txt")

count_4fold = {}
count_0fold = {}

gene_dict = {}
for line in annotation_file:
    # print(line)
    line = line.rstrip().split()
    # if not line[0].startswith("chr"):
    scaf_gene = line[0]
    if scaf_gene in gene_dict:
        gene_dict[scaf_gene].append(line[3:])
    else:
        gene_dict[scaf_gene] = [line[3:]]


print(gene_dict)

gene_dict_annot = {}
for element in gene_dict.keys():
    fourfold = 0
    zerofold = 0
    c = 0
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

print(gene_dict_annot)

variant_dict_annot = {}
for element in gene_dict.keys():
    fourfold_snp = 0
    zerofold_snp = 0
    for position in gene_dict[element]:
        syn = position[0]
        nonsyn = position[1]
        if syn == "1.0" and position[7] != "0":
            fourfold_snp += float(position[7])
        else:
            pass
        if nonsyn == "1.0" and position[7] != "0":
            zerofold_snp += float(position[7])
        else:
            pass
    four_zero_snp = [fourfold_snp, zerofold_snp]
    # print(four_zero)
    if not element in variant_dict_annot:
        variant_dict_annot[element] = four_zero_snp
    else:
        raise Exception("Duplicated gene ID, something is not right!")

print(variant_dict_annot)

# print("chr"+"\t"+"gene"+"\t"+"ps"+"\t"+"pn"+"\t"+"pnps")
# for gene in gene_dict_annot:
#     if not 0 in variant_dict_annot[gene]:
#         ps = variant_dict_annot[gene][0]/gene_dict_annot[gene][0]
#         pn = variant_dict_annot[gene][1]/gene_dict_annot[gene][1]
#         # print(gene, gene_dict_annot[gene], variant_dict_annot[gene], ps, pn, pn/ps)
#         print(gene.split(":")[0]+"\t"+gene.split(":")[1].split("-")[0]+"\t"+str(ps)+"\t"+str(pn)+"\t"+str(pn/ps))


# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt \
# /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone.frq.count
# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone.frq.count  > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone_annotated_pi.txt
# python fourfold_zerofold_diversity.py > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crotone_pnps.txt
# grep -f italy_testis_overdominant_gene_list.txt crotone_pnps.txt > overdominant_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/corsica.frq.count 40 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/corsica_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/corsica_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/corsica_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crete.frq.count 18 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crete_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crete_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/crete_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/guglionesi.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/guglionesi_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/guglionesi_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/guglionesi_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/rimini.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/rimini_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/rimini_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/rimini_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/malta.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/malta_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/malta_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/malta_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/sicily.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/sicily_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/sicily_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/sicily_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/house.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/house_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/house_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/house_pnps.txt

# python count_total_num_syn_nonsyn.py annotate_old_2/house_sparrow_annotation.txt /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/spanish.frq.count 20 > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/spanish_annotated_pi.txt
# python fourfold_zerofold_diversity.py /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/spanish_annotated_pi.txt > /proj/snic2020-6-222/Projects/Pitaliae/working/Homa/projects/sparrow_popgen/data/allele_counts/spanish_pnps.txt

