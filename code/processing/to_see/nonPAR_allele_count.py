#!/use/bin/python
import sys

f = open(sys.argv[1], "r")
m = open(sys.argv[2], "r")

f_dict = {}
for line in f:
    if not line.startswith("CHROM"):
        line = line.rstrip().split()
        scaf = line[0]+":"+line[1]
        info = line[2:]
        f_dict[scaf] = info

m_dict = {}
for line in m:
    if not line.startswith("CHROM"):
        line = line.rstrip().split()
        scaf = line[0]+":"+line[1]
        info = line[2:]
        m_dict[scaf] = info

# print(f_dict)
# print(m_dict)
# superscaffold36 ['199866', '2', '10', 'A:0', 'C:10'] 0.0 5.0
print("CHROM"+"\t""POS"+"\t"+"N_ALLELES"+"\t"+"N_CHR"+"\t"+"{ALLELE:COUNT}")
for snp in f_dict.keys():
    a1 = int(f_dict[snp][2].split(":")[1])/2
    a2 = int(f_dict[snp][3].split(":")[1])/2
    a1_count_male_female = f_dict[snp][2].split(":")[0]+":"+str(a1 + int(m_dict[snp][2].split(":")[1]))
    a2_count_male_female = f_dict[snp][3].split(":")[0]+":"+str(a2 + int(m_dict[snp][3].split(":")[1]))
    line = [snp.split(":")[0], snp.split(":")[1], f_dict[snp][0], str(int(m_dict[snp][1])+a1+a2).split(".")[0], str(a1_count_male_female).split(".")[0], str(a2_count_male_female).split(".")[0]]
    print("\t".join(line))