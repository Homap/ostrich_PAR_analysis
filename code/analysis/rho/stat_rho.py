#!/usr/bin/python

#----------------------------------------------------------------
# Author: Homa Papoli - June 2022
# Takes LDhat stat output and position files of SNPs for each window
# and outputs two files, one map length for each window in a scaffold
# and the other per site rho.

# Run: 
# python summarise_rho.py ldhatstat.res.txt pos.txt scaffold.per.site scaffold.map.length
#----------------------------------------------------------------


import sys

def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)

def append_multiple_lines(file_name, lines_to_append):
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        appendEOL = False
        # Move read cursor to the start of file.
        file_object.seek(0)
        # Check if file is not empty
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Iterate over each string in the list
        for line in lines_to_append:
            file_object.write(line+"\n")


# read the stat output
stat_output = open(sys.argv[1], "r")
pos = open(sys.argv[2], "r")
per_site_rho = sys.argv[3]
length_file = sys.argv[4]
scaffold = sys.argv[5]


# read stat and position files into list
stat_list = stat_output.readlines()
pos_list = pos.readlines()[1:]


list_out = []
for line in stat_list[2:]:
    line = line.rstrip().split("\t")
    loci = int(line[0].strip())
    mean = line[1].strip()
    median = line[2].strip()
    L95 = line[3].strip()
    U95 = line[4].strip()

    locus_start = pos_list[loci-1].rstrip()
    locus_end = pos_list[loci].rstrip()

    out = str(scaffold+"\t"+locus_start)+"\t"+str(locus_end)+"\t"+mean+"\t"+median+"\t"+L95+"\t"+U95
    list_out.append(out)

append_multiple_lines(per_site_rho, list_out)

region_rho = scaffold+"\t"+pos_list[0].rstrip()+"\t"+pos_list[-1].rstrip()+"\t"+stat_list[1].rstrip().split("\t")[1]+"\t"\
    +stat_list[1].rstrip().split("\t")[2]+"\t"+stat_list[1].rstrip().split("\t")[3]+"\t"+stat_list[1].rstrip().split("\t")[4]
append_new_line(length_file, region_rho)




