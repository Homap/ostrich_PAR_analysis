#!/usr/bin/python
import sys
import matplotlib.pyplot as plt

# puts the recombination window file into correct format for R graph. 

window_rec = open(sys.argv[1], "r")

l = []
for line in window_rec:
	if not line.startswith("Scaffold"):
		line = line.strip("\n").split("\t")
		l.append(line)

#print l

for element in l:
	print(element[0]+"\t"+element[1]+"\t"+element[3])
	print(element[0]+"\t"+element[2]+"\t"+element[3])

