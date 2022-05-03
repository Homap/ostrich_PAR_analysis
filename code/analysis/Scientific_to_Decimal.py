#!/usr/bin/python
import sys

f1 = open(sys.argv[1], "r")

for line in f1:
	line = line.strip("\n").split("\t")
	line[1] = int(float(line[1]))
	line[2] = int(float(line[2]))
	print("\t".join(str(i) for i in line))