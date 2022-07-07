#!/usr/bin/python

import sys
import numpy as np

coverage = open(sys.argv[1], "r")
num_samples = int(sys.argv[2])
samples = open(sys.argv[3], "r")

samples_l = []
for line in samples:
    line = line.rstrip()
    samples_l.append(line)

lists = [[] for _ in range(num_samples)]

for line in coverage:
    if not line.startswith("#"):
        line = line.rstrip().split()
        for index, item in enumerate(lists):
            lists[index].append(int(line[index + 2]))

print("sample"+"\t"+"median_coverage")
for index, item in enumerate(lists):
    print(samples_l[index]+"\t"+str(np.median(item)))


