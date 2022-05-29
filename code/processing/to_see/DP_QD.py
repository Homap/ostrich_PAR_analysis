#!/usr/bin/python
import sys
import re
import gzip
import matplotlib.pyplot as plt
import numpy as np

dp_array = np.empty(0, dtype=float)
qd_array = np.empty(0, dtype=float)

c = 0
with gzip.open(sys.argv[1], 'rb') as vcf:
	for line in vcf:
		line = line.decode() # decode byte into string
		if not line.startswith("#"):
			line = line.strip("\n").split("\t")
			c += 1
			#print(c)
			#print(line)
			info = line[7]
			# print(info)
			pattern1 = re.compile(r'DP=([^;]+)')
			depth = re.search(pattern1, info).group(1)
			dp_array = np.append(dp_array, depth)
			#print(dp_array)
			pattern2 = re.compile(r'QD=([^;]+)')
			qd = re.search(pattern2, info).group(1)
			qd_array = np.append(qd_array, qd)
			#print(qd_array)

# Calculate the average
depth_mean = np.array(dp_array)
qd_mean = np.array(qd_array)

print("mean", depth_mean)
print("median", qd_mean)

depth_median = np.array(dp_array)
qd_median = np.array(qd_array)

print("mean", depth_median)
print("median", qd_median)

# Calculate the min and max
plt.hist(depth_mean, bins = 20, facecolor='blue', alpha=0.5)
plt.hist(qd_mean, bins = 20, facecolor='blue', alpha=0.5)
plt.savefig("dp.qd.pdf")

