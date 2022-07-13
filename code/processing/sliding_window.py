#!/usr/bin/python
import sys
import numpy as np
# Homa Papoli Yazdi
# It takes a file containing scaffold name and their length,
# with window size and step and outputs windows across the region in bed format.
# start and end of each interval is closed and open as such [start, end).

seq_l = open(sys.argv[1], 'r')
winSize = int(sys.argv[2])
step = int(sys.argv[3])


# Function for sliding window
def slidingWindow(start, sequence_l, winSize, step):
	""" Returns a generator that will iterate through
	the defined chunks of input sequence. Input
	sequence must be iterable."""

	# Verify the inputs
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > sequence_l:
		pass
	#winSize = sequence_l
	#numOfChunks = ((int(sequence_l-winSize)/step))+1
	#numOfChunks = int(numOfChunks)
	#else:
	# Pre-compute number of chunks to emit
	window_start = []
	for i in range(start, (sequence_l - step), step):
		window_start.append(i)

	window_end = np.add(window_start, winSize)

	window_end[-1] = sequence_l

	window_list = []

	for c in range(0, len(window_start)):
		window_list.append([window_start[c], window_end[c]])

	return(window_list)



print("CHROM"+"\t"+"Start"+"\t"+"End")
for line in seq_l:
    if not line.startswith("chrom"):
        start = int(line.strip("\n").split("\t")[1])
        sequence_l = int(line.strip("\n").split("\t")[2])
        if winSize > sequence_l:
            print(line.strip("\n").split("\t")[0]+"\t"+"0"+"\t"+str(sequence_l))
        else:
            windows = slidingWindow(start, sequence_l, winSize, step)
            for index, item in enumerate(windows):
                print(line.strip("\n").split("\t")[0]+"\t"+str(item[0])+"\t"+str(item[1]))
