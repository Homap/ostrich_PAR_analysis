#!/usr/bin/python
import sys

# Homa Papoli Yazdi
# It takes a file containing scaffold name and their length,
# with window size and step and outputs windows across the region in bed format.
# start and end of each interval is closed and open as such [start, end).

seq_l = open(sys.argv[1], 'r')
winSize = int(sys.argv[2])
step = winSize


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
	numOfChunks = ((int(sequence_l-winSize)/step))+1
	numOfChunks = int(numOfChunks) 

	for i in range(start, numOfChunks*step, step):
		yield i,i+winSize

print("CHROM"+"\t"+"Start"+"\t"+"End")
for line in seq_l:
    if not line.startswith("chrom"):
        start = int(line.strip("\n").split("\t")[1])
        sequence_l = int(line.strip("\n").split("\t")[2])
        if winSize > sequence_l:
            print(line.strip("\n").split("\t")[0]+"\t"+"0"+"\t"+str(sequence_l))
        else:
            for index, item in enumerate(slidingWindow(start, sequence_l, winSize, step)):
                if item[1] > sequence_l:
                    pass
                elif item[1] < sequence_l and item[1] + winSize > sequence_l:
                    print(line.strip("\n").split("\t")[0]+"\t"+str(item[0])+"\t"+str(item[1]))
                    print(line.strip("\n").split("\t")[0]+"\t"+str(item[1])+"\t"+str(sequence_l))
                else:
                    print(line.strip("\n").split("\t")[0]+"\t"+str(item[0])+"\t"+str(item[1]))
