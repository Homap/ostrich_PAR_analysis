#!/usr/bin/python
import sys


seq_l = open(sys.argv[1], 'r')
winSize = int(sys.argv[2])
step = winSize

# Function for sliding window
def slidingWindow(sequence_l, winSize, step):
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

	for i in range(0, numOfChunks*step, step):
		yield i,i+winSize

print("CHROM"+"\t"+"Start"+"\t"+"End")
for line in seq_l:
    if not line.startswith("chrom"):
        sequence_l = int(line.strip("\n").split("\t")[2])
        if winSize > sequence_l:
            print(line.strip("\n").split("\t")[0]+"\t"+"1"+"\t"+str(sequence_l))
        else:
            for index, item in enumerate(slidingWindow(sequence_l, winSize, step)):
                print(line.strip("\n").split("\t")[0]+"\t"+str(item[0]+1)+"\t"+str(item[1]))
                # For the last window when sequence_l%winSize != 0
                if (sequence_l-item[1]) < winSize and (sequence_l%winSize) != 0:
                    print(line.strip("\n").split("\t")[0]+"\t"+str(item[1]+1)+"\t"+str(sequence_l))