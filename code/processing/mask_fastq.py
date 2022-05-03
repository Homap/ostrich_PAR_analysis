#!/use/bin/python
import sys
import time
start_time = time.time()


#fastq = open("test.fq")
#../data/bed/SC.allRepeats.bed
repeat_coord = {}
with open(sys.argv[1]) as repeats:
	for line in repeats:
		line = line.strip().split()
		if line[0] in repeat_coord.keys():
			repeat_coord[line[0]].append([int(line[1]), int(line[2])])
		else:
			repeat_coord[line[0]] = [[int(line[1]), int(line[2])]]

#print(repeat_coord)
#coordinates = [[3, 6], [10, 15]]

#@scaffold53
#nnnntgatcgacagggcctctctctctctgtctaacacgcagagttgatcggc
#+
#!!!!?????????????BBBBBBBBBBBBEEEEEEEEEEEEEEEEEEEEEEEE
# out.fq
# Print every second line.
step = 4
with open(sys.argv[2]) as handle:
	for lineno, line in enumerate(handle):
		if lineno % step == 0:
			header = line.strip().replace("@", "")
			print("@"+header)
		elif lineno % step == 1:
			result = list(line.strip())
			if header in repeat_coord.keys():
				for start, end in repeat_coord[header]:
					result[start:end] = ["N"]*(end-start)
					#length = element[1] - element[0]
					#line = line.replace(line[element[0]:element[1]], 'N'*length)
				#print(line)
			res = ''.join(result)
			#print(line.strip())
			print(res)
		elif lineno % step == 2:
			print(line.strip())
		elif lineno % step == 3:
			qual = list(line.strip())
			if header in repeat_coord.keys():
				for start,end in repeat_coord[header]:
					qual[start:end] = ["!"]*(end-start)
			quality = ''.join(qual)
			#print(line.strip())
			print(quality)

#print("--- %s seconds ---" % (time.time() - start_time))