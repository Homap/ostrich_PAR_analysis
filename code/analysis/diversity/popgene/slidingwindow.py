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
    if sequence_l > numOfChunks*step:
        yield i+winSize, sequence_l
        #     for window in intervals:
#   File "/crex/proj/uppstore2017180/private/homap/ostrich_z/bin/popgene/slidingwindow.py", line 24, in slidingWindow
#     yield i+winSize, sequence_l
# UnboundLocalError: local variable 'i' referenced before assignment

# def test():
#         #check heterzygosity
#         if not slidingWindow(16, 2, 2) == 0.5:
#                 raise Exception
# test()
        