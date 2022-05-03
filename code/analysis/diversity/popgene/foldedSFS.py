import numpy as np
from collections import Counter
# **********************************
# Function for folded SFS
# **********************************
# Scaffold1 pos allele1_count allele2_count
# f = [[1, 19], [1, 19], [1, 19], [2, 18], [16, 4]]
def folded_sfs(geno_array):
    minor_array = np.min(geno_array, axis=1)
    sfs_dict = Counter(minor_array)
    return(sfs_dict)