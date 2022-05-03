import math
# **********************************
# Function for N choose r
# **********************************
def ncr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
# **********************************
# Function for Pairwise differences
# **********************************
def pi_fun(nchr, allele1, allele2):
        P_diff = allele1*allele2
        NCR = ncr(nchr, 2)
        pi = round(P_diff/NCR, 4)
        return (P_diff, NCR, pi) 