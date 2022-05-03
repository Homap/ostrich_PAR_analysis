import numpy as np
import math
# **********************************
# Function for Tajima's D # Can only calculate per window
#' R script to calculate Tajima's D based on site frequency spectrum
#' Reference:
#' Tajima F., 1983 Evolutionary relationship of DNA sequences in finite
#'   populations. Genetics 123: 437–460.
#' Tajima F., 1989 Statistical method for testing the neutral mutation 
#'   hypothesis by DNA polymorphism. Genetics 123: 585–595.
#' Justin C. Fay and Chung-I Wu, 2000 Hitchhiking Under Positive Darwinian
#'   Selection. Genetics 155: 1405-1413``
# **********************************

def TajimaD(sfs, num_chr, theta_pi): # folded SFS
    #' sfs (site frequency spectrum): number of singletons, doubletons, ..., etc
    # n = len(sfs) + 1
    n = 2*len(sfs) #+1 # length of SFS is not always equal to the number of chromosomes because if in
    # a given window, we have few SNPs we might have only few categories, for example, only SNPs of a certain frequency and not all possible frequencies.
    # In such cases, the calculation of Tajima's D with this equaltion is not correct and is inflated.
    ss = sum(sfs)
    
    a1 = sum([1/i for i in range(1,num_chr)]) #sum(1 / seq_len(n-1)) # seq_len: Harmonic mean
    a2 = sum([1/i**2 for i in range(1,num_chr)])
    
    b1 = (num_chr + 1) / (3 * (num_chr - 1))
    b2 = 2 * (num_chr**2 + num_chr + 3)/(9 * num_chr * (num_chr - 1))
    
    c1 = b1 - 1/a1
    c2 = b2 - (num_chr + 2)/(a1 * num_chr) + a2 / a1**2
    
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    
    Vd = e1 * ss + e2 * ss * (ss - 1) 
    
    #theta_pi = sum(2 * seq_len(n-1) * (n - seq_len(n-1)) * sfs)/(n*(n-1))
    # first_multiply = np.multiply([2*i for i in range(1,int(n/2)+1)], [n-i for i in range(1,int(n/2)+1)])
    # theta_pi = sum(np.multiply(first_multiply, sfs))/(n*(n-1))
    theta_w = ss / a1
    res = (theta_pi - theta_w) / math.sqrt(Vd)
    return(res, theta_pi, theta_w)