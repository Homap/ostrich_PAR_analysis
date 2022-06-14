import msprime
import numpy as np
import pandas as pd
from functions_Conv import cMtoR
import matplotlib.pyplot as plt

# Import dataframe with pairwise rho, r, Ne, estimates
cumRhoFile = pd.read_csv('../data/rho_scaffold_36_cumRho_thin.csv') # Change the path to the one in the repository
print(cumRhoFile)
thin_pos =  np.array(cumRhoFile['thin.pos'])


# pull out relevant variables from ss36File dataframe
cumRho      = np.array(cumRhoFile['thin.cumRho'][0:145]) # the cumulative Rho as you move away from the PAR boundary, up to rho = 100
Ne_smooth_r = np.array(cumRhoFile['thin.Ne.smooth'][0:145]) # the Ne obtained using smooth_rate_map_cM_Mb
Ne_real_r   = np.array(cumRhoFile['thin.Ne'][0:145]) # the Ne obtained using real_rate_map_cM_Mb
 
# Define per generation neutral mutation rate 
# -- mu for Struthioformes from Zhang et al. 2014 Fig.3
# -- generation time for ostriches from Papoli et al. 2020
# mu = 1.0*10**(-3)*10**(-6)*16.8
#mu = 1.68e-08

# Define per generation neutral mutation rate 
# -- ratite-specific dS per site per year from Yazdi & Ellegren (2014)
# 1.1e-9
# -- generation time for ostriches from Papoli et al. 2020
# 16.8
mu = 1.1e-9*16.8


# Define number of sampled chromosomes
# [TZZ: #Z's, TZZ: #W's, TZW: #Z's, TZW: #W's, TWW: #Z's, TWW: #W's ]
BlackSamp  = np.array([15, 0, 15, 5, 0, 5])
#BluRedSamp = np.array([14, 0, 14, 4, 0, 4])
#simpleSamp = np.array([2, 0, 1, 1, 0, 2])


# Define simulation function which accepts above variables
def sexCoal_2deme_PARbound(chromSamp=BlackSamp, rho=1, Ne=1, num_replicates=10**4, mu=mu):
    # Effective population size of Z and W chromosomes
    Ne_Z = (3*Ne)/4
    Ne_W = Ne/4
    # r_f is the recombination rate in females
    r_f = rho/(4*Ne)
    # # of demes (d=2 because we are modeling two chromosome classes: Xs and Ys).
    d = 2
    # Allocate the initial samples among X and Y chromosomes.
    population_configurationsZZ = [
        msprime.PopulationConfiguration(sample_size=chromSamp[0], initial_size=Ne_Z),
        msprime.PopulationConfiguration(sample_size=chromSamp[1], initial_size=Ne_W)]
    population_configurationsZW = [
        msprime.PopulationConfiguration(sample_size=chromSamp[2], initial_size=Ne_Z),
        msprime.PopulationConfiguration(sample_size=chromSamp[3], initial_size=Ne_W)]
    population_configurationsWW = [
        msprime.PopulationConfiguration(sample_size=chromSamp[4], initial_size=Ne_Z),
        msprime.PopulationConfiguration(sample_size=chromSamp[5], initial_size=Ne_W)]
    # Now we set up the migration matrix.
    recombination_matrix =  [[0  , (1/3)*r_f],
                             [r_f,         0]]
    # We pass these values to the simulate function, and ask it
    # to run the required number of replicates.
    replicatesZZ = msprime.simulate(Ne=Ne,
        population_configurations=population_configurationsZZ,
        migration_matrix=recombination_matrix,
        num_replicates=num_replicates,
        mutation_rate=mu)
    replicatesZW = msprime.simulate(Ne=Ne,
        population_configurations=population_configurationsZW,
        migration_matrix=recombination_matrix,
        num_replicates=num_replicates,
        mutation_rate=mu)
    replicatesWW = msprime.simulate(Ne=Ne,
        population_configurations=population_configurationsWW,
        migration_matrix=recombination_matrix,
        num_replicates=num_replicates,
        mutation_rate=mu)
    # And then calculate the coalescence time for each replicate tree
    TZZ    = np.zeros(num_replicates)
    TZW    = np.zeros(num_replicates)
    TWW    = np.zeros(num_replicates)
    for i, tree_sequence in enumerate(replicatesZZ):
        tree = tree_sequence.first()
        TZZ[i] = TbarWin(chromSamp=chromSamp, tree=tree, chromType="Z") / ((2 * Ne))
    print('Finished T_ZZ')
    for i, tree_sequence in enumerate(replicatesZW):
        tree = tree_sequence.first()
        TZW[i] = TbarBtw(chromSamp=chromSamp, tree=tree) / (2 * Ne)
    print('Finished T_ZW')
    for i, tree_sequence in enumerate(replicatesWW):
        tree = tree_sequence.first()
        TWW[i] = TbarWin(chromSamp=chromSamp, tree=tree, chromType="W") / ((2 * Ne))
    print('Finished T_WW')
    # Calculate quantities of interest
    Tbar_ZZ    =  np.mean(TZZ)
    Tbar_ZZ_ci =  [np.percentile(TZZ, 2.5), np.percentile(TZZ, 97.5)]
    Tbar_ZW    =  np.mean(TZW)
    Tbar_ZW_ci =  [np.percentile(TZW, 2.5), np.percentile(TZW, 97.5)]
    Tbar_WW    =  np.mean(TWW)
    Tbar_WW_ci =  [np.percentile(TWW, 2.5), np.percentile(TWW, 97.5)]
    Tw         =  (TZZ + TWW)/2
    Tbar_w     =  np.mean(Tw)
    Tw_ci      =  [np.percentile(Tw, 2.5), np.percentile(Tw, 97.5)]
    Tt         =  (TZZ + 2*TZW + TWW)/4
    Tbar_t     =  np.mean(Tt)
    Tt_ci      =  [np.percentile(Tt, 2.5), np.percentile(Tt, 97.5)]
    FstZW      =  1 - (Tw/Tt)
    FstZW      =  FstZW[np.isfinite(FstZW)]
    FstBarZW   =  np.mean(FstZW)
    FstZW_ci   =  [np.percentile(FstZW, 2.5), np.percentile(FstZW, 97.5)]
    Tff        =  TZZ
    Tbar_ff    =  np.mean(Tff)
    Tff_ci     =  [np.percentile(Tff, 2.5), np.percentile(Tff, 97.5)]
    Tmm        =  Tt
    Tbar_fm    =  np.mean(Tmm)
    Tmm_ci     =  [np.percentile(Tmm, 2.5), np.percentile(Tmm, 97.5)]
    Tfm        =  (TZZ + TZW)/2
    Tbar_fm    =  np.mean(Tfm)
    Tfm_ci     =  [np.percentile(Tfm, 2.5), np.percentile(Tfm, 97.5)]
    TwFM       =  (Tff + Tmm)/2
    Tbar_wFM   =  np.mean(TwFM)
    TwFM_ci    =  [np.percentile(TwFM, 2.5), np.percentile(TwFM, 97.5)]
    TtFM       =  (Tff + Tmm + 2*Tfm)/4
    Tbar_tFM   =  np.mean(TtFM)
    TtFM_ci    =  [np.percentile(TtFM, 2.5), np.percentile(TtFM, 97.5)]
    FstFM      =  1 - (TwFM/TtFM)
    FstFM      =  FstFM[np.isfinite(FstFM)]
    FstBarFM   =  np.mean(FstFM)
    FstFM_ci   =  [np.percentile(FstFM, 2.5), np.percentile(FstFM, 97.5)]

    # Theoretical predictions
    T_ZZ     =  (9 + (2 * (4 * Ne * r_f))) / (8 + (2 * (4 * Ne * r_f)))
    T_ZW     =  1 + (3 / (2*(4 * Ne * r_f)))
    T_WW     =  (5 + (2 * (4 * Ne * r_f))) / (8 + (2 * (4 * Ne * r_f)))
    return [Tbar_ZZ,  Tbar_ZZ_ci[0], Tbar_ZZ_ci[1],
            Tbar_ZW,  Tbar_ZW_ci[0], Tbar_ZW_ci[1],
            Tbar_WW,  Tbar_WW_ci[0], Tbar_WW_ci[1],
            Tbar_w,   Tw_ci[0],      Tw_ci[1],
            Tbar_t,   Tt_ci[0],      Tt_ci[1],
            FstBarZW, FstZW_ci[0],   FstZW_ci[1],
            Tbar_ff,  Tff_ci[0],     Tff_ci[1],
            Tbar_fm,  Tmm_ci[0],     Tmm_ci[1],
            Tbar_fm,  Tfm_ci[0],     Tfm_ci[1],
            Tbar_wFM, TwFM_ci[0],    TwFM_ci[1],
            Tbar_tFM, TtFM_ci[0],    TtFM_ci[1],
            FstBarFM, FstFM_ci[0],   FstFM_ci[1]]

# Calculate average time to coalescence for Z-W samples (constant Ne)
def TbarWin(chromSamp, tree, chromType):
    # get sample numbers to loop over samples in each population
    if chromType == "Z":
        s1 = list(range(chromSamp[0]))
    if chromType == "W":
        s1 = list(range(chromSamp[5]))
    Tbs = []
    
    # Loop over leaves, find pairwise coalescence times
    for i in s1[:(np.size(s1)-1)]:
        s2 = s1[(i+1):]
        for j in s2:
            Tbs.append(tree.tmrca(u=i, v=j))

    # calculate and return mean
    #TbarW = np.mean(list(set(Tbs)))
    return np.mean(Tbs)

# Calculate average time to coalescence for Z-W samples (constant Ne)
def TbarBtw(chromSamp, tree):
    # get sample numbers to loop over samples in each population
    p1 = list(range(chromSamp[2]))
    p2 = list((chromSamp[2]) + range(chromSamp[3]))
    Tbs = []

    # Loop over samples, find coalescence times for each pair of samples
    # from pop1 and pop2  
    for i in p1:
        for j in p2:
            Tbs.append(tree.tmrca(u=i, v=j))

    # calculate and return mean
#    TbarBtw = np.mean(list(set(Tbs)))
    TbarBtw = np.mean(Tbs)
    return TbarBtw


# Calculate average time to coalescence for Z-W samples (constant Ne)
def TbarChangNe(chromSamp, tree):
    # get sample numbers to loop over samples in each population
    p1 = list(range(chromSamp[2]))
    p2 = list((chromSamp[2] + 1) + range(chromSamp[3]))
    Tbs = []

    # Loop over samples, find coalescence times for each pair of samples
    # from pop1 and pop2  
    for i in p1:
        for j in p2:
            Tbs.append(tree.tmrca(u=i, v=j))

    Tbs = list(Tbs)
#    Tbs = list(set(Tbs))
    Tbar12 = []
    for i in Tbs:
        subdf     =  dfNe[dfNe['beginEpochGen'] < Tbs[i]]
        nEpochs   =  np.shape(subdf)[0]-1
        temp      =  np.zeros(nEpochs)
        for j in range(0, nEpochs):
            if j < nEpochs:
                temp[j] = (subdf['beginEpochGen'][j+1] - subdf['beginEpochGen'][j])  / subdf['uniqueNe'][j]
            elif j == nEpochs:
                temp[j] = (Tbs[i] - subdf['beginEpochGen'][j])  / subdf['uniqueNe'][j]
        Tbar12[i] =  np.sum(temp)       
    return np.mean(Tbar12)

#R style seq() function
def seq(from_ = None, to = None, by = None, length = None, along = None):
    #init
    falling = False

    #length
    if not length is None:
        #TODO: implement when needed
        raise ValueError("`length` is not implemented yet")

    #along object
    if not along is None:
        return list(range(1, len(along) + 1))

    #if all numeric None, set defaults
    if from_ is None and to is None and by is None:
        to = 1
        from_ = 1

    #1 argument
    if not from_ is None and to is None and by is None:
        #this is actually to
        to = from_
        from_ = 1

    #determine by and adjustment
    if from_ > to:
        adjustment = -1

        #set by
        if by is None:
            by = -1
    else:
        adjustment = 1

        #set by
        if by is None:
            by = 1

    #impossible by
    if from_ > to and by > 0:
        raise ValueError("`by` cannot be positive when falling direction")
    if from_ < to and by < 0:
        raise ValueError("`by` cannot be negative when increasing direction")
    if by == 0:
        raise ValueError("`by` cannot be zero")

    #from, to
    if not from_ is None and not to is None:

        return list(range(from_, to + adjustment, by))


    #if you get here, something is wrong!
    raise ValueError("Unknown error. Bad input?")


#convenience wrapper
def seq_along(x):
    return(seq(along = x))


# Variable names from rhoFileRev.txt
# real_rate_map_cM_Mb
# smooth_rate_map_cM_Mb
# rhoVals
# Ne_smooth_r
# Ne_real_r
def sexCoal_2deme_PARbound_makePlotData(chromSamp=BlackSamp, rhoVals=cumRho, NeVals=Ne_real_r, num_replicates=10**4, mu=mu):
    # storage
    rLen = np.size(rhoVals)
    Ts = np.zeros(shape=[rLen,36])
    rhoValsToRun = list(np.array(seq_along(rhoVals))[~np.isnan(NeVals)]-1)
    count = 1
    for i in rhoValsToRun:
        if rhoVals[i] == 0:
            count = count + 1
            continue
        res = sexCoal_2deme_PARbound(num_replicates=num_replicates, rho=rhoVals[i], Ne=NeVals[i])
        Ts[i,0:36] = res[0:36]
        print('Progress:',(count/np.size(rhoValsToRun))*100, '%')
        count = count + 1
    return Ts

