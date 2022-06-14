##################################
# Make data for plots PAR-boundary
exec(open("./sexCoal_PARboundary_plot.py").read())



# Run simulations and save results to .csv
# Results combine original rhoFile_rev.txt data and simulation results

# Simulations using real rate data
realValsData = sexCoal_2deme_PARbound_makePlotData(chromSamp=BlackSamp, rhoVals=cumRho, NeVals=Ne_real_r, num_replicates=10**4, mu=mu)
df_realValsData = pd.DataFrame(np.c_[cumRhoFile[0:145],realValsData])
df_realValsData.reset_index(drop=True, inplace=True)
df_realValsData = df_realValsData.rename(columns={0 : "thin.pos", 1 : "pos.Mb", 2 : "thin.rho", 3 : "thin.cumRho", 4 : "thin.Ne", 5 : "thin.Ne.smooth", 6 : "Tbar_ZZ", 7 : "Tbar_ZZ_ci_lo", 8 : "Tbar_ZZ_ci_hi",  9 : "Tbar_ZW",  10 : "Tbar_ZW_ci_lo",  11 : "Tbar_ZW_ci_hi",  12 : "Tbar_WW",  13 : "Tbar_WW_ci_lo",  14 : "Tbar_WW_ci_hi",  15 : "Tbar_w",  16 : "Tw_ci_lo",  17 : "Tw_ci_hi",  18 : "Tbar_t",  19 : "Tt_ci_lo",  20 : "Tt_ci_hi",  21 : "FstBarXY",  22 : "FstBarXY_ci_lo",  23 : "FstBarXY_ci_hi",  24 : "Tbar_ff",  25 : "Tff_ci_lo",  26 : "Tff_ci_hi",  27 : "Tbar_mm",  28 : "Tmm_ci_lo",  29 : "Tmm_ci_hi",  30 : "Tbar_fm",  31 : "Tfm_ci_lo",  32 : "Tfm_ci_hi",  33 : "Tbar_wFM",  34 : "TwFM_ci_lo",  35 : "TwFM_ci_hi",  36 : "Tbar_tFM",  37 : "TtFM_ci_lo",  38 : "TtFM_ci_hi",  39 : "FstBarFM",  40 : "FstFM_ci_lo",  41 : "FstFM_ci_hi"})
print(df_realValsData)
pd.DataFrame.to_csv(df_realValsData, '../../simulation_output/PARboundary_real_PlotData.csv', index = False) 

	





# Results combine original rhoFile_rev.txt data and simulation results
# Simulations using smoothed data
smoothValsData = sexCoal_2deme_PARbound_makePlotData(chromSamp=BlackSamp, rhoVals=cumRho, NeVals=Ne_smooth_r, num_replicates=10**4, mu=mu)
df_smoothValsData = pd.DataFrame(np.c_[cumRhoFile[0:145],smoothValsData])
df_smoothValsData.reset_index(drop=True, inplace=True)
df_smoothValsData = df_smoothValsData.rename(columns={0 : "thin.pos", 1 : "pos.Mb", 2 : "thin.rho", 3 : "thin.cumRho", 4 : "thin.Ne", 5 : "thin.Ne.smooth", 6 : "Tbar_ZZ", 7 : "Tbar_ZZ_ci_lo", 8 : "Tbar_ZZ_ci_hi",  9 : "Tbar_ZW",  10 : "Tbar_ZW_ci_lo",  11 : "Tbar_ZW_ci_hi",  12 : "Tbar_WW",  13 : "Tbar_WW_ci_lo",  14 : "Tbar_WW_ci_hi",  15 : "Tbar_w",  16 : "Tw_ci_lo",  17 : "Tw_ci_hi",  18 : "Tbar_t",  19 : "Tt_ci_lo",  20 : "Tt_ci_hi",  21 : "FstBarXY",  22 : "FstBarXY_ci_lo",  23 : "FstBarXY_ci_hi",  24 : "Tbar_ff",  25 : "Tff_ci_lo",  26 : "Tff_ci_hi",  27 : "Tbar_mm",  28 : "Tmm_ci_lo",  29 : "Tmm_ci_hi",  30 : "Tbar_fm",  31 : "Tfm_ci_lo",  32 : "Tfm_ci_hi",  33 : "Tbar_wFM",  34 : "TwFM_ci_lo",  35 : "TwFM_ci_hi",  36 : "Tbar_tFM",  37 : "TtFM_ci_lo",  38 : "TtFM_ci_hi",  39 : "FstBarFM",  40 : "FstFM_ci_lo",  41 : "FstFM_ci_hi"})
print(df_smoothValsData)
pd.DataFrame.to_csv(df_smoothValsData, '../../simulation_output/PARboundary_smooth_PlotData.csv', index = False)
