##################################
# Make data for plots Full-PAR
exec(open("./sexCoal_fullPAR_plot.py").read())



# Run simulations and save results to .csv
# Results combine original rhoFile_rev.txt data and simulation results

# Simulations using real rate data
realValsData = sexCoal_2deme_fullPAR_makePlotData(chromSamp=BlackSamp, rVals=real_rate_map_cM_Mb, rhoVals=rhoVals, NeVals=Ne_real_r, num_replicates=10**4, mu=mu)
df_realValsData = pd.DataFrame(np.c_[rhoFileRev,realValsData])
df_realValsData.reset_index(drop=True, inplace=True)
df_realValsData = df_realValsData.rename(columns={0 : "start_rev", 1 : "end_rev", 2 : "real_rate_map_cM_Mb_rev", 3 : "smooth_rate_map_cM_Mb_rev", 4 : "rho_rev", 5 : "Ne_smooth_r_rev", 6 : "Ne_real_r_rev", 7 : "Tbar_ZZ", 8 : "Tbar_ZZ_ci_lo", 9 : "Tbar_ZZ_ci_hi", 10 : "Tbar_ZW", 11 : "Tbar_ZW_ci_lo", 12 : "Tbar_ZW_ci_hi", 13 : "Tbar_WW", 14 : "Tbar_WW_ci_lo", 15 : "Tbar_WW_ci_hi", 16 : "Tbar_w", 17 : "Tw_ci_lo", 18 : "Tw_ci_hi", 19 : "Tbar_t", 20 : "Tt_ci_lo", 21 : "Tt_ci_hi", 22 : "FstBarXY", 23 : "FstBarXY_ci_lo", 24 : "FstBarXY_ci_hi", 25 : "Tbar_ff", 26 : "Tff_ci_lo", 27 : "Tff_ci_hi", 28 : "Tbar_mm", 29 : "Tmm_ci_lo", 30 : "Tmm_ci_hi", 31 : "Tbar_fm", 32 : "Tfm_ci_lo", 33 : "Tfm_ci_hi", 34 : "Tbar_wFM", 35 : "TwFM_ci_lo", 36 : "TwFM_ci_hi", 37 : "Tbar_tFM", 38 : "TtFM_ci_lo", 39 : "TtFM_ci_hi", 40 : "FstBarFM", 41 : "FstFM_ci_lo", 42 : "FstFM_ci_hi"})
print(df_realValsData)
pd.DataFrame.to_csv(df_realValsData, './fullPAR_Black_real_PlotData.csv', index = False)
	

# Results combine original rhoFile_rev.txt data and simulation results
# Simulations using smoothed data
smoothValsData = sexCoal_2deme_fullPAR_makePlotData(chromSamp=BlackSamp, rVals=smooth_rate_map_cM_Mb, rhoVals=rhoVals, NeVals=Ne_smooth_r, num_replicates=10**4, mu=mu)
df_smoothValsData = pd.DataFrame(np.c_[rhoFileRev,smoothValsData])
df_smoothValsData.reset_index(drop=True, inplace=True)
df_smoothValsData = df_smoothValsData.rename(columns={0 : "start_rev", 1 : "end_rev", 2 : "real_rate_map_cM_Mb_rev", 3 : "smooth_rate_map_cM_Mb_rev", 4 : "rho_rev", 5 : "Ne_smooth_r_rev", 6 : "Ne_real_r_rev", 7 : "Tbar_ZZ", 8 : "Tbar_ZZ_ci_lo", 9 : "Tbar_ZZ_ci_hi", 10 : "Tbar_ZW", 11 : "Tbar_ZW_ci_lo", 12 : "Tbar_ZW_ci_hi", 13 : "Tbar_WW", 14 : "Tbar_WW_ci_lo", 15 : "Tbar_WW_ci_hi", 16 : "Tbar_w", 17 : "Tw_ci_lo", 18 : "Tw_ci_hi", 19 : "Tbar_t", 20 : "Tt_ci_lo", 21 : "Tt_ci_hi", 22 : "FstBarXY", 23 : "FstBarXY_ci_lo", 24 : "FstBarXY_ci_hi", 25 : "Tbar_ff", 26 : "Tff_ci_lo", 27 : "Tff_ci_hi", 28 : "Tbar_mm", 29 : "Tmm_ci_lo", 30 : "Tmm_ci_hi", 31 : "Tbar_fm", 32 : "Tfm_ci_lo", 33 : "Tfm_ci_hi", 34 : "Tbar_wFM", 35 : "TwFM_ci_lo", 36 : "TwFM_ci_hi", 37 : "Tbar_tFM", 38 : "TtFM_ci_lo", 39 : "TtFM_ci_hi", 40 : "FstBarFM", 41 : "FstFM_ci_lo", 42 : "FstFM_ci_hi"})
print(df_smoothValsData)
pd.DataFrame.to_csv(df_smoothValsData, './fullPAR_Black_smooth_PlotData.csv', index = False)
