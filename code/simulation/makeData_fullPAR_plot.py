##################################
# Make data for plots Full-PAR
exec(open("./sexCoal_fullPAR_plot.py").read())

# Run simulations and save results to .csv
# Results combine original rhoFile.txt data and simulation results  ... (Use the original coordinates & remove '' from filenames)

# Simulations using real rate data
realValsData = sexCoal_2deme_fullPAR_makePlotData(chromSamp=BlackSamp, rVals=real_rate_map_cM_Mb, rhoVals=rhoVals, NeVals=Ne_real_r, num_replicates=10**4, mu=mu)
df_realValsData = pd.DataFrame(np.c_[rhoFile,realValsData])
df_realValsData.reset_index(drop=True, inplace=True)
df_realValsData = df_realValsData.rename(columns={0 : "start", 1 : "end", 2 : "real_rate_map_cM_Mb", 3 : "smooth_rate_map_cM_Mb", 4 : "rho", 5 : "Ne_smooth_r", 6 : "Ne_real_r", 7 : "Tbar_ZZ", 8 : "Tbar_ZZ_ci_lo", 9 : "Tbar_ZZ_ci_hi", 10 : "Tbar_ZW", 11 : "Tbar_ZW_ci_lo", 12 : "Tbar_ZW_ci_hi", 13 : "Tbar_WW", 14 : "Tbar_WW_ci_lo", 15 : "Tbar_WW_ci_hi", 16 : "Tbar_w", 17 : "Tw_ci_lo", 18 : "Tw_ci_hi", 19 : "Tbar_t", 20 : "Tt_ci_lo", 21 : "Tt_ci_hi", 22 : "FstBarXY", 23 : "FstBarXY_ci_lo", 24 : "FstBarXY_ci_hi", 25 : "Tbar_ff", 26 : "Tff_ci_lo", 27 : "Tff_ci_hi", 28 : "Tbar_mm", 29 : "Tmm_ci_lo", 30 : "Tmm_ci_hi", 31 : "Tbar_fm", 32 : "Tfm_ci_lo", 33 : "Tfm_ci_hi", 34 : "Tbar_wFM", 35 : "TwFM_ci_lo", 36 : "TwFM_ci_hi", 37 : "Tbar_tFM", 38 : "TtFM_ci_lo", 39 : "TtFM_ci_hi", 40 : "FstBarFM", 41 : "FstFM_ci_lo", 42 : "FstFM_ci_hi"})
print(df_realValsData)
pd.DataFrame.to_csv(df_realValsData, '../../simulation_output/fullPAR_Black_real_PlotData.csv', index = False)
	

# Results combine original rhoFile_.txt data and simulation results
# Simulations using smoothed data
smoothValsData = sexCoal_2deme_fullPAR_makePlotData(chromSamp=BlackSamp, rVals=smooth_rate_map_cM_Mb, rhoVals=rhoVals, NeVals=Ne_smooth_r, num_replicates=10**4, mu=mu)
df_smoothValsData = pd.DataFrame(np.c_[rhoFile,smoothValsData])
df_smoothValsData.reset_index(drop=True, inplace=True)
df_smoothValsData = df_smoothValsData.rename(columns={0 : "start_", 1 : "end_", 2 : "real_rate_map_cM_Mb_", 3 : "smooth_rate_map_cM_Mb_", 4 : "rho_", 5 : "Ne_smooth_r_", 6 : "Ne_real_r", 7 : "Tbar_ZZ", 8 : "Tbar_ZZ_ci_lo", 9 : "Tbar_ZZ_ci_hi", 10 : "Tbar_ZW", 11 : "Tbar_ZW_ci_lo", 12 : "Tbar_ZW_ci_hi", 13 : "Tbar_WW", 14 : "Tbar_WW_ci_lo", 15 : "Tbar_WW_ci_hi", 16 : "Tbar_w", 17 : "Tw_ci_lo", 18 : "Tw_ci_hi", 19 : "Tbar_t", 20 : "Tt_ci_lo", 21 : "Tt_ci_hi", 22 : "FstBarXY", 23 : "FstBarXY_ci_lo", 24 : "FstBarXY_ci_hi", 25 : "Tbar_ff", 26 : "Tff_ci_lo", 27 : "Tff_ci_hi", 28 : "Tbar_mm", 29 : "Tmm_ci_lo", 30 : "Tmm_ci_hi", 31 : "Tbar_fm", 32 : "Tfm_ci_lo", 33 : "Tfm_ci_hi", 34 : "Tbar_wFM", 35 : "TwFM_ci_lo", 36 : "TwFM_ci_hi", 37 : "Tbar_tFM", 38 : "TtFM_ci_lo", 39 : "TtFM_ci_hi", 40 : "FstBarFM", 41 : "FstFM_ci_lo", 42 : "FstFM_ci_hi"})
print(df_smoothValsData)
pd.DataFrame.to_csv(df_smoothValsData, '../../simulation_output/fullPAR_Black_smooth_PlotData.csv', index = False)
