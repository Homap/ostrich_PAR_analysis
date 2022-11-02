# Statistical analyses presented in the manuscript

#-----------
# Clean the workspace
rm(list = ls())
#-----------
# Load libraries
pacman::p_load(nlme,ggplot2,jtools,car)

#----------
mb = 10^6
kb = 10^3
#----------

# --------------------------
# Comparison between male and female recombination rate
# --------------------------

female_map <- read.table("../../data/geneticmap/female.1MB.window.Z.txt", header = T)
male_map <- read.table("../../data/geneticmap/male.1MB.window.Z.txt", header = T)

female_map_par <- female_map[female_map$Window_start<52000000,]
male_map_par <- male_map[male_map$Window_start<52000000,]

mean_female_par <- mean(female_map_par$CM_per_bp*mb, na.rm = T)
mean_male_par <- mean(male_map_par$CM_per_bp*mb, na.rm = T)

mean_female_first_five_par <- mean(female_map_par[female_map_par$Window_start<6000000,]$CM_per_bp*mb, na.rm = T)
mean_female_last_five_par <- mean(female_map_par[female_map_par$Window_start>46000000,]$CM_per_bp*mb, na.rm = T)

mean_male_first_filve_par <- mean(male_map_par[male_map_par$Window_start<6000000,]$CM_per_bp*mb, na.rm = T)
mean_male_last_five_par <- mean(male_map_par[male_map_par$Window_start>46000000,]$CM_per_bp*mb, na.rm = T)

wilcox.test(female_map_par$CM_per_bp, male_map_par$CM_per_bp)

#-----------
# rho vs. cM/Mb
#-----------
rho_mb <- read.table("../../data/rho/ldhat_rho/z/1Mb1Mb.rho.chrom.Z.txt", header=T)
sex_averaged_map <- read.table("../../data/geneticmap/sex_averaged.1MB.window.Z.txt", header=T)

cor.test(rho_mb$rho_per_window[2:79], sex_averaged_map$CM_per_bp[2:79])

rho_par <- rho_mb[rho_mb$Window_start<52000000,]
sex_averaged_map_par <- sex_averaged_map[sex_averaged_map$Window_start<52000000,]

plot(sex_averaged_map_par$CM_per_bp, rho_par$rho_per_window)
plot(female_map_par$CM_per_bp, rho_par$rho_per_window)
plot(male_map_par$CM_per_bp, rho_par$rho_per_window)

cor.test(sex_averaged_map_par$CM_per_bp, rho_par$rho_per_window)
cor.test(female_map_par$CM_per_bp, rho_par$rho_per_window)
cor.test(male_map_par$CM_per_bp, rho_par$rho_per_window)


#smoothed_cm_mb <- loess(sex_averaged_map$CM_per_bp*mb~sex_averaged_map$Window_start, degree = 2, span = 0.15)
#smoothed_cm_mb_predicted <- predict(smoothed_cm_mb)

# --------------------------
# Comparison of rho between the PAR, nonPAR and autosome
# --------------------------
rho_Z <- read.table("../../data/rho/ldhat_rho/z/200Kb200Kb_rho_Z.chr.coord.txt", header=F)
colnames(rho_Z) <- c("Chr", "Start", "End", "scaffold", "scaf_start", "scaf_end", "rho_per_site", "rho_per_window")
rho_chr4 <- read.table("../../data/rho/ldhat_rho/chr4/superscaffold11.200kb200Kb.rho.txt", header = T)
rho_chr5 <- read.table("../../data/rho/ldhat_rho/chr5/superscaffold8.200kb200Kb.rho.txt", header = T)

# Average rho for PAR
rho_par <- rho_Z[rho_Z$Start<52000000,]
rho_nonpar <- rho_Z[rho_Z$Start>52000000,]
rho_a <- rbind(rho_chr4, rho_chr5)
rho_per_window_par <- rho_par$rho_per_window/kb
rho_per_window_nonpar <- rho_nonpar$rho_per_window/kb
rho_per_window_a <- sample(rho_a$rho_per_window/kb, length(rho_per_window_par))

mean(rho_per_window_par)
sd(rho_per_window_par)
mean(rho_per_window_nonpar, na.rm = T)
sd(rho_per_window_nonpar, na.rm = T)
mean(rho_per_window_a)
sd(rho_per_window_a)

# Is rho of the PAR and autosome different from one another?
wilcox.test(rho_per_window_par, rho_per_window_a, alternative = "two.sided")

# --------------------------
# Relationship between rho and LD
# --------------------------
rho_data <- read.table("../../data/rho/ldhat_rho/z/200Kb50Kb_rho_Z.chr.coord.txt", header=F)
colnames(rho_data) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")

ld_data <- read.table("../../data/ld/ld_window/z/Z.LD.05-500.200kbBin.150.tab")
colnames(ld_data) <- c("chr", "Window_start", "Window_end", "Scaffold", "pos1", "pos2", "Dprime", "LOD", "r2", "CIlow", "CIhi", "meanDist", "overlap")

scaled_rho <- scale(rho_data$rho_per_window)
scaled_ld <- scale(ld_data$r2)
rho_ld <- data.frame(scaled_ld, scaled_rho)

rho_ld$obs = as.numeric(rownames(rho_ld))
rho_ld_gls <- gls(scaled_ld ~ + scaled_rho, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = rho_ld)

summary(rho_ld_gls)

# --------------------------
# LD decay
# --------------------------
# Autosome
ld_decay_autosome <- read.table("../../data/ld/ld_decay/autosome/chr4_5.bin")
colnames(ld_decay_autosome) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# PAR
ld_decay_par <- read.table("../../data/ld/ld_decay/z/par/PAR.bin")
colnames(ld_decay_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# nonPAR
ld_decay_nonpar <- read.table("../../data/ld/ld_decay/z/nonpar/nonPAR.bin")
colnames(ld_decay_nonpar) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")

# Summary statistics
autosome <- unclass(summary(ld_decay_autosome$Mean_r2))
PAR <- unclass(summary(ld_decay_par$Mean_r2))
nonPAR <- unclass(summary(ld_decay_nonpar$Mean_r2))

mean(ld_decay_autosome$Mean_r2)
sd(ld_decay_autosome$Mean_r2)
mean(ld_decay_par$Mean_r2)
sd(ld_decay_par$Mean_r2)
mean(ld_decay_nonpar$Mean_r2)
sd(ld_decay_nonpar$Mean_r2)

# LD decay for SNPs 50 Kb apart
ld_decay_autosome[ld_decay_autosome$Dist==50000,]

# Create data frame for GLM
ld_mat = matrix(ncol = 3, nrow = 6020)
for(i in seq(1,3010)){
  ld_mat[i+i-1,1] <- as.numeric(ld_decay_autosome[i,1])
  ld_mat[i+i,1] <- as.numeric(ld_decay_par[i,1])

  ld_mat[i+i-1,2] <- "A"
  ld_mat[i+i,2] <- "PAR"

  ld_mat[i+i-1,3] <- ld_decay_autosome[i,2]
  ld_mat[i+i,3] <- ld_decay_par[i,2]
}
# Data set
ld_df <- data.frame(ld_mat, stringsAsFactors = TRUE)
colnames(ld_df) <- c("Dist", "A_PAR", "r2")
ld_df$Dist<-as.numeric(as.vector(ld_df$Dist))
ld_df$r2<-as.numeric(as.vector(ld_df$r2))
ld_df$r2_events<-round(ld_df$r2*100,0)
ld_df$r2_trials<-100-ld_df$r2_events

# --------------------------
# GLM with binomial
# --------------------------
#Full model
par_autosome <- glm(cbind(r2_events,r2_trials) ~ A_PAR + log(Dist) + log(Dist):A_PAR, data = ld_df,family=binomial)

summ(par_autosome)
Anova(par_autosome, test = "Wald",type=3)

#Model fit to the data
ggplot(ld_df,aes(x=log(Dist), y=r2_events/(r2_events+r2_trials),fill=A_PAR,colour=A_PAR,succ=r2_events, fail=r2_trials))+
  geom_point()+
  geom_smooth(method="glm",method.args=list(family="binomial"),formula = cbind(succ,fail) ~ log(x))+
  ylab("r2")+
  xlab("Distance (log)")+
  labs(fill="Chromosome")+
  guides(colour="none")+
  theme_classic()

#Main effects
par_autosome_1 <- glm(cbind(r2_events,r2_trials) ~ A_PAR+ log(Dist), data = ld_df, family = binomial)

summ(par_autosome_1)
Anova(par_autosome_1, test = "Wald",type=3)

# --------------------------
# Genetic diversity and rho, GC and CDS
# --------------------------
chr4_pi <- read.table("../../data/diversity/genetic_variation/chr4/chr4.200000.sfs.txt", header=T)
mean(chr4_pi$pi/chr4_pi$win_length, na.rm = T)

chr5_pi <- read.table("../../data/diversity/genetic_variation/chr5/chr5.200000.sfs.txt", header=T)
mean(chr5_pi$pi/chr5_pi$win_length, na.rm = T)

a_pi <- rbind(chr4_pi, chr5_pi)
mean(a_pi$pi/a_pi$win_length)
sd(a_pi$pi/a_pi$win_length)

Z_pi <- read.table("../../data/diversity/genetic_variation/z/Z.200000.chr.coord.sfs.txt", header=F)
colnames(x = Z_pi) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_mid",
                       "pi", "theta", "Td", "win_length", "pi_resamp_mean", "pi_resamp_CI_low", "pi_resamp_CI_up")

par_pi <- read.table("../../data/diversity/genetic_variation/par/PAR.200000.sfs.txt", header=T)
nonpar_pi <- read.table("../../data/diversity/genetic_variation/nonpar/nonPAR.1000000.sfs.txt", header=T)

mean(Z_pi$pi/Z_pi$win_length, na.rm = T)
mean(par_pi$pi/par_pi$win_length, na.rm = T)
sd(par_pi$pi/par_pi$win_length, na.rm = T)
mean(nonpar_pi$pi/nonpar_pi$win_length, na.rm = T)
sd(nonpar_pi$pi/nonpar_pi$win_length, na.rm = T)

wilcox.test(par_pi$pi/par_pi$win_length, a_pi$pi/a_pi$win_length)
wilcox.test(par_pi$pi/par_pi$win_length, nonpar_pi$pi/nonpar_pi$win_length)

rho <- read.table("../../data/rho/ldhat_rho/z/200Kb200Kb_rho_Z.chr.coord.txt")
colnames(rho) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")
GC <- read.table("../../data/genomic_features/z_scaf.200000.intergene.overlap.density.sorted.coord.txt", header=F)
colnames(GC) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                 "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
CDS <- read.table("../../data/genomic_features/z_scaf.200000.CDS.overlap.density.sorted.coord.txt", header=F)
colnames(CDS) = c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start", "Window_end", "Window_Base_count",
                  "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

par_rho <- rho[rho$Window_start < 52000000,]
nonpar_rho <- rho[rho$Window_start > 52000000,]

par_pi <- Z_pi[Z_pi$Chr_start < 52000000,]
nonpar_pi <- Z_pi[Z_pi$Chr_start > 52000000,]

par_gc <- GC[GC$Chr_start < 52000000,]
nonpar_gc <- GC[GC$Chr_start > 52000000,]

par_cds <- CDS[CDS$Chr_start < 52000000,]
nonpar_cds <- CDS[CDS$Chr_start > 52000000,]

cds_gc_rho_pi_dataset <- data.frame( cds = c(par_cds$Feat_Base_count/par_cds$Window_Base_count, nonpar_cds$Feat_Base_count/nonpar_cds$Window_Base_count),
                                     gc = c(par_gc$Feat_GC_count/par_gc$Feat_Base_count, nonpar_gc$Feat_GC_count/nonpar_gc$Feat_Base_count),
                                     rho = c(par_rho$rho_per_window, nonpar_rho$rho_per_window),
                                     pi = c(par_pi$pi/par_pi$win_length, nonpar_pi$pi/nonpar_pi$win_length),
                                     chr = c(rep("PAR", 258), rep("nonPAR", 139)))
# ----------------------------------------------------
# GLS model of pi and cds, gc, rho for Z
# ----------------------------------------------------
scaled_pi <- scale(par_pi$pi/par_pi$win_length)
scaled_gc <- scale(par_gc$Feat_GC_count/par_gc$Feat_Base_count)
scaled_cds <- scale(par_cds$Feat_Base_count/par_cds$Window_Base_count)
scaled_rho <- scale(par_rho$rho_per_window)
scaled_dataset <- data.frame(scaled_pi, scaled_gc, scaled_cds, scaled_rho)

scaled_dataset$obs = as.numeric(rownames(scaled_dataset))
full_PAR_gls <- gls( scaled_pi ~  +  scaled_gc + scaled_cds + scaled_rho, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = scaled_dataset )

summary(full_PAR_gls)
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
nonpar_scaled_pi <- scale(nonpar_pi$pi/nonpar_pi$win_length)
nonpar_scaled_gc <- scale(nonpar_gc$Feat_GC_count/nonpar_gc$Feat_Base_count)
nonpar_scaled_cds <- scale(nonpar_cds$Feat_Base_count/nonpar_cds$Window_Base_count)
nonpar_scaled_rho <- scale(nonpar_rho$rho_per_window)
nonpar_scaled_dataset <- data.frame(nonpar_scaled_pi, nonpar_scaled_gc, nonpar_scaled_cds, nonpar_scaled_rho)

nonpar_scaled_dataset$obs = as.numeric(rownames(nonpar_scaled_dataset))
full_nonPAR_gls <- gls( nonpar_scaled_pi ~  +  nonpar_scaled_gc + nonpar_scaled_cds + nonpar_scaled_rho, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = nonpar_scaled_dataset )

summary(full_nonPAR_gls)
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
z_scaled_pi <- scale(Z_pi$pi/Z_pi$win_length)
z_scaled_gc <- scale(GC$Feat_GC_count/GC$Feat_Base_count)
z_scaled_cds <- scale(CDS$Feat_Base_count/CDS$Window_Base_count)
z_scaled_rho <- scale(rho$rho_per_window)
z_scaled_dataset <- data.frame(z_scaled_pi, z_scaled_gc, z_scaled_cds, z_scaled_rho)

z_scaled_dataset$obs = as.numeric(rownames(z_scaled_dataset))
full_Z_gls <- gls( z_scaled_pi ~  +  z_scaled_gc + z_scaled_cds + z_scaled_rho, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = z_scaled_dataset )

summary(full_Z_gls)
# ----------------------------------------------------
# ----------------------------------------------------
# Tajima's D and Fst
# ----------------------------------------------------
chr4_pi <- read.table("../../data/diversity/genetic_variation/chr4/chr4.200000.sfs.txt", header=T)
chr5_pi <- read.table("../../data/diversity/genetic_variation/chr5/chr5.200000.sfs.txt", header=T)

pi_a <- rbind(chr4_pi, chr5_pi)
td_a_mean <- mean(pi_a$Td)
td_par_mean <- mean(par_pi$Td)

wilcox.test(pi_a$Td, par_pi$Td)

#-----------------------
# Fst
#-----------------------
fst <- read.table("../../data/fst/z/black_male_female_Z_200Kb.windowed.weir.Z.coord.fst")
colnames(fst) <- c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start",
                   "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst$WEIGHTED_FST[fst$WEIGHTED_FST<0] <- 0

chr4_fst <- read.table("../../data/fst/autosome/chr4/chr4.200Kb.windowed.weir.fst", header = F)
chr5_fst <- read.table("../../data/fst/autosome/chr5/chr5.200Kb.windowed.weir.fst", header = F)

fst_a <- rbind(chr4_fst, chr5_fst)
colnames(fst_a) <- c("Scaffold", "Window_start", "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst_a$WEIGHTED_FST[fst_a$WEIGHTED_FST<0] <- 0

fst_a_mean <- mean(fst_a$WEIGHTED_FST)

par_fst <- fst[fst$Chr_start < 52000000,]
nonpar_fst <- fst[fst$Chr_start > 52000000,]

mean(par_fst$WEIGHTED_FST)
mean(nonpar_fst$WEIGHTED_FST)

wilcox.test(fst_a$WEIGHTED_FST, par_fst$WEIGHTED_FST)
wilcox.test(par_fst$WEIGHTED_FST, nonpar_fst$WEIGHTED_FST)

# ----------------------------------------------------
# Input for simulation
# ----------------------------------------------------
# Whole PAR
# ----------------------------------------------------
rho_mb_par <- read.table("../../data/rho/ldhat_rho/z/1Mb1Mb.rho.chrom.PAR.txt", header=T)
sex_averaged_map_par <- read.table("../../data/geneticmap/sex_averaged.1MB.window.PAR.txt", header=T)
Ne <- rho_mb_par$rho_per_window/(4*sex_averaged_map_par$kosambi_r_per_bp*mb)

start <- rho_mb_par$Window_start[2:52]
end <- rho_mb_par$Window_end[2:52]
real_rate_map_cM_Mb <- sex_averaged_map_par$CM_per_bp[2:52]*mb

smoothed_kosambi <- loess(sex_averaged_map_par$kosambi_r_per_bp[2:52]*mb~sex_averaged_map_par$Window_start[2:52], span = 0.25)
smoothed_kosambi_predicted <- predict(smoothed_kosambi)

smooth_rate_map_cM_Mb <- predict(loess(sex_averaged_map_par$CM_per_bp[2:52]*mb~sex_averaged_map_par$Window_start[2:52], span = 0.25))

rho <- rho_mb_par$rho_per_window[2:52]
Ne_smooth_r <- rho/(4*smoothed_kosambi_predicted)
sex_averaged_map_par$kosambi_r_per_bp[2:52][which(sex_averaged_map_par$kosambi_r_per_bp[2:52]*mb<0.001)] <- NA
Ne_real_r <- rho/(4*sex_averaged_map_par$kosambi_r_per_bp[2:52]*mb)


# start	end	real_rate_map_cM_Mb	smooth_rate_map_cM_Mb	rho	Ne_smooth_r	Ne_real_r

Ne_dataframe_fullPAR <- data.frame(start, end, real_rate_map_cM_Mb, smooth_rate_map_cM_Mb , rho,
                                   Ne_smooth_r, Ne_real_r)

write.table(Ne_dataframe_fullPAR, "../../data/simulation_input/rho_r_Ne_table.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# ----------------------------------------------------
# PAR boundary
# ----------------------------------------------------
par_boundary_rho <- read.table("../../data/simulation_input/boundary.PAR.txt", skip = 3)
colnames(par_boundary_rho) <- c("Chr", "chr_start", "chr_end", "Scaffold", "scaf_start", "scaf_end", "Mean_rho", "Median", "L95", "U95")

real_rate_map_cM_Mb = c()
smooth_rate_map_cM_Mb = c()
Ne_smooth_r = c()
Ne_real_r = c()
k = 0

for(i in par_boundary_rho$chr_start) {
  if(i >= 5*10^7 && i <= 5.1*10^7){
    k = k + 1
    real_rate_map_cM_Mb[k] = Ne_dataframe_fullPAR[50,3]
    smooth_rate_map_cM_Mb[k] = Ne_dataframe_fullPAR[50,4]
    Ne_smooth_r[k] = Ne_dataframe_fullPAR[50,6]
    Ne_real_r[k] = Ne_dataframe_fullPAR[50,7]
  }
  else if (i >= 5.1*10^7 && i <= 5.3*10^7){
    k = k + 1
    real_rate_map_cM_Mb[k] = Ne_dataframe_fullPAR[51,3]
    smooth_rate_map_cM_Mb[k] = Ne_dataframe_fullPAR[51,4]
    Ne_smooth_r[k] = Ne_dataframe_fullPAR[51,6]
    Ne_real_r[k] = Ne_dataframe_fullPAR[51,7]
  }
}

start_sc36 <- rev(par_boundary_rho$scaf_start)
end_sc36 <- rev(par_boundary_rho$scaf_end)
chr_start <- rev(par_boundary_rho$chr_start)
chr_end <- rev(par_boundary_rho$chr_end)
rho_sc36 <- rev((par_boundary_rho$chr_end - par_boundary_rho$chr_start)*(par_boundary_rho$Mean_rho/1000)) # Mean_rho is per kb.

rho_scaffold_36 <- data.frame(start_sc36, end_sc36, chr_start, chr_end, rho_sc36, real_rate_map_cM_Mb = rev(real_rate_map_cM_Mb), smooth_rate_map_cM_Mb = rev(smooth_rate_map_cM_Mb), Ne_smooth_r = rev(Ne_smooth_r), Ne_real_r = rev(Ne_real_r))

write.table(rho_scaffold_36, "../../data/simulation_input/rho_scaffold_36.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

