#!/usr/bin/R

source("../../software_module/R/vcftools_figure_summary.R")
#--------------------------
# Variant based statistics
#--------------------------
print("Variant statistics")
# Autosome
A_var_stat <- variant_statistics("../../data/vcf/a_vcf/stats/A.lqual", "../../data/vcf/a_vcf/stats/A.ldepth.mean", 
                                 "../../data/vcf/a_vcf/stats/A.lmiss", "../../data/vcf/a_vcf/stats/A.frq")
png("../../data/vcf/a_vcf/stats/A.var_qual.png"); print(A_var_stat[[1]]); dev.off()
png("../../data/vcf/a_vcf/stats/A.ldept.png"); print(A_var_stat[[2]]); dev.off()
png("../../data/vcf/a_vcf/stats/A.lmiss.png"); print(A_var_stat[[3]]); dev.off()
png("../../data/vcf/a_vcf/stats/A.frq.png"); print(A_var_stat[[4]]); dev.off()
sink("../../data/vcf/a_vcf/stats/A_summary_site_depth")
print(A_var_stat[[5]])
sink()
# PAR
par_var_stat <- variant_statistics("../../data/vcf/par_vcf/stats/PAR.lqual", "../../data/vcf/par_vcf/stats/PAR.ldepth.mean", 
                                 "../../data/vcf/par_vcf/stats/PAR.lmiss", "../../data/vcf/par_vcf/stats/PAR.frq")
png("../../data/vcf/par_vcf/stats/PAR.var_qual.png"); print(par_var_stat[[1]]); dev.off()
png("../../data/vcf/par_vcf/stats/PAR.ldept.png"); print(par_var_stat[[2]]); dev.off()
png("../../data/vcf/par_vcf/stats/PAR.lmiss.png"); print(par_var_stat[[3]]); dev.off()
png("../../data/vcf/par_vcf/stats/PAR.frq.png"); print(par_var_stat[[4]]); dev.off()
sink("../../data/vcf/par_vcf/stats/PAR_summary_site_depth")
print(par_var_stat[[5]])
sink()
# nonPAR
nonpar_var_stat <- variant_statistics("../../data/vcf/nonpar_vcf/stats/nonPAR.lqual", "../../data/vcf/nonpar_vcf/stats/nonPAR.ldepth.mean", 
                                 "../../data/vcf/nonpar_vcf/stats/nonPAR.lmiss", "../../data/vcf/nonpar_vcf/stats/nonPAR.frq")
png("../../data/vcf/nonpar_vcf/stats/nonPAR.var_qual.png"); print(nonpar_var_stat[[1]]); dev.off()
png("../../data/vcf/nonpar_vcf/stats/nonPAR.ldept.png"); print(nonpar_var_stat[[2]]); dev.off()
png("../../data/vcf/nonpar_vcf/stats/nonPAR.lmiss.png"); print(nonpar_var_stat[[3]]); dev.off()
png("../../data/vcf/nonpar_vcf/stats/nonPAR.frq.png"); print(nonpar_var_stat[[4]]); dev.off()
sink("../../data/vcf/nonpar_vcf/stats/nonPAR_summary_site_depth")
print(nonpar_var_stat[[5]])
sink()
# #-----------------------------
# # Individual based statistics
# #-----------------------------
print("Individual statistics")
# Autosome
A_ind_stat <- individual_statistics("../../data/vcf/a_vcf/stats/A.idepth", "../../data/vcf/a_vcf/stats/A.imiss", "../../data/vcf/a_vcf/stats/A.het")
png("../../data/vcf/a_vcf/stats/A.idepth.png"); print(A_ind_stat[[1]]); dev.off()
png("../../data/vcf/a_vcf/stats/A.imiss.png"); print(A_ind_stat[[2]]); dev.off()
png("../../data/vcf/a_vcf/stats/A.het.png"); print(A_ind_stat[[3]]); dev.off()
# PAR
par_ind_stat <- individual_statistics("../../data/vcf/par_vcf/stats/PAR.idepth", "../../data/vcf/par_vcf/stats/PAR.imiss", "../../data/vcf/par_vcf/stats/PAR.het")
png("../../data/vcf/par_vcf/stats/PAR.idepth.png"); print(par_ind_stat[[1]]); dev.off()
png("../../data/vcf/par_vcf/stats/PAR.imiss.png"); print(par_ind_stat[[2]]); dev.off()
png("../../data/vcf/par_vcf/stats/PAR.het.png"); print(par_ind_stat[[3]]); dev.off()
# nonPAR
nonpar_ind_stat <- individual_statistics("../../data/vcf/nonpar_vcf/stats/nonPAR.idepth", "../../data/vcf/nonpar_vcf/stats/nonPAR.imiss", "../../data/vcf/nonpar_vcf/stats/nonPAR.het")
png("../../data/vcf/nonpar_vcf/stats/nonPAR.idepth.png"); print(nonpar_ind_stat[[1]]); dev.off()
png("../../data/vcf/nonpar_vcf/stats/nonPAR.imiss.png"); print(nonpar_ind_stat[[2]]); dev.off()
png("../../data/vcf/nonpar_vcf/stats/nonPAR.het.png"); print(nonpar_ind_stat[[3]]); dev.off()

