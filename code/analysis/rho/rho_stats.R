#*************************************************************************************************
#* rho estimation
#*************************************************************************************************
r_Ne_table <- read.table("../rho_r_Ne_table.txt", header = T)
head(r_Ne_table)
attach(r_Ne_table)
cor.test(rho, smooth_rate_map_cM_Mb)


plot((r_Ne_table$start+r_Ne_table$end)/(2*mb), Ne_smooth_r, pch = 20)

range(r_Ne_table$Ne_smooth_r)
range(r_Ne_table$Ne_real_r, na.rm = T)