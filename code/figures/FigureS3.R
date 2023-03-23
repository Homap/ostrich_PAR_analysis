#*************************************************************************************************************
# Change-point analysis
#*************************************************************************************************************
# Load library for change-point analysis
library(segmented)
#*************************************************************************************************************
mb = 10^6
kb = 10^3

rho_data <- read.table("../../data/rho/ldhat_rho/z/200Kb50Kb_rho_Z.chr.coord.txt", header=F)
colnames(rho_data) <- c("Chr", "Window_start", "Window_end", "scaffold", "start", "end", "rho_per_site", "rho_per_window")
rho_dataset <- data.frame(pos = (rho_data$Window_start+rho_data$Window_end)/(2*mb), rho_per_window=rho_data$rho_per_window/kb)

x = rho_data$Window_start
y = rho_data$rho_per_window/1000

df_Z <- data.frame(x, y)
na.omit(df_Z)

fit_lm = lm(y ~ 1 + x, data = df_Z)  # intercept-only model
fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi = 3)  # Two change points along x

summary(fit_segmented)

pdf("../../figures/FigureS3.pdf", height = 5, width = 8)
plot(fit_segmented, xlab = "Position (Mb)", ylab = "rho/Kb")
points(df_Z, pch = 20, col = "grey")
lines.segmented(fit_segmented)
points.segmented(fit_segmented, col = "red")
dev.off()


