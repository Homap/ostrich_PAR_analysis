# Statistical analyses presented in the manuscript

# Homa Papoli Yazdi

# --------------------------
# H1: LD is lower where rho is higher
# --------------------------
library(nlme)
scaled_ld <- scale()
scaled_rho <- scale()
d <- data.frame(scaled_LD, scaled_rho)

d$obs = as.numeric(rownames(d))
full_PAR <- gls(scaled_ld ~ + scaled_rho, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = d )

# --------------------------
# H1: LD decays faster in the PAR and autosomal regions than in the nonPAR
# --------------------------

pacman::p_load(nlme,ggplot2,jtools,car)

# Autosome
ld_decay_autosome <- read.table("./ld/ld_decay_data/chr4_5.bin")
colnames(ld_decay_autosome) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# PAR
ld_decay_par <- read.table("./ld/ld_decay_data/PAR.bin")
colnames(ld_decay_par) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")
# nonPAR
ld_decay_nonpar <- read.table("./ld/ld_decay_data/nonPAR.bin")
colnames(ld_decay_nonpar) = c("Dist", "Mean_r2", "Mean_Dprime", "Sum_r2", "Sum_Dprime", "NumberPairs")

# Summary statistics
autosome <- unclass(summary(ld_decay_autosome$Mean_r2))
PAR <- unclass(summary(ld_decay_par$Mean_r2))
nonPAR <- unclass(summary(ld_decay_nonpar$Mean_r2))

ld_decay_3_cat <- cbind(data.frame(autosome, check.names = FALSE, stringsAsFactors = FALSE),
  data.frame(PAR, check.names = FALSE, stringsAsFactors = FALSE), 
  data.frame(nonPAR, check.names = FALSE, stringsAsFactors = FALSE))


write.table(ld_decay_3_cat, file=paste("ld_decay.tab"), sep='\t', quote=FALSE)

# Question: Number of pairs within bins: downsample to keep the same across categories or is fine as it is?

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

head(ld_df)
# GLM with binomial 

#Full model
par_autosome <- glm(cbind(r2_events,r2_trials) ~ A_PAR+ log(Dist) + log(Dist):A_PAR, data = ld_df,family=binomial)

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

# Next to do: 
# Make the bins smaller, specially for the first part containing the shortest distance and re-run the model.

# --------------------------
# H1: Genetic diveristy is impacted by recombination and genomic features. 
# Genetic diversity is lower in places of low rho and high gene density. 
# --------------------------

