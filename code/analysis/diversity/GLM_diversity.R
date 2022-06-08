#*************************************************************************************************
#*************************************************************************************************
# GLS
#*************************************************************************************************
#*************************************************************************************************
# PAR section 1
plot((gc_density_500Kb$Window_GC_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
cor.test((gc_density_500Kb$Window_GC_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
plot((gc_density_500Kb$Feat_Base_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]), Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])
cor.test((gc_density_500Kb$Feat_Base_count[1:39]/gc_density_500Kb$Window_Base_count[1:39]),  Z_500Kb$pi[1:39]/Z_500Kb$winsize[1:39])

# PAR section 2
plot((gc_density_500Kb$Window_GC_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
cor.test((gc_density_500Kb$Window_GC_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
plot((gc_density_500Kb$Feat_Base_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])
cor.test((gc_density_500Kb$Feat_Base_count[40:106]/gc_density_500Kb$Window_Base_count[40:106]), Z_500Kb$pi[40:106]/Z_500Kb$winsize[40:106])

# nonPAR
cor.test(gc_density_500Kb$Window_GC_count[106:167]/gc_density_500Kb$Window_Base_count[106:167], Z_500Kb$pi[106:167]/Z_500Kb$winsize[106:167])
cor.test(gc_density_500Kb$Feat_Base_count[106:167]/gc_density_500Kb$Window_Base_count[106:167], Z_500Kb$pi[106:167]/Z_500Kb$winsize[106:167])

# GLS
library(nlme)
scaled_Z <- scale(Z_500Kb$pi/Z_500Kb$winsize)
scaled_gc <- scale(gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count)
scaled_gene <- scale(gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count)
d <- data.frame(scaled_Z, scaled_gc, scaled_gene, scaled_rec)

cor_strut <- pcr_df 
d$obs = as.numeric(rownames(d))
full_PAR <- gls( scaled_Z ~  +  scaled_gc + scaled_gene + scaled_rec, na.action="na.exclude", method = "ML", correlation = corAR1(form =~ obs), data = d )

summary(full_PAR)

scaled_rec <- scale(rec)
diversity <- Z_500Kb$pi/Z_500Kb$winsize
gc <- gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count
gene_dense <- gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count
rec <- test

cor.test(diversity, gc)
cor.test(diversity, gene_dense)
cor.test(diversity, rec)

cor.test(gc, rec)

pcr_df <- data.frame(diversity, gc, gene_dense, rec)

library(pls)
set.seed(1)
model <- pcr(diversity~gc+gene_dense+rec, data=pcr_df, scale=TRUE, validation="CV")

summary(model)

validationplot(model, val.type="R2", cex.axis=0.7)
axis(side = 1, at = c(8), cex.axis=0.7)
abline(v = 8, col = "blue", lty = 3)

validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")

# The model improves by adding three PCAs

train <- pcr_df[1:130, c("diversity", "gc", "gene_dense", "rec")]
y_test <- pcr_df[131:nrow(pcr_df), c("diversity")]
test <- pcr_df[131:nrow(pcr_df), c("gc", "gene_dense", "rec")]
model <- pcr(diversity~gc+gene_dense+rec, data=train, scale=TRUE, validation="CV")
pcr_pred <- predict(model, test, ncomp=3)
sqrt(mean((pcr_pred - y_test)^2))

library(factoextra)
d <- prcomp(~Z_500Kb$pi/Z_500Kb$winsize+gc_density_500Kb$Window_GC_count/gc_density_500Kb$Window_Base_count+
              gc_density_500Kb$Feat_Base_count/gc_density_500Kb$Window_Base_count+test, scale = TRUE)
fviz_eig(d)
fviz_pca_ind(d,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
