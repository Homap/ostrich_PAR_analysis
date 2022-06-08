# Dos and HKA
dnds <- read.table("../divergence/ostrich_dnds_pariwise_txt", header = T)
head(dnds)


boxplot(black_A_pnps$ps, black_PAR_pnps$ps, black_nonPAR_pnps$ps, outline = F, ylab = "Ps", names = c("Autosome Ps", "PAR Ps", "nonPAR Ps"))
boxplot(black_A_pnps$pn, black_PAR_pnps$pn, black_nonPAR_pnps$pn, outline = F, ylab = "Pn", names = c("Autosome Pn", "PAR Pn", "nonPAR Pn"))
boxplot(black_A_pnps$pnps, black_PAR_pnps$pnps, black_nonPAR_pnps$pnps, outline = F, ylab = "PnPs", names = c("Autosome PnPs", "PAR PnPs", "nonPAR PnPs"))
median(black_A_pnps$ps)
median(black_PAR_pnps$ps)
median(black_nonPAR_pnps$ps)

median(black_A_pnps$pn)
median(black_PAR_pnps$pn)
median(black_nonPAR_pnps$pn)

median(black_A_pnps$pnps)
median(black_PAR_pnps$pnps)
median(black_nonPAR_pnps$pnps)

# dnds
dnds <- read.table("../divergence/ostrich_dnds_pariwise_txt", header = T)
head(dnds)
PAR_scaffolds = c("superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36")
nonPAR_scaffolds = c("superscaffold69-1", "superscaffold67", "superscaffold62", "superscaffold63", "superscaffold88", "superscaffold93")
dnds_PAR <- dnds[dnds$Scaffold %in% PAR_scaffolds,]
dnds_nonPAR <- dnds[dnds$Scaffold %in% nonPAR_scaffolds,]

boxplot(dnds_PAR$dS, dnds_nonPAR$dS, outline = F, ylab = "ds", names = c("PAR ds", "nonPAR ds"))
boxplot(dnds_PAR$dN, dnds_nonPAR$dN, outline = F, ylab = "dn", names = c("PAR dn", "nonPAR dn"))
boxplot(black_A_pnps$pnps, black_PAR_pnps$pnps, black_nonPAR_pnps$pnps, dnds_PAR$dNdS, dnds_nonPAR$dNdS, outline = F, ylab = "dnds", names = c("PAR dnds", "nonPAR dnds"))

median(dnds_PAR$dS)
median(dnds_nonPAR$dS)
median(dnds_PAR$dN)
median(dnds_nonPAR$dN)
median(dnds_PAR$dNdS)
median(dnds_nonPAR$dNdS)


DoS = Dn/(Dn + Ds) - Pn/(Pn + Ps)

PAR <- read.table("black.PAR.pnps.number_syn_numer_nonsyn.txt", header = T)
nonPAR <- read.table("black.nonPAR.pnps.number_syn_numer_nonsyn.txt", header = T)
PAR[,2] <- sub("_XP_rna", "", PAR[,2])
dnds[,2] <- sub("XM_", "", dnds[,2])
nonPAR[,2] <- sub("_XP_rna", "", nonPAR[,2])

pnps_dnds_list <- dnds[dnds$geneID %in% PAR$gene,]
pnps_pnps_list <- PAR[PAR$gene %in% pnps_dnds_list$geneID,]

nonPAR_pnps_dnds_list <- dnds[dnds$geneID %in% nonPAR$gene,]
nonPAR_pnps_pnps_list <- nonPAR[nonPAR$gene %in% nonPAR_pnps_dnds_list$geneID,]

DoS <- function(Dn, Ds, Pn, Ps){
  DoS_cal = Dn/(Dn + Ds) - Pn/(Pn + Ps)
  return(DoS_cal)
}

alpha <- function(Dn, Ds, Pn, Ps){
  alpha_cal <- 1 - (Ds*Pn)/(Dn*Ps) 
  return(alpha_cal)
}

DoS_PAR <- c()
c = 0
Pn_t = 0
Ps_t = 0
Dn_t = 0
Ds_t = 0
positive_DoS = c()
for(gene in pnps_dnds_list$geneID){
  c = c + 1
  Pn = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NN
  Ps = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NS
  N = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$N
  S = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$S
  dN = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dN
  dS = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dS
  Dn = N * dN
  Ds = S * dS
  
  DoS_PAR[c] = DoS(Dn, Ds, Pn, Ps)
  if(DoS(Dn, Ds, Pn, Ps) > 0){
    print(c(gene, DoS(Dn, Ds, Pn, Ps)))
    positive_DoS[c] = gene
  }
  Pn_t = Pn_t + Pn
  Ps_t = Ps_t + Ps
  Dn_t = Dn_t + Dn
  Ds_t = Ds_t + Ds
  #print(c(Pn_t, Ps_t, Dn_t, Ds_t))
}

c = 0
gene_id = c()
L = c()
Seg = c()
n = c()
D = c()
startingtheta = c()
inheritance = c()
for(gene in pnps_dnds_list$geneID){
  c = c + 1
  gene_id[c] = gene
  L[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$whole_syn
  Seg[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$NS
  n[c] = 20
  S = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$S
  dS = pnps_dnds_list[pnps_dnds_list$geneID==gene,]$dS
  D[c] = S * dS
  startingtheta[c] = pnps_pnps_list[pnps_pnps_list$gene==gene,]$ps
  inheritance[c] = 1
}

HKA_dataframe <- data.frame(gene_id, L, Seg, n, D, startingtheta, inheritance)
write.table(HKA_dataframe, file="HKA_dataframe", sep = "\t", quote=FALSE, row.names = FALSE)


PAR[PAR$gene %in% positive_DoS,]

alpha <- function(Dn, Ds, Pn, Ps){
  alpha_cal <- 1 - (Ds*Pn)/(Dn*Ps) 
  return(alpha_cal)
}

alpha(4508.9460, 7627.0692, 94.5464, 60.6771)


alpha_PAR[c] = alpha(Dn, Ds, Pn, Ps)

hist(DoS_PAR, breaks = 10)

median(DoS_PAR)

length(which(DoS_PAR<0))/length(DoS_PAR)

DoS_nonPAR <- c()
c = 0
for(gene in nonPAR_pnps_dnds_list$geneID){
  c = c + 1
  Pn = nonPAR_pnps_pnps_list[nonPAR_pnps_pnps_list$gene==gene,]$NN
  Ps = nonPAR_pnps_pnps_list[nonPAR_pnps_pnps_list$gene==gene,]$NS
  N = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$N
  S = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$S
  dN = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$dN
  dS = nonPAR_pnps_dnds_list[nonPAR_pnps_dnds_list$geneID==gene,]$dS
  Dn = N * dN
  Ds = S * dS
  print(c(Pn, Ps, Dn, Ds))
  print(DoS(Dn, Ds, Pn, Ps))
  DoS_nonPAR[c] = DoS(Dn, Ds, Pn, Ps)
}

median(DoS_nonPAR)

length(which(DoS_nonPAR<0))/length(DoS_nonPAR)

hist(DoS_nonPAR, breaks = 10)




#********************************************************************************************************************************************
#********************************************************************************************************************************************
#********************************************************************************************************************************************
#********************************************************************************************************************************************
