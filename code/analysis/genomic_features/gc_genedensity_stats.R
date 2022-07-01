gc_density_100Kb <- read.table("../genomic_features/black.Z.genedensity.500Kb.txt")
colnames(x = gc_density_100Kb) = c("Chr", "Scaffold", "Window_start", "Window_end", "Chr_start", "Chr_end", "Window_Base_count",
                                   "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
# Genomic features
gene_GC_0.1Mb_PAR <- read.table("../genomic_features/black.par_scaf.100Kb.CDS.overlap.density.sorted.txt")
gene_GC_0.5Mb_PAR <- read.table("../genomic_features/black.par_scaf.500Kb.CDS.overlap.density.sorted.txt")
gene_GC_1Mb_PAR <- read.table("../genomic_features/black.par_scaf.1000Kb.CDS.overlap.density.sorted.txt")
colnames(x = gene_GC_0.1Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                    "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_0.5Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                    "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_1Mb_PAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                  "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
gene_GC_0.1Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.100Kb.CDS.overlap.density.sorted.txt")
gene_GC_0.5Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.500Kb.CDS.overlap.density.sorted.txt")
gene_GC_1Mb_nonPAR <- read.table("../genomic_features/black.nonpar_scaf.1000Kb.CDS.overlap.density.sorted.txt")
colnames(x = gene_GC_0.1Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                       "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_0.5Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                       "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")
colnames(x = gene_GC_1Mb_nonPAR) = c("Scaffold", "Window_start", "Window_end", "Window_Base_count", 
                                     "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")

GC_PAR <- read.table("../genomic_features/black.Z.density.txt", header=T)
gc_density_100Kb <- read.table("../genomic_features/black.Z.genedensity.100Kb.txt")

gc_density_1000Kb <- read.table("../genomic_features/black.Z.genedensity.1000Kb.txt")

gc_density_500Kb <- read.table("../genomic_features/black.Z.genedensity.500Kb.txt")
colnames(x = gc_density_500Kb) = c("Chr", "Scaffold", "Window_start", "Window_end", "cstart", "cend", "Window_Base_count", 
                                   "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count")