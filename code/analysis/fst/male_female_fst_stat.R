mfFST <- read.table("black_male_female_Z_100Kb.windowed.weir.Z.coord.fst", header=F)
mfFST$V8[mfFST$V8 < 0] <- 0
par_mfFST <- mfFST[(mfFST$V5<52193205),]
nonpar_mfFST <- mfFST[(mfFST$V5>52193205),]
mean(par_mfFST$V8)
mean(nonpar_mfFST$V8)
t.test(par_mfFST$V8, nonpar_mfFST$V8)
