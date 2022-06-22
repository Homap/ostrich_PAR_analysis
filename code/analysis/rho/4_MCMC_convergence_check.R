#!/usr/bin/R

library(coda)

mcmc_summary_table <- function(mcmc_input) {
  mcmc <- read.table(mcmc_input, header = T)
  mcmc_obj <- mcmc(mcmc[-seq(1, 200),c(2,3,4)])
  summary_mcmc <- summary(mcmc_obj)
  ess <- effectiveSize(mcmc_obj)
  likelihood <- as.vector(c(summary_mcmc$statistics[1,1:2], summary_mcmc$quantiles[1,], ess[1]))
  blocks <- as.vector(c(summary_mcmc$statistics[2,1:2], summary_mcmc$quantiles[2,], ess[2]))
  map_length <- as.vector(c(summary_mcmc$statistics[3,1:2], summary_mcmc$quantiles[3,], ess[3]))
  return(list(mcmc_obj, likelihood, blocks, map_length))
}

mcmc_summary_gelman <- function(chrom, scaf_vector, path_to_folder, outputdir) {
  
  output_list <- list(c(), c(), c(), c(), c(), c(), c())
  
  for(scaffold in scaf_vector){
    print(scaffold)
    for(window in seq(1, length(list.files(path_to_folder, pattern = paste("^", scaffold, sep="")))/3)){
      print(window)
      filename1 <- paste(scaffold,window,"1.chainreport.txt", sep=".")
      filename2 <- paste(scaffold,window,"2.chainreport.txt", sep=".")
      filename3 <- paste(scaffold,window,"3.chainreport.txt", sep=".")
      
      mcmc1 <- mcmc_summary_table(paste(path_to_folder, filename1, sep = "/"))
      mcmc2 <- mcmc_summary_table(paste(path_to_folder, filename2, sep = "/"))
      mcmc3 <- mcmc_summary_table(paste(path_to_folder, filename3, sep = "/"))
      
      rowname1 = paste(scaffold,window,"1", sep=".") 
      
      likelihood <- rbind(c(rowname1, mcmc1[[2]]))
      block <- rbind(c(rowname1, mcmc2[[3]]))
      map_length <- rbind(c(rowname1, mcmc1[[4]]))
      
      output_list[[1]] <- rbind(likelihood, output_list[[1]])
      output_list[[2]] <- rbind(block, output_list[[2]])
      output_list[[3]] <- rbind(map_length, output_list[[3]])
      
      mcmc_list <- mcmc.list(mcmc1[[1]], mcmc2[[1]], mcmc3[[1]])
      gelman_diag <- gelman.diag(mcmc_list)
      
      likelihood_diag <- c(paste(scaffold,window, sep="."), as.vector(gelman_diag$psrf[1,]))
      blocks_diag <- c(paste(scaffold,window, sep="."), as.vector(gelman_diag$psrf[2,]))
      map_length_diag <- c(paste(scaffold,window, sep="."), as.vector(gelman_diag$psrf[3,]))
      
      output_list[[4]] <- rbind(likelihood_diag, output_list[[4]])
      output_list[[5]] <- rbind(blocks_diag, output_list[[5]])
      output_list[[6]] <- rbind(map_length_diag, output_list[[6]])
    }
  }
  
  likelihood_summary <- data.frame(output_list[[1]])
  colnames(likelihood_summary) <- c("scaffold", "mean", "SD", "2.5", "25", "50", "75", "97.5", "ESS")
  
  block_summary <- data.frame(output_list[[2]])
  colnames(block_summary) <- c("scaffold", "mean", "SD", "2.5", "25", "50", "75", "97.5", "ESS")
  
  map_summary <- data.frame(output_list[[3]])
  colnames(map_summary) <- c("scaffold", "mean", "SD", "2.5", "25", "50", "75", "97.5", "ESS")
  
  likelihood_gelman <- data.frame(output_list[[4]])
  colnames(likelihood_gelman) <- c("scaffold", "point_estimate", "UPPER_CI")
  
  block_gelman <- data.frame(output_list[[5]])
  colnames(block_gelman) <- c("scaffold", "point_estimate", "UPPER_CI")
  
  map_gelman <- data.frame(output_list[[6]])
  colnames(map_gelman) <- c("scaffold", "point_estimate", "UPPER_CI")

  
  write.table(likelihood_summary, paste(outputdir, paste(chrom,"mcmc.likelihood.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(block_summary, paste(outputdir, paste(chrom,"mcmc.block.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(map_summary, paste(outputdir, paste(chrom,"mcmc.map.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(likelihood_gelman, paste(outputdir, paste(chrom,"gelman.likelihood.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(block_gelman, paste(outputdir, paste(chrom,"gelman.block.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(map_gelman, paste(outputdir, paste(chrom,"gelman.map.txt", sep ="."), sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(NULL)
}


plot_gelman_point_estimate <- function(gelman_likelihood, gelman_block, gelman_map, plot_out) {
    gelman_likelihood_table <- read.table(gelman_likelihood, header = TRUE)
    gelman_block_table <- read.table(gelman_block, header = TRUE)
    gelman_map_table <- read.table(gelman_map, header = TRUE)

    pdf(plot_out)
    boxplot(gelman_likelihood_table$point_estimate,
            gelman_block_table$point_estimate,
            gelman_map_table$point_estimate, names = c("Likelihood", "Block_number", "Map_length"), ylab = "Gelman diag point estimate")

    dev.off()
    return(NULL)
}


print("Chr4")
chr4_scaf_list = c("superscaffold11")
mcmc_summary_gelman(chrom = "chr4", scaf_vector = chr4_scaf_list, path_to_folder = "../../../data/rho/ldhat_mcmc/chr4", outputdir = "../../../data/rho/ldhat_mcmc/")

plot_gelman_point_estimate("../../../data/rho/ldhat_mcmc/chr4.gelman.likelihood.txt",
                           "../../../data/rho/ldhat_mcmc/chr4.gelman.block.txt",
                           "../../../data/rho/ldhat_mcmc/chr4.gelman.map.txt", "../../../data/rho/ldhat_mcmc/chr4_gelman.pdf")

print("Chr5")
chr5_scaf_list = c("superscaffold8")
mcmc_summary_gelman(chrom = "chr5", scaf_vector = chr5_scaf_list, path_to_folder = "../../../data/rho/ldhat_mcmc/chr5", outputdir = "../../../data/rho/ldhat_mcmc/")

plot_gelman_point_estimate("../../../data/rho/ldhat_mcmc/chr5.gelman.likelihood.txt",
                           "../../../data/rho/ldhat_mcmc/chr5.gelman.block.txt",
                           "../../../data/rho/ldhat_mcmc/chr5.gelman.map.txt", "../../../data/rho/ldhat_mcmc/chr5_gelman.pdf")

print("PAR")
par_scaf_list = c("superscaffold26", "superscaffold54", "superscaffold35", "superscaffold35")
mcmc_summary_gelman(chrom = "par", scaf_vector = par_scaf_list, path_to_folder = "../../../data/rho/ldhat_mcmc/z/par", outputdir = "../../../data/rho/ldhat_mcmc/")

plot_gelman_point_estimate("../../../data/rho/ldhat_mcmc/par.gelman.likelihood.txt",
                           "../../../data/rho/ldhat_mcmc/par.gelman.block.txt",
                           "../../../data/rho/ldhat_mcmc/par.gelman.map.txt", "../../../data/rho/ldhat_mcmc/par_gelman.pdf")

print("nonPAR")
nonpar_scaf_list = c("superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92")
mcmc_summary_gelman(chrom = "nonpar", scaf_vector = nonpar_scaf_list, path_to_folder = "../../../data/rho/ldhat_mcmc/z/nonpar", outputdir = "../../../data/rho/ldhat_mcmc/")

plot_gelman_point_estimate("../../../data/rho/ldhat_mcmc/nonpar.gelman.likelihood.txt",
                           "../../../data/rho/ldhat_mcmc/nonpar.gelman.block.txt",
                           "../../../data/rho/ldhat_mcmc/nonpar.gelman.map.txt", "../../../data/rho/ldhat_mcmc/nonpar_gelman.pdf")