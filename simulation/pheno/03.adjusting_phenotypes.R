#!/usr/bin/env Rscript

library(tidyverse)

#Define helper functions

inverse_normal_transformation <- function(r){
  qnorm((rank(r,na.last="keep")-0.5)/sum(!is.na(r)))
}

pheno_transform <- function(phenotype, .data, .control_vars, .log_list = list()){

  log <- phenotype %in% .log_list
  if(log){
    str <- str_c("log10(", phenotype, " + 1) ~ " ,.control_vars)
  } else {
    str <- paste(phenotype, .control_vars, sep = "~")
  }
  f <- formula(str)
  r <- inverse_normal_transformation(resid(lm(f, data = .data, na.action = "na.exclude")))
  return(r)
}


#Read in data
adjust_pheno_by_id <- function(sim_name, id_type, rep){
  ids <- read_tsv(sprintf("../../data/simulated_phenotypes/%s/COSI_phenos_nloci_1/COSI_sims_nloci_1_%s_ids", 
                          sim_name, id_type), 
                  col_names = c("FID", "IID"))
  
  pheno <- read_tsv(sprintf("../../data/simulated_phenotypes/%s/COSI_phenos_nloci_1/COSI_phenos_nloci_1_rep_%s.pheno", 
                            sim_name, rep)) %>%
    inner_join(ids, by = c("FID", "IID"))
  
  pcs <- read_delim(sprintf("../../data/simulated_genotypes/causal_genotypes/%s/COSI_causal_genotypes_nloci_1_%s_LD_pruned.eigenvec", 
                            sim_name, id_type), 
                    col_names = c("FID", "IID", paste0("genotype_pc", 1:10)), 
                    delim = " ")
  
  adjust_data <- left_join(pheno, pcs, 
                           by = c("FID", "IID"))
  #Create formula string for adjustment
  control_vars <- paste(paste0("genotype_pc", 1:10), collapse = "+")
  
  # run analysis
  
  pheno_list <- colnames(pheno)[str_detect(colnames(pheno), "pheno")]
  names(pheno_list) <- pheno_list
  
  adjusted_phenotypes <- bind_cols(adjust_data %>% select(-any_of(pheno_list)), 
                                       map_dfc(pheno_list, ~pheno_transform(.x, adjust_data))) %>% 
    select(FID, IID, all_of(pheno_list))
  
  print("writing phenotype files")
  write_tsv(adjusted_phenotypes, sprintf("../../data/simulated_phenotypes/%s/COSI_phenos_nloci_1/COSI_phenos_nloci_1_rep_%s_%s_adjusted.pheno", 
                                             sim_name, rep, id_type))
  
}

adjust_pheno <- function(sim_name, rep){
  
  pheno <- read_tsv(sprintf("../../data/simulated_phenotypes/%s/COSI_phenos_nloci_1/COSI_phenos_nloci_1_rep_%s.pheno", 
                            sim_name, rep)) 
  
  pcs <- read_delim(sprintf("../../data/simulated_genotypes/causal_genotypes/%s/COSI_causal_genotypes_nloci_1_LD_pruned.eigenvec", 
                            sim_name), 
                    col_names = c("FID", "IID", paste0("genotype_pc", 1:10)), 
                    delim = " ")
  
  adjust_data <- left_join(pheno, pcs, 
                           by = c("FID", "IID"))
  #Create formula string for adjustment
  control_vars <- paste(paste0("genotype_pc", 1:10), collapse = "+")

  # run analysis
  
  pheno_list <- colnames(pheno)[str_detect(colnames(pheno), "pheno")]
  names(pheno_list) <- pheno_list
  
  adjusted_phenotypes <- bind_cols(adjust_data %>% select(-any_of(pheno_list)), 
                                   map_dfc(pheno_list, ~pheno_transform(.x, adjust_data, control_vars))) %>% 
    select(FID, IID, all_of(pheno_list))
  
  print("writing phenotype files")
  write_tsv(adjusted_phenotypes, sprintf("../../data/simulated_phenotypes/%s/COSI_phenos_nloci_1/COSI_phenos_nloci_1_rep_%s_adjusted.pheno", 
                                         sim_name, rep))
  
}

sim_name <- "ideal_FL"
set.seed(13)
map(1:10, ~adjust_pheno(sim_name = sim_name, rep = .x))

