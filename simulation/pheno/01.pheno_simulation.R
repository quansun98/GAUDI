#This script is for the simulated phenotypes using the model:
# y_i = x_ij * b_j + e_i 
#
#
# Heritability: proportion of phenotypic variance explained by SNPs.
# Proportion of Heritability Explained by Ancestry Shared Effects (self explanatory). 
#     As this variable decreases, the variance of the ancestry specific betas increases, 
#     which makes the two ancestral populations less similar. This is why ancestry specific effects
#     can explain more of the phenotypic variance. 
# Proportion of SNPs which are causal
# Total number of SNPs (controlled by # of COSI loci)

library(tidyverse)
library(Hmisc)
library(Matrix)
library(splitTools)
library(matrixStats)

sim_prefix <- "ideal_FL_v2"
r2 <- "0.5"
sim_name <- sprintf("%s_r2%s", sim_prefix, r2)

#Defining helper functions. 
#Scale a matrix such that a column of constats returns 0's, not NaN.
custom_scale_matrix <- function(.matrix){
  sds <- colSds(.matrix)
  sds[which(sds == 0)] <- 1
  scale(.matrix, center = T, scale = sds) %>% 
    return()
}

assign_loci <- function(.sim_setup){
  
  distinct_loci <- .sim_setup %>% distinct(n_loci)   
  
  loci_vector <- distinct_loci %>% pull(n_loci)
  
  .loci_list <- c(1:sum(loci_vector))
  
  assigned_loci <- distinct_loci %>% 
    mutate(partitions = partition.vector(
      sample(.loci_list), 
      sep = loci_vector
    ), 
    list_index = row_number())
  
  left_join(.sim_setup, 
            assigned_loci, 
            by = c("n_loci")) %>% 
    return()
}

read_genotype_matrices <- function(.sim_setup){
  distinct_partitions <- .sim_setup %>% 
    distinct(n_loci, partitions)
  
  partition_list <- distinct_partitions %>% pull(partitions)
  
  l_matrix <- vector("list", 
                     length = nrow(distinct_partitions))
  l_varfile <- vector("list", 
                       length = nrow(distinct_partitions))
  l_colnames <- vector("list", 
                       length = nrow(distinct_partitions))
  
  for(i in 1:length(l_matrix)){
    l_read <- vector("list", length = length(partition_list[[i]]))
    l_names <- vector("list", length = length(partition_list[[i]]))
    l_vars <- vector("list", length = length(partition_list[[i]]))
    for(j in 1:length(l_read)){
      print(str_c("Reading from locus num ", partition_list[[i]][[j]]))
      #Read in LA matrix
      l_read[[j]] <- read_tsv(str_c("../../data/simulated_genotypes/AFR_10656_EUR_6624/locus", partition_list[[i]][[j]],
                                    "/locus", partition_list[[i]][[j]], "_admixed_r2", r2,".la_dosage.gz"), 
                              col_names = F, 
                              col_types = cols(
                                .default = col_double()
                              )) %>% 
        as.matrix() 
      
      #Get column names from var file. 
      l_vars[[j]] <- read_tsv(str_c("../../data/simulated_genotypes/AFR_10656_EUR_6624/locus", partition_list[[i]][[j]],
                                "/locus", partition_list[[i]][[j]], "_admixed_r2", r2, ".var"), 
                              col_types = cols(
                                var_id = col_character(),
                                locus_number = col_double(),
                                position = col_double(),
                                minor_allele = col_double(),
                                afr_col_num = col_double(),
                                maf = col_double(),
                                afr_maf = col_double(),
                                eur_maf = col_double()
                              ))
      
      var_ids <- l_vars[[j]] %>% pull(var_id)
      l_names[[j]] <- rbind(
        str_c(var_ids, "_AFR"), 
        str_c(var_ids, "_EUR")
      ) %>% c() 
      
    }
    
    l_matrix[[i]] <- custom_scale_matrix(do.call(cbind, l_read))
    #add in allele categories for sampling
    l_varfile[[i]] <- bind_rows(l_vars) %>% 
      rowwise() %>% 
      mutate(af_cat = case_when(
        afr_maf >= 0.05 & eur_maf >= 0.05 ~ "common in both", 
        afr_maf >= 0.05 & eur_maf < 0.05 ~ "uncommon in eur", 
        afr_maf < 0.05 & eur_maf >= 0.05 ~ "uncommon in afr", 
        afr_maf < 0.05 & eur_maf < 0.05 ~ "rare in both"
      )) %>% 
      ungroup()
    l_colnames[[i]] <- do.call(c, l_names)
    
    rm(l_read)
  }
  return(list(matrices = l_matrix, 
              var_file = l_varfile, 
              col_names = l_colnames))
}

sample_shared_betas <- function(total_snps, n_shared_snps,
                                shared_beta_indices, 
                                shared_beta_sd, 
                                list_index){
  
  shared_beta_vec <- vector("double", 
                            length = 2 * total_snps)
  shared_betas <- rnorm(n = n_shared_snps, mean = 0, sd = shared_beta_sd)
  
  index <- c(2 * shared_beta_indices - 1, 2 * shared_beta_indices)
  
  shared_beta_vec[index] <- rep(shared_betas, 2)
  names(shared_beta_vec) <- l_genotype_matrices$col_names[[1]]
  return(shared_beta_vec)
}

sample_specific_betas <- function(total_snps, n_diff_snps,
                                  diff_beta_indices, 
                                  specific_beta_sd){
  specific_beta_vec <- vector("double", 
                              length = 2 * total_snps)
  specific_betas <- rnorm(n = 2 * n_diff_snps, mean = 0, sd = specific_beta_sd)
  
  index <- c(2 * diff_beta_indices - 1, 2 * diff_beta_indices)
  
  specific_beta_vec[index] <- specific_betas
  return(specific_beta_vec)
}

pheno_num_key <- tibble(
  pheno_num = c(1,2,3), 
  af_cat = c("common in both", 
             "uncommon in afr", 
             "uncommon in eur")
) 

sample_causal_indicies_by_af_cat <- function(.pheno_num, 
                                             .index, .p_causal){
  
  inner_join(
    pheno_num_key %>% 
      filter(pheno_num == .pheno_num),
    l_genotype_matrices$var_file[[.index]] %>% 
      mutate(snp_num = row_number()),
    by = "af_cat"
  ) %>% 
    slice_sample(n = max(ceiling(.p_causal * nrow(.)), 2)) %>%
    pull(snp_num)
}


sample_causal_indices_by_diff_cat <- function(.n_causal_snps,
                                              .prop_small_diff_snps, 
                                              .index, 
                                              .var_file){
  bind_rows(
    .var_file %>% 
      mutate(snp_num = row_number()) %>% 
      filter(str_detect(diff_cat, "\\[0,")) %>% 
      slice_sample(n = .n_causal_snps * .prop_small_diff_snps),
    .var_file %>% 
      mutate(snp_num = row_number()) %>% 
      filter(str_detect(diff_cat, ",1]")) %>% 
      slice_sample(n = .n_causal_snps * (1-.prop_small_diff_snps))
  ) %>% 
    pull(snp_num)
}


sample_causal_indices <- function(.n_causal_snps, .index){
  l_genotype_matrices$var_file[[.index]] %>% 
    mutate(snp_num = row_number()) %>% 
    group_by(af_cat) %>% 
    slice_sample(n = .n_causal_snps/3) %>% 
    pull(snp_num)
}

get_causal_betas <- function(diff_betas, 
                             shared_betas, 
                             total_snps, 
                             .index){
  
  causal_betas <- vector("double", 
                         length = 2 * total_snps)
  
  i_shared <- which(shared_betas != 0)
  i_diff <- which(diff_betas != 0)
  
  causal_betas[i_shared] <- shared_betas[i_shared]
  causal_betas[i_diff] <- diff_betas[i_diff]
  
  names(causal_betas) <- l_genotype_matrices$col_names[[.index]]
  
  causal_betas
}


write_causal_genotypes <- function(n_loci, list_index, 
                                   out_root){

  out_matrix <- str_c(out_root, "/COSI_causal_genotypes_nloci_", n_loci, "_r2", r2, ".tsv.gz")
  out_col <- str_c(out_root, "/COSI_causal_genotypes_nloci_", n_loci, "_r2", r2, ".colnames")
  
  #Output column names
  write_tsv(as.data.frame(l_genotype_matrices$col_names[[list_index]]), 
            file = out_col, 
            col_names = F)
  
  #Output matrix 
  write_tsv(as.data.frame(l_genotype_matrices$matrices[[list_index]]),
          file = out_matrix, 
          col_names = F)
}



write_nested_lists <- function(n_loci, list, list_name, out_root){
  write_tsv(file = str_c(out_root, "/COSI_phenos_nloci_", 
                         n_loci, "/COSI_sims_nloci_", n_loci, "_", list_name), 
            x = as.data.frame(list), 
            col_names = F)
}


#Simulation parameters to form a grid
heritability <- c(0.2,0.6)
pheno_num <- c(1, 2, 3)
p_causal <- c(1, .5, .05)
shared_prop <- c(1, 0.8, 0.5)
n_loci <- 1
rep <- c(1:10) # change this in the future if I want repeated runs per set up. 

#Get subject IDs from one file because they're the same across all simulations. 
subject_ids <- read_tsv("../../data/simulated_genotypes/AFR_10656_EUR_6624/locus1/locus1_admixed.id") %>% 
  distinct(sample_id) %>% pull()

set.seed(13)
sim_setup <- expand_grid(
  heritability, 
  pheno_num,
  p_causal,
  shared_prop, 
  n_loci,
  rep
) %>%
  assign_loci()

l_genotype_matrices <- read_genotype_matrices(sim_setup)

#Used to loop through the locus level data outside R, for example, with PLINK. 
locus_key <- sim_setup %>% 
  distinct(n_loci, rep, partitions) %>% 
  unnest(c(partitions)) %>% 
  rename(locus_num = partitions)

set.seed(13)
sims <- sim_setup %>% 
  rowwise() %>% 
  mutate(n_snps = ncol(l_genotype_matrices$matrices[[list_index]])/2,
         n_subjects = nrow(l_genotype_matrices$matrices[[list_index]]),
         causal_beta_indices = list(sample_causal_indicies_by_af_cat(pheno_num, 
                                                                     .index = list_index, 
                                                                     p_causal)),
         n_causal_snps = length(causal_beta_indices), 
         shared_beta_indices = list(sample(causal_beta_indices, 
                                           size = floor(
                                             length(causal_beta_indices) * 
                                               shared_prop
                                           ))), 
         diff_beta_indices = list(causal_beta_indices[which(!(causal_beta_indices %in% shared_beta_indices))]), 
         shared_betas = list(sample_shared_betas(total_snps = n_snps,
                                      n_shared_snps = length(shared_beta_indices), 
                                      shared_beta_indices = shared_beta_indices, 
                                      shared_beta_sd = 1, 
                                      list_index = list_index)), 
         diff_betas = list(sample_specific_betas(total_snps = n_snps, 
                                                 n_diff_snps = length(diff_beta_indices), 
                                                diff_beta_indices = diff_beta_indices, 
                                                specific_beta_sd = 1)), 
         causal_betas = list(get_causal_betas(diff_betas = diff_betas, 
                                         shared_betas = shared_betas, 
                                         total_snps = n_snps, 
                                         .index = list_index)), 
         var_xb = var(as.numeric(l_genotype_matrices$matrices[[list_index]] %*% 
                                   causal_betas)),
         env_var = (var_xb * (1-heritability))/heritability,
         env_noise = list(rnorm(n = n_subjects, 
                                mean = 0, 
                                sd = sqrt(env_var))), 
         sim_data = list(tibble(FID = subject_ids, 
                                IID = subject_ids,
                                pheno = as.numeric(l_genotype_matrices$matrices[[list_index]] %*% causal_betas + env_noise)))) %>% 
  ungroup() 


#Output files
#Output the "c-bound" matrices for each unique value of n_loci and rep. 
causal_genotypes_out_root <- sprintf("../../data/simulated_genotypes/causal_genotypes/%s", 
                                     sim_name)
dir.create(path = causal_genotypes_out_root, showWarnings = F)


sims %>% 
  distinct(n_loci, list_index) %>% 
  pmap(.f = write_causal_genotypes, out_root = causal_genotypes_out_root)

#Output a plink pheno file for each unique value of n_loci and rep. 
phenotypes_out_root <- sprintf("../../data/simulated_phenotypes/%s", 
                               sim_name)
dir.create(path = phenotypes_out_root, showWarnings = F)

write_plink_pheno <- function(.n_loci, .rep, out_root){
  out_dir <- str_c(out_root, "/COSI_phenos_nloci_", .n_loci)
  dir.create(out_dir, showWarnings = F)
  out_file <- str_c(out_dir, "/", "COSI_phenos_nloci_", .n_loci, "_rep_", .rep, "_r2", r2, ".pheno")
  
  sims %>% 
    filter(n_loci == .n_loci, 
           rep == .rep) %>% 
    dplyr::select(heritability:n_loci, sim_data) %>% 
    unnest(sim_data) %>% 
    pivot_wider(id_cols = c(FID, IID), 
                names_from = c(heritability, pheno_num, p_causal, shared_prop), 
                names_sep = "_", 
                names_prefix = "pheno_",
                values_from = pheno) %>%
    write_tsv(file = out_file)
}


sims %>% 
  distinct(n_loci, rep)  %>% 
  rename(.n_loci = n_loci, 
         .rep = rep) %>% 
  pmap(.f = write_plink_pheno, out_root = phenotypes_out_root)

set.seed(2)
train_test_splits <- sims %>% 
  distinct(n_loci) %>% 
  rowwise() %>% 
  mutate(ids = list(subject_ids), 
         tt_split = list(partition(ids, 
                                   p = c(train = 0.8, 
                                         test = 0.2), 
                                   type = "basic")), 
         training_ids = list(
           tibble(
             FID = ids[tt_split[[1]]], 
             IID = ids[tt_split[[1]]]
           )), 
         testing_ids = list(
           tibble(FID = ids[tt_split[[2]]], 
                  IID = ids[tt_split[[2]]])
         ))



#Output training IDs. 
pmap(train_test_splits %>% 
       select(n_loci, list = training_ids), .f = write_nested_lists, 
     out_root = phenotypes_out_root, 
     list_name = "training_ids") %>% 
  invisible()

#Output testing IDs. 
pmap(train_test_splits %>% 
       select(n_loci, list = testing_ids), .f = write_nested_lists, 
     out_root = phenotypes_out_root, 
     list_name = "testing_ids") %>% 
  invisible()

#Output the locus key
write_tsv(locus_key, sprintf("../../data/COSI_%s_locus_key", sim_name))

#Save nested data frame as .RDS object. 
saveRDS(sims, file = sprintf("../../data/COSI_%s_results.RDS", sim_name))

#Output table to loop through for bash scripting of simulation parameters
sims %>% 
  select(heritability:rep) %>% 
  write_tsv(sprintf("../../data/COSI_%s_parameters", sim_name))

