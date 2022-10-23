library(tidyverse)

for(r2 in c(0.2, 0.5, 1)){
  for(ancestry in c("AFR", "EUR")){
    sim_name <- sprintf("ideal_FL_v2_r2%s", r2)
    
    sims <- readRDS(sprintf("../../data/COSI_%s_results.RDS", 
                            sim_name)) %>% 
      select(heritability:rep,causal_betas)
    
    ref <- read_delim(sprintf("../../data/simulated_genotypes/AFR_10656_EUR_6624/locus1/%s_ref_LDPruned%s_admixedFiltered.raw", 
                              str_to_lower(ancestry), r2), 
                      delim = ' ')
    
    var_file <- read_tsv(sprintf("../../data/simulated_genotypes/AFR_10656_EUR_6624/locus1/locus1_admixed_r2%s.var", 
                                 r2))
    
    
    #Function to go from ref to ref_matrix
    get_matrix <- function(.ref){
      .ref %>% 
        select(IID, all_of(contains("locus"))) %>% 
        column_to_rownames("IID") %>% 
        as.matrix()
    }
    
    #Function to get the summary statistics to flip. 
    flip_vars <- function(.ref_matrix, .var){
      cols <- colnames(.ref_matrix)
      
      var_file_cols <- var_file %>% 
        mutate(allele_call = if_else(minor_allele == 1, 
                                     "C", "A"), 
               col_name = str_c(var_id, allele_call, sep = "_")) %>% 
        pull(col_name)
      
      which(cols != var_file_cols)
    }
    
    get_pred_betas <- function(.causal_betas, .anc){

      b <- .causal_betas[str_detect(names(.causal_betas), .anc)]
      names(b) <- str_remove(names(b), pattern = sprintf("_%s", .anc))
      
      
      flip <- flip_vars(ref_m, var_file)
      b[flip] <- b[flip] * -1
      b
    }
    
    
    ref_m <- get_matrix(ref)
  
    set.seed(1)
    sim_df <- sims %>% 
      rowwise() %>% 
      mutate(pred_betas = list(get_pred_betas(causal_betas, 
                                              .anc = ancestry)), 
             var_xb = var(as.numeric(ref_m %*% 
                                       pred_betas)),
             env_var = (var_xb * (1-heritability))/heritability,
             env_noise = list(rnorm(n = nrow(ref_m), 
                                    mean = 0, 
                                    sd = sqrt(env_var))), 
             sim_data = list(tibble(FID = rownames(ref_m), 
                                    IID = rownames(ref_m),
                                    pheno = as.numeric( as.numeric(ref_m %*% pred_betas)+ env_noise))))
    
    
    #Output a plink pheno file for each unique value of n_loci and rep. 
    phenotypes_out_root <- sprintf("../../data/simulated_phenotypes/%s/ref_%s", 
                                   sim_name, ancestry)
    dir.create(path = phenotypes_out_root, showWarnings = F)
    
    write_plink_pheno <- function(.n_loci, .rep, out_root){
      out_dir <- str_c(out_root, "/COSI_phenos_nloci_", .n_loci)
      dir.create(out_dir, showWarnings = F)
      out_file <- str_c(out_dir, "/", "COSI_phenos_nloci_", .n_loci, "_rep_", .rep, "_r2", r2, ".pheno")
      
      sim_df %>% 
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
    
    
    sim_df %>% 
      distinct(n_loci, rep)  %>% 
      rename(.n_loci = n_loci, 
             .rep = rep) %>% 
      pmap(.f = write_plink_pheno, out_root = phenotypes_out_root)
  }
}

