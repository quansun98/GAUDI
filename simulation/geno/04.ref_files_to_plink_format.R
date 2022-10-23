#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--cosiSimName"))
parser <- add_option(parser, c("--locusNum"))
opt <- parse_args(parser)

haps_to_plink_ped <- function(.hap){
  .hap %>% 
    mutate(across(starts_with("locus"), 
                  ~case_when(
                    #ref allele is always A, alt allele is always C. 
                    . == 1 ~ "C", 
                    . == 2 ~ "A"
                  ))) %>% 
    group_by(sample_id) %>% 
    mutate(chr_num = row_number()) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = sample_id, 
                names_from = chr_num, 
                values_from = starts_with("locus")) %>% 
    mutate(FID = sample_id,
           IID = sample_id,
           father_id = 0, 
           mother_id = 0, 
           sex = 0, 
           phenotype = 0) %>% 
    select(-sample_id) %>% 
    relocate(FID, IID, father_id, mother_id, sex, phenotype)
}

afr_haps <- read_delim(sprintf("../../data/simulated_genotypes/%s/locus%s/sim_haplotypes_ref_afr", 
                               opt$cosiSimName,
                               opt$locusNum), 
                       delim = " ", 
                       trim_ws = T,
                       col_types = cols(
                         .default = col_double(),
                         chr_id = col_character(),
                         left_ancestry = col_character(),
                         switch_point_pos = col_logical()
                       ))
#Generate Ref IDs and shuffle them. 
afr_sample_ids <- rep(paste0("AFR_REF_", 1:(nrow(afr_haps)/2)), 2) %>% sample()


afr_ped <- afr_haps %>% 
  mutate(sample_id = afr_sample_ids) %>% 
  haps_to_plink_ped()

write_tsv(afr_ped, sprintf("../../data/simulated_genotypes/%s/locus%s/afr_ref.ped", 
                           opt$cosiSimName,
                           opt$locusNum), 
          col_names = F)

#EUR
eur_haps <- read_delim(sprintf("../../data/simulated_genotypes/%s/locus%s/sim_haplotypes_ref_eur", 
                               opt$cosiSimName,
                               opt$locusNum), 
                       delim = " ", 
                       trim_ws = T,
                       col_types = cols(
                         .default = col_double(),
                         chr_id = col_character(),
                         left_ancestry = col_character(),
                         switch_point_pos = col_logical()
                       ))
#Generate Ref IDs and shuffle them. 
eur_sample_ids <- rep(paste0("EUR_REF_", 1:(nrow(eur_haps)/2)), 2) %>% sample()


eur_ped <- eur_haps %>% 
  mutate(sample_id = eur_sample_ids) %>% 
  haps_to_plink_ped()


write_tsv(eur_ped, sprintf("../../data/simulated_genotypes/%s/locus%s/eur_ref.ped", 
                           opt$cosiSimName,
                           opt$locusNum), 
          col_names = F)





