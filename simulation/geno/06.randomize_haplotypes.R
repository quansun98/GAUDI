library(tidyverse)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser,
                     opt_str = c("-o", "--out-dir"))
parser <- add_option(parser,
                     opt_str = c("-l", "--locus-num"))
parser <- add_option(parser,
                     opt_str = c("-s", "--seed"))
opt <- parse_args(parser)
opt

#opt <- list(
#  o = "../../data/simulated_genotypes/AFR_10656_EUR_6624/locus1/",
#  l = 1,
#  s = 1, 
#  n = "AFR_10656_EUR_6624"
#)

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


#Read in data
raw_eur_haps <- read_delim(str_c(opt$o, "/sim_haplotypes_eur"), 
                           delim = " ", 
                           trim_ws = T, 
                           col_types = cols(
                             .default = col_double(),
                             chr_id = col_character(),
                             left_ancestry = col_character(),
                             switch_point_pos = col_logical()
                           ))

raw_afr_haps <- read_delim(str_c(opt$o, "/sim_haplotypes_afr"), 
                           delim = " ", 
                           trim_ws = T, 
                           col_types = cols(
                             .default = col_double(),
                             chr_id = col_character(),
                             left_ancestry = col_character(),
                             switch_point_pos = col_logical()
                           )) 

raw_aa_haps <- read_delim(str_c(opt$o, "/sim_haplotypes_aa"),
                          delim = " ", 
                          col_types = cols(
                            .default = col_double(),
                            chr_id = col_character(),
                            left_ancestry = col_character()
                          )) 


joint_raw_haps <- bind_rows(raw_aa_haps, raw_eur_haps, 
                            raw_afr_haps)

#get chr ids
chr_ids <- joint_raw_haps %>% 
  distinct(chr_id) %>% 
  pull()

#Simulation IDs. Put each in twice. Then randomize chromosome IDs.
#Each row is a haplotype. 
set.seed(opt$s)
sims <- tibble(
  sample_id = str_c("SIM_", sort(rep(1:(nrow(joint_raw_haps)/2), 2)))
) %>% 
  mutate(chr_id = sample(chr_ids))

# Only need this chunk of code if we're doing LD pruning on the admixed samples, which is probably not a good idea. 
#
# sim_haps_ped <- left_join(sims, joint_raw_haps, 
#           by = "chr_id") %>% 
#   haps_to_plink_ped()
# 
# write_tsv(sim_haps_ped, sprintf("%s/sim_haps_all_subjects.ped", opt$o), 
#           col_names = F)

left_join(sims, joint_raw_haps, 
          by = "chr_id") %>% 
  write_tsv(sprintf("%s/sim_haps_all_subjects.tsv", opt$o))

### End Script here. 









