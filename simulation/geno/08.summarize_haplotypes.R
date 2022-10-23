#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser,
                     opt_str = c("-o", "--out-dir"))
parser <- add_option(parser,
                     opt_str = c("-l", "--locus-num"))
parser <- add_option(parser,
                     opt_str = c("-s", "--seed"))
parser <- add_option(parser,
                     opt_str = c("-r", "--r2"))
opt <- parse_args(parser)

# opt <- list(
#   o = "../../data/simulated_genotypes/AFR_10656_EUR_6624/locus1/",
#   l = 1,
#   s = 1,
#   r2 = 0.2
# )


### Functions to transform the long genotype file ####
long_to_id_file <- function(long_df){
  long_df %>% 
    distinct(sample_id, chr_id, left_ancestry, switch_point_pos) %>% 
    rowwise() %>% 
    mutate(id_num = as.numeric(str_split(sample_id, "_")[[1]][2])) %>% 
    arrange(id_num) %>% 
    select(-id_num) %>% 
    return()
}

name_list_col <- function(vec, name){
  l <- list(vec)
  names(l) <- name
  return(l)
}



long_grouped_to_h_matrix <- function(long_df_grouped){
  la_long_filter_grouped_by_id_and_position %>% 
    group_by(sample_id, chr_id)%>% 
    summarise(row_vec = list(c(rbind(hap_minor_afr, hap_minor_eur)))) %>% 
    rowwise() %>% 
    mutate(row_vec = name_list_col(row_vec, sample_id), 
           id_num = as.numeric(str_split(sample_id, "_")[[1]][2])) %>% 
    arrange(id_num) %>% 
    select(-id_num) %>% 
    ungroup() %>% 
    pull(row_vec) %>% 
    do.call(rbind, .) 
}

long_grouped_to_g_matrix <- function(long_df_grouped){
  long_df_grouped %>%  
    summarise(g_afr = sum(minor_afr), 
              g_eur = sum(minor_eur)) %>% 
    select(sample_id, g_afr, g_eur) %>% 
    group_by(sample_id) %>% 
    summarise(row_vec = list(c(rbind(g_afr, g_eur)))) %>% 
    rowwise() %>% 
    mutate(row_vec = name_list_col(row_vec, sample_id), 
           id_num = as.numeric(str_split(sample_id, "_")[[1]][2])) %>% 
    arrange(id_num) %>% 
    select(-id_num) %>% 
    ungroup() %>% 
    pull(row_vec) %>% 
    do.call(rbind, .) %>% 
    return()
}

long_grouped_to_dosage_matrix <- function(long_df_grouped){
  long_df_grouped %>% 
    summarise(dosage = sum(minor_allele_indicator)) %>%
    select(-position) %>% 
    summarise(dosage = list(dosage)) %>% 
    rowwise() %>% 
    mutate(dosage = name_list_col(dosage, sample_id), 
           id_num = as.numeric(str_split(sample_id, "_")[[1]][2])) %>% 
    arrange(id_num) %>% 
    select(-id_num) %>% 
    ungroup() %>% 
    pull(dosage)%>% 
    do.call(rbind, .) %>% 
    return()
}

long_to_var_file <- function(long_df){
  long_df %>% 
    group_by(locus_number, position, minor_allele) %>% 
    summarise(maf = mean(minor_allele_indicator == 1), 
              afr_maf = sum(minor_afr, na.rm = T)/sum(afr_local_ancestry, na.rm = T), 
              eur_maf = sum(minor_eur, na.rm = T)/sum(eur_local_ancestry, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(afr_col_num = 2 * row_number() - 1, 
           var_id = str_c("locus_", locus_number, "_pos_", 
                          position)) %>% 
    relocate(var_id) %>% 
    relocate(afr_col_num, .before = maf) %>% 
    return()
}

long_to_map_file <- function(long_df){
  long_df %>% 
    distinct(locus_number, position) %>% 
    rowwise() %>% 
    mutate(var_id = str_c("locus_", locus_number, "_pos_", 
                          position), 
           cm = 0) %>% 
    select(locus_number, var_id, cm, position) %>% 
    return()
}

long_grouped_to_plink_ped <- function(long_df_grouped) {
  long_df_grouped %>% 
    mutate(
      chr_num = row_number()
    ) %>% 
    ungroup() %>% 
    select(sample_id, chr_num,  locus_number,position, allele) %>% 
    #Make sure the chromosomes are grouped together correctly. 
    arrange(sample_id, position, chr_num) %>% 
    pivot_wider(id_cols = sample_id, 
                names_from = c(locus_number, position, chr_num), 
                names_prefix = "locus_", 
                names_sep = "_" ,
                values_from = allele) %>% 
    mutate(FID = sample_id,
           IID = sample_id,
           father_id = 0, 
           mother_id = 0, 
           sex = 0, 
           phenotype = 0) %>% 
    select(-sample_id) %>% 
    relocate(FID, IID, father_id, mother_id, sex, phenotype)
}
### End of Functions ###
joint_haps <- read_tsv(sprintf("%s/sim_haps_all_subjects.tsv", 
                               opt$o), 
                       col_types = cols(
                         .default = col_double(),
                         sample_id = col_character(),
                         chr_id = col_character(),
                         left_ancestry = col_character()
                       ))

ld_pruned_vars <- read_tsv(sprintf("%s/locus%s_ld_pruned_snplist%s", 
                                   opt$o,
                                   opt$l,
                                   opt$r2)) %>% 
  pull(SNP)

joint_haps_id <- joint_haps %>% 
  select(1:4, all_of(ld_pruned_vars)) %>% 
  pivot_longer(cols = starts_with("locus"), 
               names_to = c("locus_number", 
                            "position"), 
               names_pattern = "locus_(.*)_pos_(.*)", 
               values_to = "allele") %>% 
  mutate(locus_number = parse_number(locus_number), 
         position = parse_number(position))

#Put phased genotype data in long format and determine local
# ancestry by switchpoint and left ancestry.
la_long <- joint_haps_id %>% 
  group_by(locus_number, position) %>% 
  mutate(allele_one_frequency = sum(allele != 2)/n()) %>% 
  ungroup() %>% 
  mutate(minor_allele = if_else(allele_one_frequency < 0.5, 
                                1, 2), 
         minor_allele_indicator = if_else(minor_allele == allele, 
                                          1, 0),
         afr_local_ancestry = case_when(
           is.na(switch_point_pos) & left_ancestry=="AFR" ~ 1, 
           !is.na(switch_point_pos) & left_ancestry=="AFR" & 
             position <= switch_point_pos ~ 1, 
           !is.na(switch_point_pos) & left_ancestry=="EUR" & 
             position > switch_point_pos ~ 1,
           TRUE ~ 0
         ), 
         eur_local_ancestry = case_when(
           is.na(switch_point_pos) & left_ancestry=="EUR" ~ 1, 
           !is.na(switch_point_pos) & left_ancestry=="EUR" & 
             position <= switch_point_pos ~ 1, 
           !is.na(switch_point_pos) & left_ancestry=="AFR" & 
             position > switch_point_pos ~ 1, 
           TRUE ~ 0
         ), 
         minor_afr = minor_allele_indicator * afr_local_ancestry, 
         minor_eur = minor_allele_indicator * eur_local_ancestry, 
         hap_minor_afr = if_else(afr_local_ancestry == 0, 
                                 NA_real_, 
                                 minor_afr),
         hap_minor_eur = if_else(eur_local_ancestry == 0, 
                                 NA_real_, 
                                 minor_eur))


#Generate var file, use it for filtering. 
#Make variant summary file which gives allele frequency, and ancestry specific allele frequencies.
var_file_unfiltered <- long_to_var_file(la_long)


la_long_filtered <- inner_join(var_file_unfiltered  %>% 
                                 filter(afr_maf > 0.0025, eur_maf > 0.0025) %>% 
                                 select(locus_number, position), 
                               la_long, 
                               by = c("locus_number", "position")) 

#Group for later processing. 
la_long_filter_grouped_by_id_and_position <- la_long_filtered %>% 
  group_by(sample_id, position)

#Output files in PLINK format for ease of data processing in GWAS and LD pruning.
# Also, this should give us access to more standard data input for other methods. 
ped_file <- long_grouped_to_plink_ped(la_long_filter_grouped_by_id_and_position)
map_file <- long_to_map_file(la_long_filtered)

write_tsv(ped_file, 
          file = sprintf("%s/locus%s_admixed_r2%s.ped", 
                         opt$o, 
                         opt$l, 
                         opt$r2), 
          col_names = F)

write_tsv(map_file, 
          file = sprintf("%s/locus%s_admixed_r2%s.map", 
                         opt$o, 
                         opt$l, 
                         opt$r2), 
          col_names = F)

#Make dosage matrix which accounts for which ancestral haplotype the risk alleles are on.
g_matrix <- long_grouped_to_g_matrix(la_long_filter_grouped_by_id_and_position)
# Currently, FL doesn't support missing data. 
# h_matrix <- long_grouped_to_h_matrix(la_long_filter_grouped_by_id_and_position)
# 
# write_tsv(as.data.frame(h_matrix), 
#           file = sprintf("%s/locus%s_admixed_r2%s.hap_la_dosage.gz", 
#                          opt$o, 
#                          opt$l, 
#                          opt$r2),
#           col_names = F)


write_tsv(as.data.frame(g_matrix), 
          file = sprintf("%s/locus%s_admixed_r2%s.la_dosage.gz", 
                         opt$o, 
                         opt$l, 
                         opt$r2),
          col_names = F)


#Make file specifying local ancestry and chromosomes for subjects at the locus. 
id_file <- long_to_id_file(la_long)

write_tsv(id_file, 
          str_c(opt$o, "/locus",opt$l,"_admixed_r2", opt$r2, ".id"))

#Write var file.
var_file <- long_to_var_file(la_long_filtered)
write_tsv(var_file, 
          str_c(opt$o, "/locus",opt$l,"_admixed_r2", opt$r2, ".var"))

