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
                     opt_str = c("-n", "--name"))
opt <- parse_args(parser)

#Function to process the recombination map


std_recomb_map <- function(.recomb){
  .recomb %>% 
    mutate(end_pos = lead(start_pos, default = 500000) - 1) %>%  
    rowwise() %>%  
    mutate(expand = list(start_pos:end_pos)) %>% 
    ungroup()%>% 
    select(expand, prob) %>%  
    unnest(expand) %>% 
    mutate(std_prob = prob/sum(prob)) %>% 
    select(pos = expand, std_prob) %>% 
    return()
}

#and return the 
#position of the admixture switchpoint. 
get_recomb_position <- function(.std_recomb, .n){
  .std_recomb %>% 
    summarise(sample(pos, prob = std_prob, size = .n, 
                     replace = T)) %>% 
    pull() %>% 
    return()
}

#Function to determine left_ancestry
get_left_ancestry <- function(pop1, pop2, 
                               .n){
  sample(c(pop1, pop2), size = .n, replace = T) %>% 
    return()
}

#Function to assemble admixed haplotypes. 
assemble_switch_point_hap <- function(afr_cosi_chr_id, 
                                      eur_cosi_chr_id, 
                                      switch_point, 
                                      left_ancestry){

  aa_eur_hap <- eur_haps %>% 
    filter(chr_id == eur_cosi_chr_id) %>% 
    pull(haplotype)
  
  aa_afr_hap <- afr_haps %>% 
    filter(chr_id == afr_cosi_chr_id) %>% 
    pull(haplotype)
  
  #check that haplotypes are same length
  if(nchar(aa_eur_hap) != nchar(aa_afr_hap)){
    print("Ancestral haplotypes are not the same length, please try again.")
    break
  }
  
  #Get the column number of the greatest position less than or equal
  # to the switch point
  left_ancestry_col_num <- pos_key %>% 
    filter(pos <= switch_point) %>% 
    filter(pos == max(pos)) %>% 
    pull(col_number)
  
  #Get the column number of the smallest position greater than the switch
  #point. 
  right_ancestry_col_num <- pos_key %>% 
    filter(pos > switch_point) %>% 
    filter(pos == min(pos)) %>% 
    pull(col_number)
  
  #If left admixture is AFR, then assemble haplotype with AFR on the left.
  if(left_ancestry == "AFR"){
    #Go from 1 to the switch point. 
    left_haplotype <- str_sub(aa_afr_hap, 1, left_ancestry_col_num)
    #Go from one past the switch point to the end. 
    right_haplotype = str_sub(aa_eur_hap, right_ancestry_col_num, 
                              str_length(aa_eur_hap))
  } 
  #If left admixture is EUR, then assemble haplotype with EUR on the left. 
  else if (left_ancestry == "EUR"){
    left_haplotype <- str_sub(aa_eur_hap, 1, left_ancestry_col_num)
    right_haplotype <- str_sub(aa_afr_hap, right_ancestry_col_num, 
                               str_length(aa_afr_hap))
  }
  switch_point_hap <- str_c(left_haplotype, right_haplotype)
  
  return(switch_point_hap)
}

#Recombination map for locus
recomb <- read_delim(str_c(opt$o, "/recomb_model_locus", opt$l), 
                     delim = "\t", 
                     col_names = c("start_pos", "prob"))
std_recomb <- std_recomb_map(recomb)

#Recombination chromosome IDs
cosi_chr_ids <- read_tsv(str_c(opt$o, "/locus", opt$l, "_COSI_hap_IDs_switchpoint"))

#AFR haplotypes
afr_haps <- read_tsv(str_c(opt$o, "/", opt$n, "_locus", opt$l, ".hap-5"), 
                     col_names = c("chr_id", 
                                   "ancestry_code", 
                                   "haplotype")) %>% 
  mutate(haplotype = haplotype %>% 
           str_replace_all(" ", ""))
head(afr_haps)
#EUR haplotypes
eur_haps <- read_tsv(str_c(opt$o, "/", opt$n, "_locus", opt$l, ".hap-1"), 
                     col_names = c("chr_id", 
                                   "ancestry_code", 
                                   "haplotype")) %>% 
  mutate(haplotype = haplotype %>% 
           str_replace_all(" ", ""))


head(eur_haps)
#Read in one COSI .pos file to convert genetic position to column number
pos_key <- read_tsv(str_c(opt$o, "/", opt$n, "_locus", opt$l, ".pos-1"), 
                    skip = 1, 
                    col_names = c("snp", "chrom", "chrom_pos", 
                                  "allele1", "freq1", "allele2", "freq2")) %>% 
  mutate(col_number = row_number()) %>% 
  select(col_number, pos = chrom_pos)
head(pos_key)

final_col_names <- pos_key %>% 
  mutate(col_name = str_c("locus_", opt$l, 
                          "_pos_", pos)) %>% 
  pull(col_name)


set.seed(opt$s)
switch_point_haps <- cosi_chr_ids %>% 
  mutate(switch_point_pos = get_recomb_position(std_recomb, 
                                          .n = nrow(.)), 
         left_ancestry = get_left_ancestry("AFR", "EUR",
                                             .n = nrow(.))) %>% 
  rowwise() %>% 
  mutate(hap = assemble_switch_point_hap(COSI_hap_ID_AA_haplotypes_AFR_mosaic, 
                                         COSI_hap_ID_AA_haplotypes_EUR_mosaic, 
                                         switch_point = switch_point_pos, 
                                         left_ancestry = left_ancestry)) %>% 
  ungroup() %>% 
  mutate(
    chr_id = str_c("AA_", row_number())
  ) %>% 
  select(chr_id, left_ancestry, switch_point_pos, hap)
         
#First, we have to re-insert spaces back into our long string of haplotypes. 
# Then convert the string with spaces to our desired columns before writing. 
#
# There's probably a better way to do this, but this works for now. 
sim_haplotypes_aa <- switch_point_haps %>% 
  rowwise() %>% 
  mutate(hap_list = str_extract_all(hap, boundary("character")),
         hap = str_c(hap_list, collapse = " ")) %>% 
  ungroup() %>% 
  select(-hap_list) %>% 
  separate(hap, into = final_col_names, sep = " ", 
           convert = T)


write_delim(sim_haplotypes_aa, str_c(opt$o, "/sim_haplotypes_aa"), 
            delim = " ")




