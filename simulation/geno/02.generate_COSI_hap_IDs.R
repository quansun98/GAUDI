#!/usr/bin/env Rscript

library(optparse)
suppressMessages(library(tidyverse))

#Compile option list for command line evaluation. 
parser <- OptionParser()
parser <- add_option(parser, 
                     opt_str = c("-o", "--out"), 
                     help = "Out prefix of COSI chromosome
              ID lists.")
parser <- add_option(parser, 
                     opt_str = c("-s", "--seed"),
                     type = "integer", 
                     help = "Random seed.")
opt<-parse_args(parser)

set.seed(opt$seed)

subj_ids_eur <- 0:6623
subj_ids_afr <- 6624:17279

eur_ref_ids <- sample(x = subj_ids_eur, size = 5000, replace = F)
subj_ids_eur_noRef <- subj_ids_eur[!(subj_ids_eur %in% eur_ref_ids)]

ref_ids_eur <- tibble(
  COSI_hap_ID_EUR_ref_haplotypes = eur_ref_ids
)

afr_ref_ids <- sample(x = subj_ids_afr, size = 5000, replace = F)
subj_ids_afr_noRef <- subj_ids_afr[!(subj_ids_afr %in% afr_ref_ids)]

ref_ids_afr <- tibble(
  COSI_hap_ID_AFR_ref_haplotypes = afr_ref_ids
)

#Simulation parameters, can be changed based on the number of recombination
# events assumed in the genomic window. 
total <- 7000 
n_recomb <- total * 0.04 #expected proportion of recombination events in 500kb window. 
n_aa_afr_chr <- (total - n_recomb) * 0.8
n_aa_eur_chr <- total - n_aa_afr_chr - n_recomb

#Sample 1 European and 1 African chromosome for each recombination event. 
# 
# First column is European ID, second column is AFR ID
recomb_ids <- tibble(
  COSI_hap_ID_AA_haplotypes_EUR_mosaic = sample(x = subj_ids_eur_noRef, size = n_recomb),
  COSI_hap_ID_AA_haplotypes_AFR_mosaic = sample(x = subj_ids_afr_noRef, size = n_recomb)
  )

#Vector the matrix for housekeeping. 
subj_ids_recomb <- unlist(recomb_ids)
length(subj_ids_recomb)

subj_ids_eur_noRecomb <- subj_ids_eur_noRef[!(subj_ids_eur_noRef %in% subj_ids_recomb)]
subj_ids_afr_noRecomb <- subj_ids_afr_noRef[!(subj_ids_afr_noRef %in% subj_ids_recomb)]


#Get ID's of chromosomes that we will use as AA chromosomes, but do not
# experience a recombination event. Used for indexing the COSI haplotype
# files
subj_ids_aa_eur <- tibble(
  COSI_hap_ID_EUR_haplotypes = subj_ids_eur_noRecomb
  )
subj_ids_aa_afr <- tibble(
  COSI_hap_ID_AFR_haplotypes = subj_ids_afr_noRecomb
)


#Write out files

write_tsv(recomb_ids, str_c(opt$o, "_switchpoint"))
write_tsv(subj_ids_aa_eur, str_c(opt$o, "_eur"))
write_tsv(subj_ids_aa_afr, str_c(opt$o, "_afr"))
write_tsv(ref_ids_eur, str_c(opt$o, "_ref_eur"))
write_tsv(ref_ids_afr, str_c(opt$o, "_ref_afr"))


