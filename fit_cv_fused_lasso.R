suppressMessages(library(tidyverse)) # Needed for pipe and other dplyr fxns.
library(Matrix)%>% suppressMessages() #Needed for sparseMatrix class
library(optparse)%>% suppressMessages() # Needed for command line ops
library(data.table)%>% suppressMessages() #Needed for fread
library(genlasso)%>% suppressMessages() #Needed for model fitting. 
library(splitTools) %>% suppressMessages() # Needed for cv function
library(caret) #Needed for the findCorrelation funciton
library(qlcMatrix) #Needed for the corSparse function.

#Creates option object
parser <- OptionParser()
#Path to the software, this is useful when it's run from not the source directory.
parser <- add_option(parser,
			opt_str = c("--gaudi-path"))

parser <- add_option(parser,
                     opt_str = c("--gwas"),
		     help = "A file path to either a set of GWAS results by choromosome or a single file.
			These files are needed for cross-validation of the p-value threshold.
			If multiple files, chromosome number should be replaced in the file path with #.
			We assume that the file has at least two columns: the variant ID given by option --gwas-col-id,
			and and the P-value given by option --gwas-col-p")

#Useful to be flexible with different software generating GWAS summary statistics.
parser <- add_option(parser,
		     opt_str = c("--gwas-col-id"),
		     help = "Column name for variant ID in the GWAS files.")
parser <- add_option(parser,
		     opt_str = c("--gwas-col-p"),
		     help = "Column name for the GWAS p-value in the GWAS files.")
#True or false as to whether the gwas p-value is on the -log10 scale or not.
parser <- add_option(parser,
			opt_str = c("--gwas-log-p"),
			action = "store_true",
			default = F)

#Options to accept the local-ancestry aware dosage files. Output of the python script.
parser <- add_option(parser,
                     opt_str = c("--la"))
parser <- add_option(parser,
                     opt_str = c("--col"))
parser <- add_option(parser,
                     opt_str = c("--row"))


parser <- add_option(parser,
                     opt_str = c("--train"))
parser <- add_option(parser,
                     opt_str = c("--pheno"))
parser <- add_option(parser,
                     opt_str = c("--pheno-name"))
parser <- add_option(parser,
                     opt_str = c("--pheno-fid"))
parser <- add_option(parser,
                     opt_str = c("--pheno-iid"))


#Starting and ending exponents of the p-value grid.
parser <- add_option(parser,
                     opt_str = c("--start-p-exp"),
                     default = -4)
parser <- add_option(parser,
                     opt_str = c("--end-p-exp"),
                     default = -8)

#Random seed
parser <- add_option(parser,
                     opt_str = c("--seed"))

#Whether or not to include a sparsity component for the GAUDI model.
parser <- add_option(parser,
                     opt_str = c("--sparsity"),
                     action = "store_true",
                     default = FALSE)

#Out path name.
parser <- add_option(parser,
                     opt_str = c("--out"))

#If true, the model is more verbose with the fitting output.
parser <- add_option(parser,
                     opt_str = c("--verbose",
                                 "-v"))

#Parameters for the gamma parameter grid.
parser <- add_option(parser, opt_str = "--gamma-start")
parser <- add_option(parser, opt_str = "--gamma-stop")
parser <- add_option(parser, opt_str = "--gamma-by")
opt <- parse_args(parser)

#Get GAUDI fitting file.
source(sprintf("%s/cv_fused_lasso.R", opt$`gaudi-path`))

# opt <- list(
#   gwas = "/proj/yunligrp/users/bryce/prs/admixed_prs/data/real_data_analysis/regenie_results/AA_adjusted_phenotypes_WHI_MEGA_AA_#_trainingSamples_fold1_regenie_hematocrit_adj.regenie",
#   `gwas-col-p` = "LOG10P",
#   `gwas-col-id` = "ID",
#   `gwas-log-p` = TRUE,
#   la = "/pine/scr/b/r/bryce38/gaudi/real_data_analysis/apply_prs_methods/gaudi/la_dosage_files/AA_WHI_MEGA_AA_fold1.allChr_hematocrit_adj_p0.05.la_dosage.mtx.gz",
#   col = "/pine/scr/b/r/bryce38/gaudi/real_data_analysis/apply_prs_methods/gaudi/la_dosage_files/AA_WHI_MEGA_AA_fold1.allChr_hematocrit_adj_p0.05.la_dosage.colnames",
#   row = "/pine/scr/b/r/bryce38/gaudi/real_data_analysis/apply_prs_methods/gaudi/la_dosage_files/AA_WHI_MEGA_AA_fold1.allChr_hematocrit_adj_p0.05.la_dosage.rownames",
#   train = "/proj/yunligrp/users/bryce/prs/admixed_prs/data/real_data_analysis/phenotype_data/training_samples_lists/AA_training_samples_WHI_MEGA_AA_fold1.tsv",
#   pheno = "/proj/yunligrp/users/bryce/prs/admixed_prs/data/real_data_analysis/phenotype_data/adjusted_phenos_by_cohort/AA_adjusted_phenotypes_WHI_MEGA_AA.tsv.gz",
#   `pheno-name` = "hematocrit_adj",
#   seed = 13,
#   out = "hem_test",
#   sparsity = F,
#   verbose = F,
#   `start-p-exp` = -5,
#   `end-p-exp` = -50
# )
# source("cv_fused_lasso.R")

#Reads in GWAS data and filters based on p-value cutoff.
read_and_filter <- function(path, p_name, filter_p, 
                            log_p){
  
  if(log_p){
    fread(path) %>% 
      as_tibble() %>% 
      mutate(P = 10^(-.data[[p_name]])) %>% 
      filter(P < filter_p)
  } else {
    fread(path) %>% 
      as_tibble() %>% 
      rename(P = all_of(p_name))%>% 
      filter(P < filter_p)
  }
  
  
}

#Parse a given GWAS file, or set of files.
parse_gwas <- function(gwas_path, 
                       p_name = "P",
                       id_name = "ID",
                       log_p = F, 
                       filter_p  = 5e-4){
  
  if(str_detect(gwas_path, pattern = "#")){
    
    p <- str_replace(basename(gwas_path), pattern = "#", replacement = "[0-9]*") %>% 
      str_c("$")
    
    f <- list.files(path = dirname(gwas_path), 
                    pattern = p, 
                    full.names = T)
    print(f)
    print(length(f))
    
    l <- vector("list", length = length(f))
    
    for(i in 1:length(f)){
      print(basename(f[[i]]))
      l[[i]] <- read_and_filter(f[[i]], p_name, 
                                filter_p, 
                                log_p)
      
    }
    .gwas <- bind_rows(l)
  }
  else{
    .gwas <- read_and_filter(gwas_path, p_name, filter_p, 
                             log_p)
  }
  .gwas
}

#Removes columns from the la matrix that are constant. This causes problems with the model fitting if not.
filter_fixed_cols <- function(la){
  
  passing_vars <- apply(la, 2, sd) != 0
  
  la[,passing_vars]
  
}

#Despite doing LD pruning, correlation may reappear in rare variants in the local ancestry dosage files. This does a cosine similarity measure to remove highly similar columns in the sparse matrix while maintaining sparsity (so it's still fast). 
prune_matrix <- function(cor_la, 
                         cor_thresh = .95){
  
  if("sparseMatrix" %in% is(cor_la)){
    print("Sparse Matrix found; using cosine similarity measure to remove highly similar columns.")
    cor_index <- findCorrelation(as.matrix(cosSparse(cor_la)), cutoff = cor_thresh, 
                             verbose = T)
  }
  else {
    cor_index <- findCorrelation(cor(cor_la), cutoff = cor_thresh, 
                                 verbose = T)
  }
  if(length(cor_index) > 0){
    cor_la[,-cor_index]
  } else{
    cor_la
  }
}



#Read in la_matrix, which is stored as a sparse matrix. We read in reference files for column and row names. 
la_matrix <- readMM(opt$la)

col_names <- read_tsv(opt$col,
                      col_names = F) %>% pull()

row_names <- read_tsv(opt$row, col_names = F) %>%
  distinct() %>% pull()

if(length(col_names) != ncol(la_matrix)){
  stop("Column names do not match dimensions of the LA matrix.\nCheck files and try again")
}

if(length(row_names) != nrow(la_matrix)){
  stop("Row names do not match dimensions of the LA matrix.\nCheck files and try again")
}

colnames(la_matrix) <- col_names
rownames(la_matrix) <- row_names

#Then transpose for subjects in rows and variants in columns.
la_matrix <- t(la_matrix)

#Read in pheno data
pheno_string <- opt$`pheno-name`

pheno <- read_tsv(opt$pheno) %>%
  select(c(FID = opt$`pheno-fid`, all_of(pheno_string))) %>% 
  mutate(FID = as.character(FID)) %>%
  inner_join(tibble(FID = rownames(la_matrix)), .) %>% 
  filter(complete.cases(.))

#Read in all IDs from training data
training_ids <- read_tsv(opt$train) %>% pull(IID)


#Subset pheno on training ids
# pheno_fold_train <- 
pheno_train <- pheno %>% filter(FID %in% training_ids) 

pheno_train_ids <- pheno_train %>% pull(FID)
pheno_train_m <- pheno_train %>% pull(-FID) %>% as.matrix()


#Initialize p-grid.
p_list <- 5 * 10^(opt$`start-p-exp`:opt$`end-p-exp`)
p_fits <- vector("list", length = length(p_list))
r2_vec <- rep(0, length(p_list))


#Read in GWAS results
print("Reading in GWAS file.")
gwas <- parse_gwas(opt$gwas, 
                   p_name = opt$`gwas-col-p`, 
                   id_name = opt$`gwas-col-id`, 
                   log_p = opt$`gwas-log-p`, 
                   filter_p = max(p_list))
print("Done.")


#Remove suffix from column names for "which" matching.
cols_no_suffix <- str_remove(colnames(la_matrix), "_AFR|_EUR")


#Subset la_matrix on training ids and maximum variants.
la_matrix_train <- la_matrix[which(rownames(la_matrix) %in% pheno_train_ids),
                             which(cols_no_suffix %in% gwas[[opt$`gwas-col-id`]])]

#Filter variants for MAC threshold in training sample. 
la_matrix_train <- filter_fixed_cols(la_matrix_train)

la_matrix_train <- prune_matrix(la_matrix_train)

#Remove suffix from column names for "which" matching.
cols_no_suffix <- str_remove(colnames(la_matrix_train), "_AFR|_EUR")

# #Clear la_matrix from memory.
# rm(la_matrix)

#for p in p_list
set.seed(13)
# for(i in 1:length(p_list)){
for(i in 1:length(p_list)){
  p <- p_list[i]
  
  sig_vars <- gwas %>% 
      filter(P < p) %>%
      pull(opt$`gwas-col-id`)
  
  #Subset on variants passing significance threshold.
  la_matrix_train_p <- la_matrix_train[,which(cols_no_suffix %in% sig_vars)]
  print(dim(la_matrix_train_p))
  
  #Need more than 2 SNPs to fit fused Lasso with GenLasso penalty matrix.
  if(ncol(la_matrix_train_p) <= 2 || is.null(ncol(la_matrix_train_p))){
    print(sprintf("Not enough variants achieve marginal significance threshold at p=%s.",p))
    p_fits[[i]] <- NA
    r2_vec[[i]] <- NA_real_
    next
  }

  #Fit the CV model
  p_fits[[i]] <- cv_fused_lasso(scale(as.matrix(la_matrix_train_p)),
                                pheno_train_m,
                                n_folds = 5,
                                verbose = T,
                                sparsity = opt$sparsity, 
				gamma_start = opt$`gamma-start`,
				gamma_stop = opt$`gamma-stop`,
				gamma_by = opt$`gamma-by`)
  
  if(!is.na(p_fits[[i]])){
    r2_vec[i] <- p_fits[[i]]$cv_r2
  }else{
    r2_vec[i] <- NA
  }
  
}


best_index <- which(r2_vec == max(r2_vec, na.rm = T))
best_fit <- p_fits[[best_index]]
best_p <- p_list[best_index]

best_list <- list("fit_model" = best_fit, 
                  "best_p" = best_p)

print(best_list$fit_model$cv_r2)
best_list$best_p

# #Output best_list as .Rds
best_list %>% saveRDS(file = str_c(opt$out, ".best_list.RDS"))

