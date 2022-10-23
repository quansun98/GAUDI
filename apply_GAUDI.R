library(tidyverse)
library(data.table)
library(genlasso)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, opt_str = c("--gaudi-path"),
		     help = "GAUDI path, where cv_fused_lasso.R is located")
parser <- add_option(parser, opt_str = c("--model"), 
		     help = "GAUDI model, *.best_list.RDS")
parser <- add_option(parser, opt_str = c("--target-la-dir"),
		     help = "Directory containing .la_dosage.tsv.gz files for target set (could be chunk-separated)")
parser <- add_option(parser, opt_str = c("--out"),
		     help = "Output file name")
opt <- parse_args(parser)

source(paste0(opt$`gaudi-path`,"/cv_fused_lasso.R"))

fits <- readRDS(opt$model)

apply_prs <- function(.chunk_file, .prs_snps, 
                      .prs_weights, 
                      .fit_vec){
  
  m <- fread(.chunk_file, header = T)  %>% 
    column_to_rownames("SNP") %>% 
    as.matrix() %>% t()

  #Subset M on PRS SNPS
  snps_in_prs <- which(colnames(m) %in% .prs_snps)
  
  
  
  if(length(snps_in_prs) > 1) { 
    pred_m <- m[,which(colnames(m) %in% .prs_snps)]
  }
  
  else if(length(snps_in_prs) == 1) { 
    pred_m <- as.matrix(m[,which(colnames(m) %in% .prs_snps)], ncol = 1)
    colnames(pred_m) <- colnames(m)[which(colnames(m) %in% .prs_snps)]
  }
  
  else{
    print(sprintf("No fitted snps in partial PRS for file %s", basename(.chunk_file)))
    return(.fit_vec)
  }
  
  names(.prs_weights) <- .prs_snps
  update <- pred_m %*% .prs_weights[colnames(pred_m)]
  .fit_vec <- .fit_vec + update
  rownames(.fit_vec) <- rownames(pred_m)
  return(.fit_vec)
}

testing_chunk_files <- list.files(opt$`target-la-dir`, pattern = "*la_dosage.tsv.gz", full.names = T)

for(i in 1:length(testing_chunk_files)){
  
  if(i == 1){
    test_prs_fits = 0
  }
  
  test_prs_fits <- apply_prs(testing_chunk_files[[i]], 
            fits$fit_model$snps, 
            .prs_weights = as.matrix(get_cv_fl_best_betas(fits$fit_model), ncol = 1), 
            .fit_vec = test_prs_fits)
}


test_prs_fits_out = data.frame(ID = rownames(test_prs_fits), prs = test_prs_fits)
write_tsv(test_prs_fits_out, opt$out)


