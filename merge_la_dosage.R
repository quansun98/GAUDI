library(tidyverse)
library(data.table)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)

# args[1]: directory containing la_dosage.tsv.gz files 
#all files under this directory with suffix ".la_dosage.tsv.gz" will be considered
#the combination is separately for each chromosome

la_files <- list.files(args[1], 
                       full.names = T, 
		       pattern = "*la_dosage.tsv.gz")

out_name <- str_replace(la_files, "chr[0-9]*_chunk[0-9]*", "allChr") %>% 
  str_remove(".tsv.gz") %>% unique()

print(length(la_files))

l <- vector("list", length = length(la_files))

for(i in 1:length(la_files)){
 
  print(sprintf("Reading file %s", basename(la_files[[i]])))
  l[[i]] <- fread(la_files[[i]], sep = "\t", 
                  header = T) %>% as_tibble()
  
}

s_m_list <- vector("list", length = 22)
for(i in 1:22){
  print(i)
  m <- bind_rows(l[which(str_detect(la_files, pattern = sprintf("chr%s_", i)))]) %>%
    column_to_rownames("SNP")  %>%
    as.matrix()	
  
  s_m_list[[i]] <- as(m, "dgCMatrix")
  rm(m)
}

s_m <- do.call(rbind, s_m_list)

writeMM(s_m, paste0(out_name, ".mtx"))
#Compress file
system(sprintf("gzip %s",paste0(out_name, ".mtx")))
write_tsv(tibble(c = colnames(s_m)), 
            paste0(out_name, ".colnames"), 
          col_names = F)
write_tsv(tibble(r = rownames(s_m)), 
            paste0(out_name, ".rownames"), 
            col_names = F)



