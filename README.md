# GAUDI

## Introduction

GAUDI is a novel polygenic risk score (PRS) method designed specifically for admixed individuals.
It takes the mosaic structure of the admixed genomes into consideration, leveraging local ancestry information to jointly estimate ancestry-shared effects and allowing ancestry-specific effects.
GAUDI model is based on a modified fussed lasso framework, balancing between fusion and sparsity.
The sparsity penalty performs variant selection, and the fusion penalty encourages similarity of the effect sizes of the same variant between two populations.
See our manuscript for more details: [GAUDI preprint](<https://www.biorxiv.org/content/10.1101/2022.10.06.511219v1.abstract>)

This folder also contains codes used in the simulation studies. See `simulation` for more details.


## Instructions on model training

GAUDI requires individual level data for model training and testing.
The below steps assume one already have the individual genotypes (or imputed haplotype dosages [HDS] from minimac4) and local ancestry information (by RFMix) for training data.
For both performance and computation considerations, we recommend performing stringent QC on the imputed dataset (e.g. imputation quality > 0.8)

1. GWAS

GWAS could be performed directly on the training samples, or seeking external GWAS results. Note that GWAS results are only used as a screening step for input variant set.
Only p-value information will be utilized. Effect sizes for variants and LD information will not be used from the GWAS results.
To reduce computational burden, we suggest keeping only variants with GWAS p-value < 0.05 for further variant selection and weight estimation.

2. Combine HDS (or GT) from VCFs with local ancestry information

The script `py_vcf_to_la.py` extracts haplotype dosages (HDS) or genotypes (GT) from VCFs to incorporate local ancestry information. 
Local ancestry inferred by RFMix is required. The below codes provide an example to run (by chunk to require less memory):

		while chunk; do
			start_pos=$(( ($chunk - 1)*10*1000000 ))
			end_pos=$(( ($chunk)*10*1000000 ))
			./py_vcf_to_la.py \
				--local-ancestry data/local_anc_chr22.msp.tsv.gz \
				--vcf data/test_chr22.vcf.gz \
				--include data/chr22_minMAF0.001_Rsq0.9_GWASp0.05.txt \
				--la-dosage-threshold 10 \
				--chr 22 --pos-start $start_pos --pos-stop $end_pos \
				--out data/test_out_chr22_chunk${chunk}
		done < <( awk -v chr=$chr '$1==chr { print $2 }' data/chunk_list

3. Combine files created in the last step

The script `merge_la_dosage.R` combines the by chunk `.la_dosage.tsv.gz` files together. 
All files under the input directory with suffix `.la_dosage.tsv.gz` will be considered to be combined.
The combine is separately for each chromosome. Just run

		Rscript merge_la_dosage.R YOUR_OUT_DIRECTORY

4. Run GAUDI for model training

The main function to run GAUDI is `fit_cv_fused_lasso.R`. 
It provides many customized options. Please see the script for details.
It also calls `cv_fused_lasso.R` internally.
An example code to train GAUDI model is provided below:


		Rscript fit_cv_fused_lasso.R \
			--gaudi-path . \
			--gwas data/GWAS_chr#_wbc.regenie \
			--gwas-col-id ID --gwas-col-p LOG10P --gwas-log-p TRUE \
			--la data/test.la_dosage.mtx.gz \
			--col data/test.la_dosage.colnames --row data/test.la_dosage.rownames \
			--pheno data/test.pheno --pheno-name wbc_adj --pheno-fid FID \
			--start-p-exp -5 --end-p-exp -50 \
			--seed 2022 --sparsity FALSE \
			out data/test_p5_50_model

## Instructions on model testing (applying GAUDI model)

1. Combine HDS (or GT) from VCFs with local ancestry information

This is the same as step 2 in model training part. 
Note that it's not required to combine by-chunk la-dosage files.

2. Apply GAUDI model

Applying GAUDI model is very simple. See `apply_GAUDI.R`.
Just need to provide the trained models and target la-dosage information.
Example code to run:

		Rscript apply_GAUDI.R \
			--model data/test_p5_50_model.best_list.RDS \
			--target-la-dir TEST_DATA_LA_DIR \
			--out OUTPUT_FILE_NAME


