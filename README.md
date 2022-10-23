# GAUDI

## Introduction

GAUDI is a novel polygenic risk score (PRS) method designed specifically for admixed individuals.
It takes the mosaic structure of the admixed genomes into consideration, leveraging local ancestry information to jointly estimate ancestry-shared effects and allowing ancestry-specific effects.
GAUDI model is based on a modified fussed lasso framework, balancing between fusion and sparsity.
The sparsity penalty performs variant selection, and the fusion penalty encourages similarity of the effect sizes of the same variant between two populations.
See our manuscript for more details: [GAUDI preprint](<https://www.biorxiv.org/content/10.1101/2022.10.06.511219v1.abstract>)

This folder also contains codes used in the simulation studies. See `simulation` for more details.

## Requirement

The below R packages and version are required:

* R >= v4.1.0
* tidyverse >= v1.3.1
* Matrix >= 1.4.0
* optparse >= v1.7.1
* data.table >= v1.14.2
* genlasso >= v1.6.1
* splitTools >= v0.3.2
* caret >= v6.0.90

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
Note that VCFs need to be block-zipped (bgzip) and have index (.tbi file).
Local ancestry inferred by RFMix is required. The below codes provide an example to run (by chunk to require less memory):

		while read chunk; do
			start_pos=$(( ($chunk - 1)*10*1000000 ))
			end_pos=$(( ($chunk)*10*1000000 ))
			./py_vcf_to_la.py \
				--local-ancestry data/local_anc_chr22.msp.tsv.gz \
				--vcf data/test_chr22.vcf.gz \
				--include data/chr22_minMAF0.001_Rsq0.9_GWASp0.05.txt \
				--la-dosage-threshold 5 \
				--chr 22 --pos-start $start_pos --pos-stop $end_pos \
				--out data/test_out_chr22_chunk${chunk}
		done < <( awk -v chr=$chr '$1==chr { print $2 }' data/chunk_list

Note that the `--la-dosage-threshold` flag is created to filter variants. At least one ancestry group must have k individuals with the minor allele.
For example, for a specific variant, if AFR-MAC = 50, EUR-MAC = 1, this variant will be kept using `--la-dosage-threshold 10`. 
The above example uses only 1 because of the small sample size in the toy dataset. We used 10 in our real data analysis.

3. Combine files created in the last step

The script `merge_la_dosage.R` combines the by chunk `.la_dosage.tsv.gz` files together. 
All files under the input directory with suffix `.la_dosage.tsv.gz` will be considered to be combined.
It will combine all chunk files together (genome-wide).
The chromosome without any chunk will be ignored.

		Rscript merge_la_dosage.R data

A sparse matrix output with column names and row names will be generated.

4. Run GAUDI for model training

The main function to run GAUDI is `fit_cv_fused_lasso.R`. 
It provides many customized options. Please see the script for details.
It also calls `cv_fused_lasso.R` internally.
An R object `OUTPUT.best_list.RDS` will be generated, which stores a list containing both fitted models and the best p-value threshold in the grid search process.
An example code to train GAUDI model is provided below:


		Rscript fit_cv_fused_lasso.R \
			--gaudi-path . \
			--gwas data/GWAS_chr22_sim.regenie \
			--gwas-col-id ID --gwas-col-p LOG10P --gwas-log-p TRUE \
			--la data/test_out_allChr.la_dosage.mtx.gz \
			--col data/test_out_allChr.la_dosage.colnames \
			--row data/test_out_allChr.la_dosage.rownames \
			--pheno data/test.pheno --pheno-name pheno --pheno-iid FID \
			--start-p-exp -1 --end-p-exp -5 \
			--seed 2022 --sparsity FALSE \
			--out data/test_p5_50_model

Note that: 

(1) You can specify chromosome-separated GWAS results together using `#` to replace the chromosome numbers.
For example, `--gwas data/GWAS_chr#_sim.regenie`.

(2) The `--start-p-exp` and `--end-p-exp` specify starting and ending p-value thresholds for grid search to maximize the PRS performance.

(3) GAUDI will automatically remove highly correlated variants (LD > 0.95) in training samples.
LD pruning or clumping could be performed before training models, but is not required.

## Instructions on model testing (applying GAUDI model)

1. Combine HDS (or GT) from VCFs with local ancestry information

This is the same as step 2 in model training part. 
Note that it's not required to combine by-chunk la-dosage files.

2. Apply GAUDI model

Applying GAUDI model is very simple. See `apply_GAUDI.R`.
Just need to provide the trained models and target la-dosage information.
Example code to run:

		Rscript apply_GAUDI.R \
			--gaudi-path . \
			--model data/test_p5_50_model.best_list.RDS \
			--target-la-dir data/ \
			--out data/test_self_fit

Note that this is just for code illustration purpose.
Please do not fit PRS on the same samples that were used in model training.

## More information on GAUDI model output

As mentioned in step 4 in the model training part, an R object (in `list` type) will be generated.
The main output is the `fit_model`, which contains the below objects.

1. `fit`: this is an output format from the `gen_lasso` package.
We here present some basic interpretations of the outputs.
Please refer to the `genlasso` package for more details.
The first three objects in `fit` are

(1). `lambda`: the lambda vectors under the best-fit p-value threshold
Note that all the information stored is everything under this best-fit p-value threshold.

(2). `beta`: an `m \times k` matrix, where m is the number of predictors and `k` is the number of different \lambda values.
Note that the predictors are not equal to variants.
Each predictor has names `variantID_ancestry` to allow ancestry-specific effects.
One variant may only have one ancestry effect presented in the final model.
That is one unique feature of GAUDI and what distinguishes GAUDI with others.

(3). `fit`: an `n \times k` matrix, where `n` is the number of individuals.
It contains individual level fitted PRS values under each \lambda value.

2. `D`: the penalty matrix in sparse matrix format

3. `snps`: predictors in the final model (not equivalent to variants)

4. `best_gamma`: the \gamma values that achieves the best PRS performance (in terms of correlation with phenotypes)

5. `best_lambda`: the \lambda values that achieves the best PRS performance (in terms of correlation with phenotypes)

6. `cv_r2`: the cross-validation R2 under the best parameter

7. `full_cv_r2`: all the cross-validation R2 under the parameter grids.
Each row is a unique \lambda value, and each column is the \gamma value.


