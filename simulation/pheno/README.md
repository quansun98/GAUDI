### Phenotype Simulation

This directory contains codes used for simulating phenotypes in the GAUDI manuscript. Detailed explanations are as following.

1. `01.pheno_simulation.R`: This is the script used to generate phenotypes given simulation parameters. It may need to run other files in this subdirectory before all the prerequisite files are generated. 
Basically, each row is a full set of simulated data. 
Rows are different simulation hyperparamters: hertiability, proportion of SNPs with ancestry shared effects, proportion of causal  SNPs, etc. 

2. `02.plink_preprocessing.sh`: Generates important covariates based on the phenotypes generated in 01, like genotype PCs. Also does LD pruning.

3. `03.adjusting_phenotypes.R`: Adjusts phenotypes based on genotype PCs like one will do in a normal study and then uses inverse-normal transformation on the residuals of the lm model.

4. `04.apply_betas_to_ref_files.R`: Simple way of creating phenotypes for reference individuals (non-admixed single-ancestry). Simulates noise too. 

