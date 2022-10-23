#!/bin/bash

module add plink

#This script generates genotype PC's for training and testing data seperately for phenotype adjustment. 

prs_sim_name=ideal_FL
cosi_sim_name=AFR_2600_EUR_1400
sim_key=../../data/COSI_${prs_sim_name}_locus_key

while read n_loci; do

	out_name=../../data/simulated_genotypes/causal_genotypes/${prs_sim_name}/COSI_causal_genotypes_nloci_${n_loci}
	if [ $n_loci -gt 1 ]; then
		#clear previous merge list
		rm -f temp_merge_list 
		while read n_loci locus_num; do
			locus_dir=/proj/yunligrp/users/bryce/prs/admixed_prs/data/simulated_genotypes/${cosi_sim_name}/locus${locus_num}
			echo ${locus_dir}/locus${locus_num}.ped ${locus_dir}/locus${locus_num}.map >> temp_merge_list
	
		done < <( awk -v n_loci=$n_loci '$1==n_loci { print $1, $3}' $sim_key  | sort -u )	
		#After I've gotten all the loci, merge them. 
		plink --merge-list temp_merge_list --make-bed --out ${out_name} 

	else
		while read n_loci locus_num; do 
			locus_dir=/proj/yunligrp/users/bryce/prs/admixed_prs/data/simulated_genotypes/${cosi_sim_name}/locus${locus_num}
			plink --file ${locus_dir}/locus${locus_num} --make-bed --out ${out_name}
		done < <( awk -v n_loci=$n_loci '$1==n_loci { print $1, $3}' $sim_key  | sort -u )
	fi

done < <( tail -n +2 $sim_key | cut -f 1 | sort -u ) 


#For loop to do LD pruning and PC generation for the 5 folds in the training data. 
while read n_loci rep; do
	data_dir=../../data/simulated_genotypes/causal_genotypes/${prs_sim_name}
	pheno_dir=../../data/simulated_phenotypes/${prs_sim_name}/COSI_phenos_nloci_${n_loci}
	GWAS_dir=../../data/cosi_sim_results/GWAS_marginal_screening/${prs_sim_name}/COSI_sims_nloci_${n_loci}

	mkdir -p $GWAS_dir
	mkdir -p ${data_dir}/cv_training_data

	echo $n_loci $rep 
	#Subset to just the training data. 
	plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci} \
		--keep ${pheno_dir}/COSI_sims_nloci_${n_loci}_training_ids \
		--make-bed \
		--out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training

	#Run for full training data, so that we can get lambda sequence for fused-lasso.
	#LD Pruning for genotype PCs
	plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training \
		--indep-pairwise 50 5 0.1 \
		--out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training

	#Estimate genotype PCs
	plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training \
                --extract ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training.prune.in \
		--pca 10 \
		--out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_training_LD_pruned


	#Run same procedure for testing data
	plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci} \
                --keep ${pheno_dir}/COSI_sims_nloci_${n_loci}_testing_ids \
                --make-bed \
                --out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing

	
	#LD Pruning for genotype PCs
        plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing \
                --indep-pairwise 50 5 0.1 \
                --out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing

	
	#Estimate genotype PCs
        plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing \
                --extract ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing.prune.in \
                --pca 10 \
                --out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_testing_LD_pruned

	#And also for the full dataset.
        #LD Pruning for genotype PCs
        plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci} \
                --indep-pairwise 50 5 0.1 \
                --out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}


        #Estimate genotype PCs
        plink --bfile ${data_dir}/COSI_causal_genotypes_nloci_${n_loci} \
                --extract ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}.prune.in \
                --pca 10 \
                --out ${data_dir}/COSI_causal_genotypes_nloci_${n_loci}_LD_pruned


done < <( tail -n +2 $sim_key | cut -f 1,2 | sort -u ) 


