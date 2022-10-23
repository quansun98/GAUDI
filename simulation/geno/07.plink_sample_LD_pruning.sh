#!/bin/bash

module add plink

locus_num=$1 # same as in step 1
sim_name=$2 # same as in step 1
max_locus_num=$3 # maximum locus number

locus_dir=../../data/simulated_genotypes/${sim_name}/locus${locus_num}

awk -v locus_num=$locus_num '{ print locus_num, "locus_"locus_num"_pos_"$3, 0, $3 }' <( tail -n +2 ${locus_dir}/${sim_name}_locus${locus_num}.pos-1 ) > ${locus_dir}/sim_haps_all_subjects.map

#Get file for LD pruning. For this simulation, we want to be sure that we're excluding all variants in LD. So set the LD window for --indep-pairwise to the entire region (500kb).
plink --file ${locus_dir}/sim_haps_all_subjects --chr-set ${max_locus_num} --indep-pairwise 500 5 0.1 --out ${locus_dir}/sim_haps_all_subjects

rm ${locus_dir}/*.nosex ${locus_dir}/*.log



