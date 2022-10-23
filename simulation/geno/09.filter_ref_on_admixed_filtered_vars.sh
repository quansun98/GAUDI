#!/bin/bash

module add plink

pop=$1 # population, EUR or AFR
locus_num=$2 # same as in step 1
sim_name=$3 # same as in step 1
r2=$4 # r2 threshold for LD pruning

locus_dir=../../../data/simulated_genotypes/${sim_name}/locus${locus_num}

#Filter vars from admixed samples
cut -f 1 ${locus_dir}/locus${locus_num}_admixed_r2${r2}.var | tail -n +2 > ${locus_dir}/filter_vars

plink --bfile ${locus_dir}/${pop}_ref_LDPruned${r2} --extract ${locus_dir}/filter_vars --make-bed --out ${locus_dir}/${pop}_ref_LDPruned${r2}_admixedFiltered

plink --bfile ${locus_dir}/${pop}_ref_LDPruned${r2} --extract ${locus_dir}/filter_vars --recode A --out ${locus_dir}/${pop}_ref_LDPruned${r2}_admixedFiltered

#rm ${locus_dir}/*.nosex ${locus_dir}/*.log filter_vars



