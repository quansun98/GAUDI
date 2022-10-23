#!/bin/bash

locus_num=$1 # same as in step 1, locus number (for multiple loci)
sim_name=$2 # same as in step 1, unique name to distinguish different simulation settings

locus_dir=../../data/simulated_genotypes/${sim_name}/locus${locus_num}

# Note that the output prefix from step 2 is ${locus_dir}/locus${locus_num}_COSI_hap_IDs. Please changed as needed.

echo Locus number: $locus_num

# Subset .hap file - EUR Haplotypes
awk 'NR==FNR {a[$1]; next} $1 in a' ${locus_dir}/locus${locus_num}_COSI_hap_IDs_eur ${locus_dir}/${sim_name}_locus${locus_num}.hap-1 > ${locus_dir}/subset_eur

#Reformat file and add informative column names for the SNPs. Then remove trailing whitespace. 
cat <( paste -d ' ' <( echo chr_id left_ancestry switch_point_pos ) <( tail -n +2 ${locus_dir}/${sim_name}_locus${locus_num}.pos-1 | cut -f 3 | sed "s/^/locus_${locus_num}_pos_/g" | tr '\n' ' ' ) ) \
  <( paste -d ' ' <( cut -f -2 ${locus_dir}/subset_eur| awk '{ print "EUR_"NR, "EUR", "NA" }' ) <( cut -f 3- ${locus_dir}/subset_eur ) )  | sed 's/ *$//g' > ${locus_dir}/sim_haplotypes_eur


#Subset .hap file - EUR  ref
awk 'NR==FNR {a[$1]; next} $1 in a' ${locus_dir}/locus${locus_num}_COSI_hap_IDs_ref_eur ${locus_dir}/${sim_name}_locus${locus_num}.hap-1 > ${locus_dir}/subset_ref_eur

#
cat <( paste -d ' ' <( echo chr_id left_ancestry switch_point_pos ) <( tail -n +2 ${locus_dir}/${sim_name}_locus${locus_num}.pos-1 | cut -f 3 | sed "s/^/locus_${locus_num}_pos_/g" | tr '\n' ' ' ) ) \
  <( paste -d ' ' <( cut -f -2 ${locus_dir}/subset_ref_eur | awk '{ print "EUR_"NR, "EUR", "NA" }' ) <( cut -f 3- ${locus_dir}/subset_ref_eur ) )  | sed 's/ *$//g' > ${locus_dir}/sim_haplotypes_ref_eur

#Same process for AFR
awk 'NR==FNR {a[$1]; next} $1 in a' ${locus_dir}/locus${locus_num}_COSI_hap_IDs_afr ${locus_dir}/${sim_name}_locus${locus_num}.hap-5 > ${locus_dir}/subset_afr

cat <( paste -d ' ' <( echo chr_id left_ancestry switch_point_pos ) <( tail -n +2 ${locus_dir}/${sim_name}_locus${locus_num}.pos-5 | cut -f 3 | sed "s/^/locus_${locus_num}_pos_/g" | tr '\n' ' ' ) ) \
 <( paste -d ' ' <( cut -f -2 ${locus_dir}/subset_afr | awk '{ print "AFR_"NR, "AFR", "NA" }' ) <( cut -f 3- ${locus_dir}/subset_afr ) ) | sed 's/ *$//g' > ${locus_dir}/sim_haplotypes_afr

#Subset .hap file - AFR  ref
awk 'NR==FNR {a[$1]; next} $1 in a' ${locus_dir}/locus${locus_num}_COSI_hap_IDs_ref_afr ${locus_dir}/${sim_name}_locus${locus_num}.hap-5 > ${locus_dir}/subset_ref_afr

#
cat <( paste -d ' ' <( echo chr_id left_ancestry switch_point_pos ) <( tail -n +2 ${locus_dir}/${sim_name}_locus${locus_num}.pos-1 | cut -f 3 | sed "s/^/locus_${locus_num}_pos_/g" | tr '\n' ' ' ) ) \
  <( paste -d ' ' <( cut -f -2 ${locus_dir}/subset_ref_afr | awk '{ print "AFR_"NR, "AFR", "NA" }' ) <( cut -f 3- ${locus_dir}/subset_ref_afr ) )  | sed 's/ *$//g' > ${locus_dir}/sim_haplotypes_ref_afr

#Clean Up
rm ${locus_dir}/subset_eur ${locus_dir}/subset_afr ${locus_dir}/subset_ref_afr ${locus_dir}/subset_ref_eur
