#!/bin/bash

# This script provides a wrapper of the runCosi.pl script if you are not familiar with pearl script.

cosi_path=$1 # path to cosi
locus_number=$2 # locus number (for multiple loci)
sim_setup=$3 # unique name to distinguish different simulation settings
in_recParam=$4 # input file containing recombination parameters
in_param=$5 # input file containing parameters

out_dir=../../data/simulated_genotypes/${sim_setup}/locus${locus_number}

#Run Cosi with the following positional parameters:
#  0. Cosi path
#  1. Locus Number
#  2. Cosi Ancestry Index - pop1 
#  3. Cosi Ancestry Index - pop2
#  4. recParams file
#  5. params file
#  6. output prefix

# 1 Corresponds to European ancestry and 5 corresponds to African ancestry. 
./runCosi.pl $cosi_path $locus_number 1 5 $in_recParam $in_param $out_dir

