#!/usr/bin/env python3

#################
# Version 0.3 ###
#################

# Added option to include a set of sample ID's as an inclusion criteria.
# Added option to only consider chunks of a given chromosome.

#################

#################
# Version 0.2 ###
#################

# Changed output file format from rows:subjects columns:variants to rows:variants and columns:subjects #
# Implemented writing output file in chunks. #
# Implemented local-ancestry minor allele count filtering. # 

#################
import pysam
import pandas as pd
import time
import numpy as np
import itertools as it
import sys
import argparse

def parse():
	parser = argparse.ArgumentParser(
		description = "Combine haplotype dosage values with local ancestry estimates from RFMix.")

	parser.add_argument("--local-ancestry", action = "store", type=str, required = True,
				help = "The RFMix generated .msp.tsv.gz file subset to only analysis subjects.")
	parser.add_argument("--vcf", action = "store", type=str, required = True,
				help = "VCF file containing Minimac4 style HDS dosages.")
	parser.add_argument("--keep", action = "store", type = str, required = False, 
				help = "Tab delimited file of sample ID's matching those in the .vcf and RFMix files.")
	parser.add_argument("--remove", action = "store", type = str, required = False, 
				help = "Tab delimited file of sample ID's to remove.")
	parser.add_argument("--include", action = "store", type = str, required = False,
				help = "Tab delimited file of variant ID's to include.")
	parser.add_argument("--chr", action = "store", type = str, required = False, default = None)
	parser.add_argument("--pos-start", action = "store", type = int, required = False, default = None)
	parser.add_argument("--pos-stop", action = "store", type = int, required = False, default = None)
	parser.add_argument("--out", action = "store", type=str, required = True, 
				help = "Out file name for LA Dosage file (no suffix needed)")
	parser.add_argument("--la-dosage-threshold", action = "store", type = int, default = 10, 
				help = "Threshold for filtering variants. At least one ancestry group must have k individuals with the minor allele.")

	return parser.parse_args()

args = parse()


#Read in RFMix local ancestry data.
rf = pd.read_csv(args.local_ancestry, sep = "\t")

#Create pysam VCF object.
vcf = pysam.VariantFile(args.vcf)

print("Creating sample vectors")
samples = list(vcf.header.samples)

if args.keep is not None:
	keep_ids = [line.rstrip() for line in open(args.keep)]
	samples = [s for s in samples if s in keep_ids ]

if args.remove is not None:
	remove_ids = [line.rstrip() for line in open(args.remove)]
	samples = [s for s in samples if s not in remove_ids ]

if args.include is not None:
	include_snps = [line.rstrip() for line in open(args.include)]


vcf_sample_index = [ i for i,s in enumerate(list(vcf.header.samples)) if s in samples]
rf_sample_index =[ i for i,col in enumerate(rf.columns.values) if col.split(".")[0] in samples ]


print("Subsetting rf file on samples")
#Subset rf file on samples.
rf = rf.iloc[:, [0,1,2] + rf_sample_index ]


la_var_df_list = []

start = time.time()
i=0
chunk_out=1
print("Starting to read vcf file")
for variant in vcf.fetch(contig = args.chr, start = args.pos_start, stop = args.pos_stop):
	#If include argument is included, check if variant in include_snps
	if args.include is not None:
		if not variant.id in include_snps:
			print("Variant %s not in --include snps, moving to next variant" % variant.id)
			continue

	#Check to see if the variant is covered by the RFMix intervals
	is_in=( rf["spos"] <= variant.pos ) & ( variant.pos < rf["epos"] )

	#If not, move on to the next variant	
	if not any(is_in): 
		print("Variant %s not in RFMix output, moving to next variant" % variant.id)
		continue

	#Otherwise, create a pandas dataframe
	row_names = [ variant.id + "_AFR", variant.id + "_EUR" ] 

	#Create output dataframe with 2 rows and columns corresponding to the number of samples.
	la_var = pd.DataFrame(np.zeros((2, len(samples))), 
				columns = samples, 
				index = row_names)
	#Check if variant has HDS dosage values
	if not ("HDS" in variant.format.keys()):
		#If we only have Genotypes, use those. Assume missing is 0.
		print("Variant %s does not have HDS dosage values, only considering GT values" % variant.id)
		#For each sample, check if HDS value is present. If so, return the hap1 HDS value. Otherwise, return the GT value.
		hds_hap1 = [ x["GT"][0] if (x["GT"][0] is not None) else 0 for x in variant.samples.values() ]
		#Then subset to just samples we care about.
		hds_hap1 = [ hds_hap1[i] for i in vcf_sample_index ]

		#Same process for hap2.
		hds_hap2= [ x["GT"][1] if (x["GT"][1] is not None) else 0 for x in variant.samples.values() ]
		hds_hap2 = [ hds_hap2[i] for i in vcf_sample_index ]
	else:
		#Use HDS if available. Assume missing is 0. Use GT if HDS is missing. Use 0 if both are missing. 
		#For each sample, check if HDS value is present. If so, return the hap1 HDS value. Otherwise, return the GT value.
		hds_hap1 = [ x["HDS"][0] if (len(x['HDS'])==2) else x["GT"][0] if (x["GT"][0] is not None) else 0 for x in variant.samples.values() ]
		#Then subset to just samples we care about.
		hds_hap1 = [ hds_hap1[i] for i in vcf_sample_index ]
	
		#Same process for hap2.
		hds_hap2= [ x["HDS"][1] if (len(x['HDS'])==2) else x["GT"][1] if (x["GT"][1] is not None) else 0 for x in variant.samples.values() ]
		hds_hap2 = [ hds_hap2[i] for i in vcf_sample_index ]


	#Get corresponding row for RFMix data.
	rf_interval = rf[( rf["spos"] <= variant.pos ) & ( variant.pos < rf["epos"] )]

	#Get column names corresponding to sample values
	rf_cols = list(rf_interval.columns)[3:]

	#Get values for Hap1 estimates and Hap2 estimates.
	rf_cols_hap1 = [col for col in rf_cols if ".0" in col ]
	rf_cols_hap2 = [col for col in rf_cols if ".1" in col ]

	rf_vals_hap1 = rf_interval[rf_cols_hap1].iloc[0,:]
	rf_vals_hap2 = rf_interval[rf_cols_hap2].iloc[0,:]

	#Add the HAP1 HDS dosage values from AFR haplotypes to the AFR output column.
	la_var.loc[row_names[0], list(rf_vals_hap1 == 0)] = np.add(la_var.loc[row_names[0],list(rf_vals_hap1 == 0)],
								list(it.compress(hds_hap1, list(rf_vals_hap1 == 0))))

	#Add the HAP2 HDS dosage values from AFR haplotypes to the AFR output column.
	la_var.loc[row_names[0], list(rf_vals_hap2 == 0)] = np.add(la_var.loc[row_names[0], list(rf_vals_hap2 == 0)],list(it.compress(hds_hap2, list(rf_vals_hap2 == 0))))	

	#Add the HAP1 HDS dosage values from EUR haplotypes to the EUR output column.
	la_var.loc[row_names[1], list(rf_vals_hap1 == 1)] = np.add(la_var.loc[row_names[1], list(rf_vals_hap1 == 1)],list(it.compress(hds_hap1, list(rf_vals_hap1 == 1))))

	#Add the HAP2 HDS dosage values from EUR haplotypes to the EUR output column.
	la_var.loc[row_names[1], list(rf_vals_hap2 == 1)] = np.add(la_var.loc[row_names[1], list(rf_vals_hap2 == 1)],list(it.compress(hds_hap2, list(rf_vals_hap2 == 1))))

	#If a variant isn't observed in enough local-ancestry haplotypes, drop the variant.

	#Check if variant is observed in enough local-ancestry haplotypes in at least one ancestry group.
	if any(la_var.astype(bool).sum(axis = 1) >= args.la_dosage_threshold):
		#But remove any rows with no observations. 
		keep_rows=(la_var.astype(bool).sum(axis = 1) != 0)
		la_var = la_var.loc[keep_rows,:]
		la_var_df_list.append(la_var)
	else:
		print("Variant %s does not pass la_dosage threshold" % variant.id )
		continue

	# If i is divisible by 100
	if i % 500 == 0 and i > 0:	
		
		#Combine all the dataframes so far
		out_df = pd.concat(la_var_df_list, axis = 0)
		#if it's the first time I'm doing this, write to a file with header
		if chunk_out == 1:
			print("Writing chunk " + str(chunk_out))
			out_df.to_csv(args.out + ".la_dosage.tsv.gz", sep = "\t", 
				index_label = "SNP", header = True)
			chunk_out += 1
		#Otherwise, append to the file and don't include headers.
		else: 
			print("Writing chunk " + str(chunk_out))
			out_df.to_csv(args.out + ".la_dosage.tsv.gz", sep = "\t",
				index_label = "SNP", header = False, mode = "a")
			chunk_out+=1
				
		#Then reset the la_var_df_list for the next chunk. 
		la_var_df_list=[]
	i += 1

if la_var_df_list == []:
	print("No chunks to write. Ending.")
	end = time.time()
	print(end - start)
	sys.exit()

#Then, write the last bit.
out_df = pd.concat(la_var_df_list, axis = 0)
if chunk_out == 1:
	print("Writing chunk " + str(chunk_out))
	out_df.to_csv(args.out + ".la_dosage.tsv.gz", sep = "\t",
		index_label = "SNP", header = True)
else:
	print("Writing chunk " + str(chunk_out))
	out_df.to_csv(args.out + ".la_dosage.tsv.gz", sep = "\t",
		index_label = "SNP", header = False, mode = "a")

end = time.time()
print(end - start)



