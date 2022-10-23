### Genotype Simulation

This directory contains codes for simulating genotypes in the GAUDI manuscript. Detailed notations for each script are as following.

1. `01.runCosi_sim_setup.sh` is used to simulate genotypes (European and African in this example). This bash script is a wrapper of the `runCosi.pl` script. Details about COSI are described in the paper with PMID 16251467.

2. `02.generate_COSI_hap_IDs.R`: Since COSI sets chromosome IDs as numbers between 0 and (2n - 1), where n is the number of subjects, this script decides which chromosomes will be fully African Ancestry, fully European Ancestry, and admixed chromosomes. Chromosomes which are saved to serve as reference chromosomes for partial PRS are also selected. It does not determine which subject receives each chromosome (this will be done in step6).

3. `03.get_AFR_EUR_COSI_haps.sh`: Subsets the COSI output files on the IDs generated in step 2.

4. `04.ref_files_to_plink_format.R`: Converts COSI files to .map and .ped files which can be read by PLINK. This makes working with a common set of tools easier rather than having to reinvent the wheel at every step. 

5. `05.create_switch_point_genotypes.R`: This creates admixed chromosomes using a random procedure. Given that a chromosome is admixed, the "left" ancestry is determined randomly between AFR and EUR, then a switch point is chosen along the chromosome segment. Then, the two contributing chromosomes are combined with respect to the haploytype ordering. 

6. `06.randomize_haplotypes.R`: For each subject, randomly select two chromsomes. 

7. `07.plink_sample_LD_pruning.sh`: Run LD pruning for a given r2 cutoff. Ties are broken based off in-sample MAF.  

8. `08.summarize_haplotypes.R`: Creates a host of useful downstream files based off of the genotypes that were generaated in step 6. A big utility script for this procedure. 

9. `09.filter_ref_on_admixed_filtered_vars.sh`: We only want to consider the reference genotypes at SNPs which survived filtering in our training subjects. 

A note: 
`pop_num_key`: Confusingly, COSI assigns populations as numbers, this is a key to remember which is which between AFR and EUR.
 
