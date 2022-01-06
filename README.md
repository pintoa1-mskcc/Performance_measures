







Tool to calculate performance metrics on different slices of data. Has additional functionality to genotype the selected samples.





ARGUMENT PARSING:
Must submit a 'ground' and 'test' MAF. Optional flags include providing an output directory (defaults to current working directory), a output prefix (defaults to current date and time to prevent overwritting previous runs), and a bed file for processing on and off target results. 





Basic Run:
If you do not wish to run fillouts the following will occur:

Identify a restrictive and permissive tag for each variant. Restrictive tags are tags identified asa: Chrom:Start:End:Ref:Alt:TumorSampleBarcode, permissive tags are: Chrom:Start:End:Ref:TumorSampleBarcode. Statisticts will be calculated for both restrictive tags and permissive tags automatically, however only restrictive statistics are plotted.

Additional identified per variant are calculated: 
	OncogenicTF = Sets variants to ONCOGENIC if the oncogenic column is "Likely Oncogenic", "Oncogenic" or "Predicted Oncogeic". 
	Clonality = Defined as INDETERMINATE if cf is NA or cf is < 0.6 * purity. Defined as CLONAL if ccf_expected_copies > 0.8 | (ccf_expected_copies > 0.7 & ccf_expected_copies_upper > 0.9, SUBCLONAL otherwise
	is_non_syn_mut = Determines if a variant is a nonsynonmous mutation or not. Variant classification must be one of the following: "Missense_Mutation", 
                                                               "Nonsense_Mutation", 
                                                               "Nonstop_Mutation", 
                                                               "Frame_Shift_Ins", 
                                                               "Frame_Shift_Del",
                                                               "In_Frame_Del",
                                                               "In_Frame_Ins",
                                                               "Translation_Start_Site",
                                                               "Splice_Site"

	substitutions = Assesses the ref to alt change and assigned it to the designated substitution

Some other more complicated variables include defining purity and t_variant_frequency as bins. We determined for the sake of this tool that the "ground" file is considered an absoulte truth, so if a variant_frequency is defined in ground for a variant, that also exists in the test set, we will set that variant frequency in the test set to the ground set (for the sake of PR metrics). This is to detemine at each VAF_bin (5% bins), that whatever is found in ground is actually found in truth, regardless of truths calculated variant frequency. 

Similarly for purity, we assign sample purity according to the ground sample purity, so we can assess at each purity bin how well the test set can recapture those variants. 

PLOTTING AND CALCULATING STATS
Basic mode will automatically create plots based on these variables: 'oncogenic_tf','clonality','substitutions','is_non_syn_mut','t_var_freq_bin', and purity_bin.  Additionally it will provide an "overview" in which it summaries the PR_metrics over SNPs/indels/all variants. When plotting overview you will see SNPS/indels and all. All represents the union results of SNPs and indels together. For each variable listed above, the tool will also break down the statistics for SNPs/indels/all. When plotting results for each variable, the tool only visualizes the results for all, results for SNPs and indels will be saved in the .txt file.  PDFs generated for these cohort level statistics will be a bar plot

Once calculations are complete on the cohort level, we go back and calculate the performance measures for each individual sample provided. PDFs for these sample level statistics will be a boxplot. 

Finishing
The last step of this script will save the annotated MAF for ground and test to your output directory. This will state all the variables that were geneated and tell you whether or not the variant was called in the other provided maf. These mafs will be saved with the suffixes "ground_PR_labeled.maf" and "test_PR_labeled.maf"


ADDING A BED FILE
If a bed file is provided, the tool will identify off and on target reads. AS of Nov 11, 2021, the bedr tool in R has a slight bug in it when working with 4.0+. This bug doesnt effect its functionality beyond not being able to return a proper R object to the environment. I found that the tool works perfectly when indicating an outputFile. Therefore, when a bed file is provided, the tool will automatically generate a ground.bed file and a test.bed file. These files represent the locations that have been found to be within the bed file you provided, with a bed_tag so we can merge them back in. 
The result of the bed file is the addition of the column on_target, which states whether or not a variant is on target within the provided bed file. 

There will be a pdf/txt file generated with on_target/off_target performance measures. However, once this file is generated, the remaining plots will ONLY be PR measures on the ON TARGET variants (you will no longer be assessing all variants, only those that are said to be on target). To get a full picture of PR measures for off target results, run the tool again without the bed file.

FILLOUTS
To run fillouts, you must do the following:
	-f flag set
	--test_fillout_mapping, provide a txt file with the Tumor/normal IDs and bam locations for the test file. MUST CORRESPOND TO THE -t test.maf!
	--ground_fillout_mapping. provide a txt file with the Tumor/normal IDs and bam locations for the ground file. MUST CORRESPOND TO THE -g ground.maf!

The first step in this fillouts command is to generate a target MAF file for fillouts to genotype. To do this, I combine the unique variants to get alist of all variants found in both test and ground. We check to ensure that all the provided BAMs in your mapping file exist in our MAFs. If a tumor/normal from your mapping file does NOT exist in the provided MAFs, a warning will be printed and that sample will not have fillouts run since no targets are provided.  HOWEVER, if a sample exists in our provided MAFs but does not exist in the mapping files provided, the tool will throw an error and stop. 

While fillouts is running, the basic analysis will continue. Within your providedd output directory, a fillouts directory will be created. Once fillouts is complete for all samples, we proceed to the next step to restructure and gather the fillouts results ( see fillouts_restruct.R) 

Fillouts_restruct restructures the fillouts results so tumor normal results per variant are on the same line. This will also merge in your annotations from the original mafs. 

The final result of fillouts_restruct.R is two MAFs that have been genotypes corresponding to your test and ground MAFs. These two mafs will then go back to the original script to have PR metrics run. 

When returning to the original script, fillouts  results will be slightly different. We will first calculate EVIDENCE (or evidence that an alt allele exists, i.el t_alt_count >= 1), then we determing detectablity (t_total_count >= 20). 
We utilize detectablity and evidence to identify our tps, fps and fns for our performance measures. TP is defined as a variant that has evidence in both ground and test. FP is defined as a variant with evidence in TEST but no evidence in GROUND while also being detectable in GROUND. FN is defined as a variant with evidence in GROUND but is detectable in TEST but has NO EVIDENCE. Variants that have no evidence and are not detectable are removed from analysis. 

Statistics and graphs are otherwise treated the same way as the original method. Final MAFs returned will have additional columns Detectable_in_Other_Run and Evidence_in_Other_Run. 



BED AND FILLOUTS
If you provided a BED file and want fillouts run, fillouts will be run on ALL locations, then filtered during performance measures. 

