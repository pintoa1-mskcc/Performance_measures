# Performance Measure Metrics
This tool is meant to calculate performance metrics on different slices of the input data. The tool expected 2 MAF files as the input to compare recall and precision between the two. The results are then subsectioned based on annotations for a clearer picture. If the FILLOUTS flag is specified, the tool will run metrics on both the CALLED (or original MAF file) and the resulting genotyped MAF files. 

NOTE: This analysis must be run on JUNO if utilizing genotyping or HTML report generation/multiqc.

## Pipeline Flowchart
<p align="center">
  <img src="./docs/performance_measure_workflow.png"/>
</p>


## Arguments for performance_measure_script.R
### Required
R 4.0+
If you want to run genotyping and/or multiqc, job must be submitted to JUNO scheduler. 
MAF file headers must match TEMPO output, at minimum *Chromosome;Start_Position;End_Position;Reference_Allele;Tumor_Seq_Allele, Tumor_Sample_Barcode; Variant_Type*


#### Basic Analysis

Analysis will only be run iff the following columns are present in the MAF: Chromosome;Start_Position;End_Position;Reference_Allele;Tumor_Seq_Allele, Tumor_Sample_Barcode; Variant_Type
  
```
-g/ --ground			Path to cohort ground .MAF file
-t/ --test			Path to cohort test .MAF file
```
#### Genotyping
```
-f/ --fillouts 			[Providing flag activates genotyping]
-r/ --ground_fillout_mapping	Path to ground tumor/normal BAM mapping .TXT file. See helper script below.
-e/ --test_fillout_mapping	Path to test tumor/normal BAM mapping .TXT file. See helper script below.
```

### Optional
```
-d/ --directory			Save results to this directory. Default is working directory
-n/ --name_ground		Name of ground cohort. Default "ground"
-s/ --name_test			Name of test cohort. Default "test"
-o/ --out_prefix		Out prefix. Default is Sys.time: %Y_%M_%D_%h:%m:%s_
-b/ --bed_file			Path to target .BED file. Activated on target analysis. See wiki below
-m/ --multiqc			  [Providing flag activates generation of a HTML report]
-v/ --additional_variables      A single or comma-seperated list of additional variables to calculate metrics on. These variables MUST be already in your .MAF files and be catagorical.
-P/ --Parallel      [Providing flag activates parallel analysis (Only recommended for large cohorts)
```

## Examples
### Basic Run:

`Rscript performance_measure_script.R --test test.maf --ground  ground.maf -n DRAGEN -s TEMPO -d base_run_dir_0`

#### Generate HTML report
`Rscript performance_measure_script.R --test test.maf --ground  ground.maf -n DRAGEN -s TEMPO -d base_mqc_run_dir_0 -m`

### Genotyping Run:

`Rscript performance_measure_script.R -t test.maf -g ground.maf -n DRAGEN -s TEMPO -f --ground_fillout_mapping ground_bam_mapping_for_fillouts.txt --test_fillout_mapping test_bam_mapping_for_fillouts.txt -d fillouts_run_dir_0  `

### On Target + Genotyping Run:

`Rscript performance_measure_script.R -t test.maf -g ground.maf -n DRAGEN -s TEMPO -f -r ground_bam_mapping_for_fillouts.txt -e test_bam_mapping_for_fillouts.txt -b IMPACT468_b37.bed -d fillouts_mqc_bed_run_dir_0`

### LSF Submission:
`bsub -o . -n 4 -R "rusage[mem=4]" -W 200:00 'Rscript performance_measure_script.R -t test.maf -g ground.maf -n DRAGEN -s TEMPO -f -r ground_bam_mapping_for_fillouts.txt -e test_bam_mapping_for_fillouts.txt -b IMPACT468_b37.bed -d fillouts_mqc_bed_run_dir_0'`

### Helper Script
[making_fillouts_mapping_file.R](./making_fillouts_mapping_file.R) 
This script can help generate a TXT mapping file. Provide the directory containing the test OR ground BAMs, a TXT file with the Tumor - Normal Pair (expected header TUMOR_ID \t NORMAL_ID), and the name of the output file. 

`Rscript  making_fillouts_mapping_file.R --bam_dir /juno/work/tempo/wes_repo/Results/v1.4.x/cohort_level/CCS_NFCXOZVP/bams/ --mapping example_inputs/CCS_NFCXOZVP.cohort.txt --output_file example_inputs/CCS_bam_mapping_for_fillouts.txt`

Sym links are not functional when genotyping. 

### System Requirements and Runtime
In a standard run, you will only need about 4G of RAM. With four threads, 200 samples runtimes around 10 minutes. When you run 200 sampels in parallel on 4 threads, you will need at least 24G of RAM. This will reduce your run time to 5 minutes. For this reason, we recommend only running parallel with cohorts > 1k with AT LEAST 36G RAM. (Note this parallel mode has not been tested with more than 1k). 

## Required R packages
- dplyr
- data.table
- stringr
- jsonlite
- binom
- ggplot2
- ggpubr
- ggsci
- tidyverse
- tidyr
- argparse
- doParallel
- cowplot
- grid
- gridExtra

If running on target analysis:
- bedR

# General Wiki
When we run this code, we made some assumptions about use cases. 
  ## 1. Ground Annotations are assumed to be the desired annotations
  The most important thing to remember when reviewing the results is when you are specifying your ground and test MAFs, the GROUND annotations are assumed to be the annotations you want to calculate recall and precision across. Therefore we set all _shared_ TEST annotations to the ground annotations. If a variant in tests exists but isn't in GROUND, the test annotation is retained. The significance of the annotations determine whether or not a variant is "recalled" when we are splitting annotations. A warning will be outputted in your log file to tell you when annotations have changed. Results can be analysed in the resulting MAFs.
  > E.G.: Say you have validated oncogenic variants in your GROUND file, but your test file was annotated prior with non-validated annotations. The script target if the variant is present or absent in the previous call, not if the annotation is correct. This allows us the assess the differences in performance only on presence or absence of a call. 
  ## 2.  Samples unique to test or ground are allowed
  Analysis will be preformed on all samples in ground and test, whether or not the sample is present in both cohorts. This will cause either recall or precision to be 0 for that sample whenever analysis is performed. A warning will be issued and the missing sample will be listed in logs/***_missing_samples.txt. NOTE IF YOU ARE RUNNING GENOTYING THIS RULE CHANGES. See [here](#place1)
  ## 3. Variables automatically run
  Data will also automatically be run for the following
  1. oncogenic_tf 
     - Performed if `oncogenic` found in MAF. 
     
     `ifelse(grepl("ncogenic", oncogenic), 'ONCOGENIC' ,'OTHER’)`
  2. Is_non_syn_mut 
     - Performed if `Variant_Classification` found in MAF
     
     ` ifelse(Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Splice_Site"), T, F))`
  3. Clonality 
     - Performed if `cf/ccf_expected_copies/ccf_expected_upper/purity` found in MAF
     
     ` mutate(clonality = ifelse(is.na(cf) | cf < (0.6 * purity), "INDETERMINATE", ifelse((ccf_expected_copies > 0.8 | (ccf_expected_copies > 0.7 & ccf_expected_copies_upper > 0.9)), "CLONAL", "SUBCLONAL"))) `
  4. Substitutions
     - Performed AUTOMATICALLY 
     ` ref_to_alt, `A>C` = "T>G", `T>G` = "T>G", `A>G` = "T>C", `T>C` = "T>C", `A>T` = "T>A", `T>A` = "T>A", `C>A` = "C>A", `G>T` = "C>A", `C>G` = "C>G", `G>C` = "C>G",`C>T` = "C>T", `G>A` = "C>T") `
   5. t_var_freq 
      - Performed automatically. Bins are formed in 5% intervals. 
   6. purity
      - Performed automatically if `purity` is present. Bins are formed in 5% intervals.  If purity is not present for a sample, sample will have a value of N/A
   ## 4. General Math
   True Positive (tps): Variant exists in Ground and Test MAF
   
   False Positive (fps): Variant exists in Test but not in Ground MAF
   
   False Negative (fns): Variant exists in Ground but not in Test MAF
   
   Recall 
   - `tps / (tps + fns)`
   - `binom.confint(tps,(tps + fns), conf.level = 0.95, method = 'wilson')`
   
   Precision 
   - `tps / (tps + fps)`
   - `binom.confint(tps,(tps + fps), conf.level = 0.95, method = 'wilson')`
 
   F1 
   - `(2 * precision * recall)/(precision + recall)`
   - ` (2 * upper_precision * upper_recall)/(upper_precision + upper_recall)` & `(2 * lower_precision * lower_recall)/(lower_precision + lower_recall)`
  ## 5. BED File addition
  If a bed file is included, analysis will FIRST be performed on on/off target results. Once off target analysis is complete, the off-target variants will beremoved from further analysis. This allows the user to have an overview of the differencecs with on-off target, and then get an in depth analysis of the on-target results. 
# Genotyping Wiki
If you run the --fillouts flag, there are some additional exceptions and use case assumptions to be made. All the general analysis will be run, but simultaneously you will be running fillouts. Once fillouts completes, the resulting MAFs will reenter the script. 
   ## 1. T_var_freq analysis
   We will use the CALLED variant frequency as the t_var_freq_bin value. When we run fillouts, we must include  chromosome, sample, start, end, **REF AND ALT ** to get any results. Called variant frequency is assumped to be correct. 
   ## 2. No Analysis of Permissive Tag Type
   Since genotyping is dependent on chromosome, sample, start, end, **REF AND ALT ** we cannot say definitively that the filled out locations will make sense with a permissive (Chrom, start, sample). therefore we do not run permissive tag type.
   ## 3.  Must have matching samples in test and ground cohorts<span id="place1">
   Your mapping files for both TEST and GROUND **must** have the same cohort. Your mapping files contain all the BAMs which you want to run fillouts on. If a sample exists in one cohort mapping file but not the other, the tool will fail. 
   ## 4. Images are combined
   Your resulting plots will only generate once you have your fillouts results completed. You will be able to see the differences in results from both the original called results and the genotyped results. 
   ## 5. General Math (definitions change)
   - A variant is DETECTABLE if `t_total_count >= 20`
   - A variant has EVIDENCE if `t_alt_count >= 1` 
   True Positive (tps): variant has EVIDENCE in both Test and Ground MAF. Detectablity status DOES NOT matter. 
   
   False Positive (fps): variant has EVIDENCE in ground, but is DETECTABLE with NO EVIDENCE in test.
   
   False Negative (fns): variant has EVIDENCE in Test, but is DETECTABLE with NO EVIDENCE in ground.
   
   Ground_set_no_ev_not_detectable: variant has no EVIDENCE and not DETECTABLE in ground
   
   Test_set_no_ev_not_detectable: variant has no EVIDENCE and not DETECTABLE in test

   Vars_with_no_evidence_in_either_test_or_ground: variant has no EVIDENCE and not DETECTABLE in test. These variants have been CALLEd by either test or ground, but genotyping did not find it. 
  
