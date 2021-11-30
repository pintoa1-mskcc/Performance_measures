############  REQUIRED LIBRARIES ############ 
suppressPackageStartupMessages({library(dplyr)
  library(data.table)
  library(stringr) 
  library(jsonlite)
  library(binom)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(tidyverse)
  library(tidyr)
  library(argparse) 
  library(doParallel)})
library(bedr)
############################################
doParallel::registerDoParallel(cores = 4)

############  REQUIRED FUNCTIONS ############ 
set.seed(123) 
oncokb <- fromJSON(readLines('http://oncokb.org/api/v1/genes', warn=F))
source("performance_measure_custom_functions.R")

############################################

############  READ IN ARGUMENTS ############ 
opt = commandArgs(TRUE)

parser=ArgumentParser()
parser$add_argument("-g", "--ground", type='character', default=NULL,
                    help="ground truth MAF file ")
parser$add_argument("-t", "--test", type='character',  default=NULL,
                    help="test MAF file ")
parser$add_argument('-d','--directory', type='character', default = getwd(),
                    help='output directory [default = current working directory ]')
parser$add_argument('-v','--additional_variables', type='character',  default = NULL,
                    help= 'additional columns from MAFs which you would like to run recall and precision. This is a single value or a comma separated list if there is more than one')
parser$add_argument("-o", "--out_prefix", type="character", default="date()",
              help="output file name basename [default= %default]")
parser$add_argument('-f','--fillouts', type='logical',
                    default = FALSE,
                    help = 'Logical stating whether or not fillouts needs to be performed ')
parser$add_argument('-e','--test_fillout_mapping', type='character',
                    default = NULL,
                    help = 'File containing BAM tumor/normal mapping for test set')
parser$add_argument('-r','--ground_fillout_mapping', type='character',
                    default = NULL, 
                    help = 'File containing BAM tumor/normal mapping for ground set')
parser$add_argument('-b','--bed_file', type = 'character', default = NULL, help='BED file representing the intersection of the BED files used to generate ground and truth MAFs. If provided, tool will return a MAF for ground and test with each unique identifier and additional columns')
parser$add_argument('-p', '--fillout_to_pr', type = 'logical',default = FALSE, help ='Logical stating whether or not fillouts has ALREADY  been performed')
parser$add_argument('-c','--called_directory', type = 'character', default = NULL, help = 'Location of performance measure results on CALLED mutations (not genotyped). If provided, will generated statistics graphs for the combined results')
parser$add_argument('-u','--called_out_prefix', type = 'character', default = NULL, help = 'Out prefix for performance measure called results If not provided, assumes the out_prefix provided is of form "fillout_%called_out_prefix%"')
opt=parser$parse_args()


test <- fread(opt$test,data.table = FALSE)
ground <- fread(opt$ground,data.table = FALSE)

if(!is.null(opt$bed)){
  bed <- fread(opt$bed, data.table = FALSE)
}

if(opt$out_prefix == "date()"){
  out_prefix <- str_replace_all(date()," ","_")

} else {
  out_prefix <- opt$out_prefix
}

if(!is.null(opt$called_directory)){
  if(is.null(opt$called_out_prefix)){
    opt$called_out_prefix <- str_replace(out_prefix,"fillout_","")
  }
}
write(paste0("Run name: ",out_prefix),stderr())


write(paste0("Ground file: ", opt$ground),stderr())
write(paste0("Test file: ", opt$test),stderr())


directory <- paste0(opt$directory,'/')
dir.create(directory)
dir.create(paste0(directory,'pdfs/'))
write(paste0("Output directory: ",directory),stderr())

if(!is.null(opt$additional_variables) ) {
  if(grepl(',',opt$additional_variables) ){
  additional_variables <- str_split(opt$additional_variables,',')
  } else {
    additional_variables <- opt$additional_variables
  }
  for (variable in additional_variables){
    if (variable %nin% colnames(test) | variable %nin% colnames(ground)){
      warning("A provided additional variable to parse does not exist in provided MAFs. This variable will be removed from analysis")
      additional_variables <- additional_variables[additional_variables != variable]
    }
  }
} else  {
  additional_variables <- NULL
}



############################################



############################################

##### REFORMAT MAF FOR EXPECTED INPUT FOR FUNCTIONS ############################################


### STEP 1: Check that the samples exist in both files (added one sample to test set for testing purposed)
missing_in_ground <- unique(test$Tumor_Sample_Barcode) %nin% unique(ground$Tumor_Sample_Barcode)
missing_in_test <- unique(ground$Tumor_Sample_Barcode) %nin% unique(test$Tumor_Sample_Barcode)

if(any(c(missing_in_test,missing_in_ground))) {
  warning("Sample(s) exists in one MAF but not in the other. These samples are retained in analysis but will have NA or zero values for either recall or precision. See 000.txt for more information")
  names_m_i_g <-unique(test$Tumor_Sample_Barcode)[missing_in_ground]
  names_m_i_t <- unique(test$Tumor_Sample_Barcode)[missing_in_test]
  if(length(names_m_i_g) > 0 ) {
    names(names_m_i_g) <- 'Missing_In_Ground_Samples' 
    }
  if(length(names_m_i_t) > 0 ) {
    names(names_m_i_t) <- 'Missing_In_Test_Samples' 
  }
  warning_return_1 <- c(names_m_i_g,names_m_i_t)
  write.table(warning_return_1,file=paste0(directory,'000.txt'),quote = FALSE)
}

#Keep test tumor samples which are not misisngi n ground 
all_samples <- unique(c(test$Tumor_Sample_Barcode, ground$Tumor_Sample_Barcode))
names(all_samples) <- all_samples

# If a BED is provided, filter all variants so only on target is analyzed 
if(!is.null(opt$bed)){
  test <- test %>% mutate(bed_tag = paste(Chromosome,Start_Position,End_Position, sep = '_'))
  ground <- ground %>% mutate(bed_tag = paste(Chromosome,Start_Position,End_Position, sep = '_'))
  
  test_bed <- test[,c('Chromosome','Start_Position','End_Position','bed_tag')] 
  colnames(test_bed) <- c('chr','start','end','names')
  
  ### BEDR CURRENTLY HAS A BUG FOR R VERSION 4.0 and BEYOND. IT WILL SAVE THE OUTPUT FILE FINE, BUT IS INCAPABLE OF SAVING AS AN R OBJECT (appears to attempt to concatinate colnames for bed and test_bed, then returns only 4 columns. erroring colnames() attempt. )
  try(bedr(engine = 'bedtools', input = list(a = test_bed, b = bed), method = "intersect", params = ' -wa ', check.chr = TRUE,
                   check.zero.based = FALSE,
                   check.valid = FALSE,
                   check.sort = FALSE,
                   check.merge = FALSE,outputFile =paste0(directory,'/test.bed'),verbose = TRUE))
  bed_test <- fread(paste0(directory,'/test.bed'),data.table = FALSE)
  
  ground_bed <- ground[,c('Chromosome','Start_Position','End_Position','bed_tag')] 
  colnames(ground_bed) <- c('chr','start','end','names')
  try(bedr(engine = 'bedtools', input = list(a = ground_bed, b = bed), method = "intersect", params = ' -wa ', check.chr = TRUE,
            check.zero.based = FALSE,
            check.valid = FALSE,
            check.sort = FALSE,
            check.merge = FALSE,outputFile = paste0(directory,'/ground.bed'),verbose = TRUE))
  
  bed_ground <- fread(paste0(directory,'/ground.bed'),data.table = FALSE)
  
  ground <- ground %>% mutate(on_target = ifelse(bed_tag %in% bed_ground$V4, TRUE,FALSE))
  test <- test %>% mutate(on_target = ifelse(bed_tag %in% bed_test$V4, TRUE,FALSE))
}



### Make tag ids, var_tag is restrictive, TAG is permissive mode
ground<- ground %>% mutate(var_tag = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2,':',Tumor_Sample_Barcode),
                           TAG = str_c(Chromosome,':',Start_Position,':',Tumor_Sample_Barcode)) %>%
  mutate(oncogenic_tf =ifelse(grepl("ncogenic", oncogenic),'ONCOGENIC' ,'OTHER')) %>%
  mutate(clonality =  ifelse(is.na(cf) | cf < (0.6 * purity), 
                            'INDETERMINATE',
                                             ifelse((ccf_expected_copies > 0.8 | (ccf_expected_copies > 0.7 & ccf_expected_copies_upper > 0.9)), 
                                                                                                                        'CLONAL', 
                                                                                                                        'SUBCLONAL')))%>%
  mutate(is_non_syn_mut = ifelse(Variant_Classification %in% c("Missense_Mutation", 
                                                               "Nonsense_Mutation", 
                                                               "Nonstop_Mutation", 
                                                               "Frame_Shift_Ins", 
                                                               "Frame_Shift_Del",
                                                               "In_Frame_Del",
                                                               "In_Frame_Ins",
                                                               "Translation_Start_Site",
                                                               "Splice_Site"), T, F)) 

### We will calcaulate clonality for test var tags in case they are false positives, but true positives must match between ground and test var tags
test <- test  %>% mutate(var_tag = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2,':',Tumor_Sample_Barcode), 
                         TAG = str_c(Chromosome,':',Start_Position,':',Tumor_Sample_Barcode)) %>% 
  mutate(oncogenic_tf = ifelse(grepl("ncogenic", oncogenic),'ONCOGENIC' ,'OTHER') ) %>%
  mutate(clonality =  ifelse(is.na(cf) | cf < (0.6 * purity), 
                             'INDETERMINATE',
                                             ifelse((ccf_expected_copies > 0.8 | (ccf_expected_copies > 0.7 & ccf_expected_copies_upper > 0.9)), 
                                                                                                                          'CLONAL', 
                                                                                                                          'SUBCLONAL'))) %>%
  mutate(is_non_syn_mut = ifelse(Variant_Classification %in% c("Missense_Mutation", 
                                                               "Nonsense_Mutation", 
                                                               "Nonstop_Mutation", 
                                                               "Frame_Shift_Ins", 
                                                               "Frame_Shift_Del",
                                                               "In_Frame_Del",
                                                               "In_Frame_Ins",
                                                               "Translation_Start_Site",
                                                               "Splice_Site"), T, F)) 
# Clonality must match to the ground values for accurate calculation of recall
shared_variants <- test$var_tag[test$var_tag %in% ground$var_tag]
test[match(shared_variants,test$var_tag),'clonality'] <- ground[match(shared_variants, ground$var_tag),'clonality']


test <- test %>% mutate(ref_to_alt = paste(Reference_Allele,Tumor_Seq_Allele2,sep=">"))  %>%
  mutate(substitutions = ifelse( Variant_Type %nin% c('SNP','SNV'), NA,sapply( ref_to_alt,function(ref_to_alt ) {
    new_sub <- switch(ref_to_alt, "A>C" = "T>G", "T>G" = "T>G","A>G" ="T>C","T>C" = "T>C","A>T" = "T>A","T>A" ="T>A",
                      "C>A" ="C>A" , "G>T" = "C>A" , "C>G" = "C>G","G>C" = "C>G", "C>T" = "C>T","G>A" = "C>T")
    return(new_sub) 
  }))) %>%  mutate_cond(Variant_Type == 'SNP', Variant_Type = 'SNV')

ground <- ground %>% mutate(ref_to_alt = paste(Reference_Allele,Tumor_Seq_Allele2,sep=">")) %>%
  mutate(substitutions = ifelse( Variant_Type %nin% c('SNP','SNV'), NA,sapply( ref_to_alt,function(ref_to_alt ) {
    new_sub <- switch(ref_to_alt, "A>C" = "T>G", "T>G" = "T>G","A>G" ="T>C","T>C" = "T>C","A>T" = "T>A","T>A" ="T>A",
                      "C>A" ="C>A" , "G>T" = "C>A" , "C>G" = "C>G","G>C" = "C>G", "C>T" = "C>T","G>A" = "C>T")
    return(new_sub) 
  }))) %>%  mutate_cond(Variant_Type == 'SNP', Variant_Type = 'SNV')


# Add purity and t_vaf_freq bin columns for parsing
#specify interval/bin labels
breaks <- seq(0,1,0.05)
tags <- c("[0-5)","[5-10)", "[10-15)", "[15-20)", "[20-25)", "[25-30)","[30-35)", "[35-40)","[40-45)", "[45-50)","[50-55)","[55-60)","[60-65)","[65-70)","[70-75)","[75-80)","[80-85)","[85-90)","[90-95)","[95-100)")
tags <- factor(tags,levels=tags)


ground$t_var_freq_bin <- cut(ground$t_var_freq, 
                             breaks=breaks, 
                             include.lowest=TRUE, 
                             right=FALSE, 
                             labels=tags)



test$t_var_freq_bin <- cut(test$t_var_freq, 
                           breaks=breaks, 
                           include.lowest=TRUE, 
                           right=FALSE, 
                           labels=tags)


# idnetify shared variants which have differing 't_var_freq_bin'
# For purposes of recall, assume that ground t_var_freq_bins are the 'correct' variables so set those test values to that
test[match(shared_variants,test$var_tag),'t_var_freq_bin'] <- ground[match(shared_variants, ground$var_tag),'t_var_freq_bin']



# Since we have two MAFs there is potential that there are two different purities between MAFs for the same samples
## To analyze recall on purity, we assume that the ground files purity is the purity which we wish to calculate recall over
# Therefore we create a new variable to get the binned samples in the same bucket for comparison

if(length(unique(ground$purity)) != 0){
  ground$purity_bin <- cut(ground$purity, 
                           breaks=breaks, 
                           include.lowest=TRUE, 
                           right=FALSE, 
                           labels=tags)
  
  
  tumor_sample_purity_mapping <- ground %>% distinct(Tumor_Sample_Barcode, purity_bin) 
  test <- left_join(test[,colnames(test) %nin% c('purity_bin')],tumor_sample_purity_mapping,by = "Tumor_Sample_Barcode")
}

#### THIS SCRIPT UTILIZES n_variant_frequency AS AN INDICATOR THAT FILLOUTS HAS BEEN RUN, If fillouts has been run, performance measures are only run on detectable reads
if(opt$fillout_to_pr){
  test <- test %>% mutate(evidence = ifelse(t_alt_count >= 1, TRUE, FALSE))
  test <- test %>% mutate(detectable = ifelse(t_total_count >= 20, TRUE, FALSE))
  
  ground <- ground %>% mutate(evidence = ifelse(t_alt_count >= 1, TRUE, FALSE))
  ground <- ground %>% mutate(detectable = ifelse(t_total_count >= 20, TRUE, FALSE))

}
# Formatting
i <- sapply(test, is.factor)
test[i] <- lapply(test[i], as.character)

i <- sapply(ground, is.factor)
ground[i] <- lapply(ground[i], as.character)

### Clarity NA as a character NA for analysis
test[is.na(test$clonality),'clonality'] <- 'N/A'
ground[is.na(ground$clonality),'clonality'] <- 'N/A'
test[is.na(test$purity_bin),'purity_bin'] <- 'N/A'
ground[is.na(ground$purity_bin),'purity_bin'] <- 'N/A'


############################################

######### IF FILLOUTS RETURN FILLOUTS MAFS ############################
if( opt$fillouts){

  fillout_maf <- rbind(ground,test)
  fillout_maf <- fillout_maf[!duplicated(fillout_maf$var_tag),]
  fillout_maf <- as.data.frame(unnest(fillout_maf, substitutions))
  
  fillout_mapping_test <- read.table(opt$test_fillout_mapping, header = TRUE, stringsAsFactors = FALSE)
  fillout_mapping_ground <- read.table(opt$ground_fillout_mapping, header = TRUE, stringsAsFactors = FALSE)
  
  fillout_output_dir <- paste0(opt$directory,'/', 'fillout/')
  fillout_combined_mafs <- paste0(opt$directory,'/', 'fillout/target_mafs/')
  fillout_results_dir <- paste0(opt$directory,'/', 'fillout/result_mafs/')
  
  all_provided_bams_tumor_ids <- unique(c(fillout_mapping_test$tumor_id,fillout_mapping_ground$tumor_id))
  if(any(all_provided_bams_tumor_ids %nin% fillout_maf$Tumor_Sample_Barcode)){
    warning('A provided Tumor BAM for fillouts does not exist within the provided MAFs. This sample will not have fillouts run as it has no mutations. See 001.txt for the sample id')
    removed_bam_sample <- all_provided_bams_tumor_ids[all_provided_bams_tumor_ids %nin% fillout_maf$Tumor_Sample_Barcode]
    write.table(removed_bam_sample,file=paste0(directory,'001.txt'),quote = FALSE) 
  }
  if (any( fillout_maf$Tumor_Sample_Barcode %nin% fillout_mapping_test$tumor_id) | any(fillout_maf$Tumor_Sample_Barcode %nin% fillout_mapping_ground$tumor_id)) {
    stop('A sample within the P/R cohort is missing from the provided fillout BAMs. Please rerun this script again with the correct bam mapping. ')
  }
  
  dir.create(fillout_output_dir)
  dir.create(fillout_combined_mafs)
  dir.create(fillout_results_dir)
  dir.create(paste0(fillout_results_dir,'logs/'))
  
  dir.create(paste0(fillout_results_dir,'ground/'))
  dir.create(paste0(fillout_results_dir,'test/'))
  

  fillout_commands <-  function(sample) {
    sample_fillout <- fillout_maf[fillout_maf$Tumor_Sample_Barcode == sample,]
    test_tumor_bam <- fillout_mapping_test$tumor_bam[fillout_mapping_test$tumor_id == sample]
    test_normal_bam <- fillout_mapping_test$normal_bam[fillout_mapping_test$tumor_id == sample]
    test_dir_norm <- dirname(fillout_mapping_test$normal_bam[fillout_mapping_test$tumor_id == sample])
    test_dir_tumor <- dirname(fillout_mapping_test$tumor_bam[fillout_mapping_test$tumor_id == sample])
    
    
    
    ground_tumor_bam <- fillout_mapping_ground$tumor_bam[fillout_mapping_ground$tumor_id == sample]
    ground_normal_bam <- fillout_mapping_ground$normal_bam[fillout_mapping_ground$tumor_id == sample]
    ground_dir_norm <- dirname(fillout_mapping_ground$normal_bam[fillout_mapping_ground$tumor_id == sample])
    ground_dir_tumor <- dirname(fillout_mapping_ground$tumor_bam[fillout_mapping_ground$tumor_id == sample])
    
    
    job_name <- paste0('fillout_',out_prefix,sample)
    sample_maf <- paste0(fillout_combined_mafs,sample,'_UNIFIED_GROUND_TEST.maf')
    write.table(sample_fillout,file=sample_maf, row.names=FALSE,quote=FALSE, sep= '\t')
    system("module load singularity")
    system("module load R/R-4.0.5")
    test_fillout_command <- paste0('bsub -J ',job_name,'_test -e ',fillout_output_dir,'/logs/ -n 2 -R rusage[mem=5] -We 0:59 " /juno/work/ccs/pintoa1/fillout_testing/maf_fillout.py -m ', sample_maf, ' -b ',test_tumor_bam,' ', test_normal_bam,
                                   ' -o ',fillout_results_dir,'test/test',sample,'_fillout.maf -mo"')
    
  
    ground_fillout_command <- paste0('bsub -J ',job_name,'_ground -e ',fillout_output_dir,'/logs/ -n 2 -R rusage[mem=5] -We 0:59  " /juno/work/ccs/pintoa1/fillout_testing/maf_fillout.py -m ', sample_maf, ' -b ',ground_tumor_bam,' ', ground_normal_bam,
                                     ' -o ',fillout_results_dir,'ground/ground_',sample,'_fillout.maf -mo"')
  write(ground_fillout_command,stderr())
  system(test_fillout_command)
     
      system(ground_fillout_command)
   
  return(c('job_name_test' = paste0(job_name,'_test'), 'job_name_ground' = paste0(job_name,'_ground')))
  }

  write(paste0("Submitting: ", length(all_samples)*2, " jobs."),stderr())
  
  queued_jobs <- plyr::adply(all_samples, 1, fillout_commands, .parallel = T)
  
  write(paste0("Queued: ", dim(queued_jobs)[1]*2, " jobs."),stderr())
  
}
############################################

#####################RUN STATISTICS AND SAVE OUTPUTS####################
variables_to_parse <- c('oncogenic_tf','clonality','substitutions','is_non_syn_mut','t_var_freq_bin',additional_variables,'purity_bin')


binned_variables <- variables_to_parse[grepl('bin',variables_to_parse)]


if(!is.null(opt$bed_file)){
  df <- parse_dataframe_on_var(ground,test,'on_target','cohort')
  df <- df[,c('permission','type','on_target','statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]
  
  write.table(df,paste0(directory,out_prefix,'_','on_target','_cohort_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')
  df1 <- df %>% filter(permission == 'restrictive')
  statistics_graphs(df1,'on_target','bar',directory,out_prefix)
  
  sample_level_df <- lapply(all_samples, function(sample){
    sample_ground <- ground[ground$Tumor_Sample_Barcode == sample,]
    sample_test <- test[test$Tumor_Sample_Barcode == sample ,]
    
    sample_level_stats <- parse_dataframe_on_var(sample_ground,sample_test,'on_target','sample')
    if(nrow(sample_level_stats) != 0){
      sample_level_stats$Tumor_Sample_Barcode <- sample
    }
    return(sample_level_stats)
  })
  
  sample_level_df <- do.call(rbind,sample_level_df)
  cols.num <- c(seq(2,6,1),seq(8,10,1))
  sample_level_df[cols.num] <- lapply(sample_level_df[cols.num], as.numeric)
  sample_level_df <- sample_level_df[,c('permission','type','Tumor_Sample_Barcode','on_target','statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]
  
  
  write.table(sample_level_df,paste0(directory,out_prefix,'_','on_target','_sample_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')
  df1 <- sample_level_df %>% filter(permission == 'restrictive')
  statistics_graphs(df1,'on_target','boxplot',directory,out_prefix)
  
  
  ### ONLY PARSING THE ON TARGET VALUES FROM HERE ON
  ground_off <- ground %>% filter(on_target == FALSE)
  ground <- ground %>% filter(on_target == TRUE)
  test_off <- test %>% filter(on_target == FALSE)
  test <- test %>% filter(on_target == TRUE)
}

overview_df <- calc_stats_by_variant_type(ground,test,'cohort')
overview_df <- overview_df[,c('permission','type','statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]
write.table(overview_df,paste0(directory,out_prefix,'_overview_all_variants_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')
if(!is.null(opt$called_directory)){
  c_df <- read.table(paste0(opt$called_directory,opt$called_out_prefix,'_overview_all_variants_performance_measures.txt'),header = TRUE)
  statistics_graphs(overview_df[overview_df$permission == 'restrictive',],'type','bar',directory,paste0('combined_',out_prefix),c_df[c_df$permission=='restrictive',])
} 
statistics_graphs(overview_df[overview_df$permission == 'restrictive',],'type','bar',directory,out_prefix)




variable_parsing_and_graph <- function(variable) {
  df <- parse_dataframe_on_var(ground,test,variable,'cohort')
  df <- df[,c('permission','type',variable,'statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]
  
  write.table(df,paste0(directory,out_prefix,'_',variable,'_cohort_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')
  if(variable %nin% binned_variables){
    df1 <- df %>% filter(permission == 'restrictive')
    if(!is.null(opt$called_directory)){
   
      
      
    } else {
      statistics_graphs(df1,variable,'bar',directory,out_prefix)
    }
    
    return(NULL)
  } else{
    assign(paste0(variable,'_df'),df)
    return(df)
  }
}

binned_vars <- plyr::adply(variables_to_parse, 1, variable_parsing_and_graph, .parallel = T)



############################################

######## SPECIAL BINNED PLOTS ########
binned_vars <- binned_vars %>% filter(type == 'all') %>% filter(permission == 'restrictive')
purity_nas <- binned_vars[grepl('^N',binned_vars$purity_bin) & !is.na(binned_vars$purity_bin),]
purity_nas$Variable_ID <- 'purity'

binned_vars <- binned_vars %>% mutate(Variable_ID =  factor(ifelse(!is.na(t_var_freq_bin),'vaf',ifelse(!is.na(purity_bin),'purity',NA)) )) %>%
  mutate(Frequency = factor(ifelse(!is.na(t_var_freq_bin),t_var_freq_bin,ifelse(!is.na(purity_bin),purity_bin,NA)),levels=tags)) %>% filter(!is.na(Frequency))
# 

vaf_mut_count <- ggplot(binned_vars[binned_vars$Variable_ID == 'vaf',], aes(x = Frequency, y = total_var_count)) +
  geom_bar(stat='identity',position=position_dodge(),width=0.75) + theme_classic() + ylab('N Mutations') +
  theme(axis.text.x = element_blank(),legend.position = "top",axis.title.x = element_blank(),legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = levels(binned_vars[binned_vars$Variable_ID == 'vaf','Frequency']),drop=FALSE) + annotate('text', label='VAF', x=Inf, y=Inf, hjust=1, vjust=1)+ labs(x=NULL)
recall_bin <- ggplot(binned_vars[binned_vars$statistic_name == 'recall',] , aes(x = Frequency, y = value,group = Variable_ID, color = Variable_ID)) + geom_line(stat='summary')  +  scale_color_jama() +
  geom_errorbar(aes(ymin =lower, ymax = upper), width=0) + theme_classic() + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),legend.position = "none",legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = levels(binned_vars[,'Frequency']),drop = FALSE) + ylim(0,1) + labs(x=NULL)+ ylab('Recall')
precision_bin <- ggplot(binned_vars[binned_vars$statistic_name == 'precision',] , aes(x = Frequency, y = value,group = Variable_ID, color = Variable_ID)) + geom_line(stat='summary')  +  scale_color_jama() +
  geom_errorbar(aes(ymin =lower, ymax = upper), width=0) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none",legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = levels(binned_vars[,'Frequency']),drop = FALSE) + ylim(0,1) + ylab('Precision')
pur_sample_count <- ggplot(binned_vars[binned_vars$Variable_ID == 'purity',], aes(x = Frequency, y = n_samples)) +
      geom_bar(stat='identity',position=position_dodge(),width=0.75) + theme_classic() +
      theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.position = "none",legend.background=element_blank(),legend.title=element_blank()) +
      scale_x_discrete(labels = levels(binned_vars[,'Frequency']),drop = FALSE) +annotate('text', label='Purity', x=-Inf, y=Inf, hjust=0, vjust=1) + ylab("N Samples")
 

na_variant <- ggplot(purity_nas[purity_nas$Variable_ID == 'vaf',], aes(x = purity_bin, y = total_var_count)) +
  geom_bar(stat='identity',position=position_dodge(),width=0.75) + theme_classic() + ylab('N Mutations') +
  theme(axis.text.x = element_blank(),legend.position = "top",axis.title.x = element_blank(),legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = c('N/A'),drop=FALSE) + annotate('text', label='VAF', x=-Inf, y=Inf, hjust=1, vjust=1)+ labs(x=NULL) + ylim(0,1)

na_sample_count <- ggplot(purity_nas,aes(x = purity_bin, y = n_samples)) + geom_bar(stat='identity',position=position_dodge(),width=0.75) + theme_classic() +
  theme(axis.title.x = element_blank(),legend.position = "top",legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = levels(purity_nas[,'purity_bin']),drop=FALSE) + annotate('text', label='Purity', x=-Inf, y=Inf, hjust=0, vjust=1) + ylab("N Samples")
na_recall_bin <- ggplot(purity_nas[purity_nas$statistic_name == 'recall',] , aes(x = purity_bin, y = value,group = Variable_ID, color = Variable_ID)) + geom_point()  +  scale_color_jama() +
  geom_errorbar(aes(ymin =lower, ymax = upper, width = 0.1)) + theme_classic() + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),legend.position = "none",legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = levels(purity_nas[,'purity_bin'])) + ylim(0,1) +  ylab('Recall')
na_precision_bin <- ggplot(purity_nas[purity_nas$statistic_name == 'precision',] , aes(x = purity_bin, y = value,group = Variable_ID, color = Variable_ID)) + geom_point()  +  scale_color_jama() +
  geom_errorbar(aes(ymin =lower, ymax = upper, width = 0.1)) + theme_classic() +theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "bottom",legend.background=element_blank(),legend.title=element_blank()) +
  scale_x_discrete(labels = unique(purity_nas[,'purity_bin'])) + ylim(0,1) +  ylab('Precision') + xlab("Frequency")





pdf(paste0(directory,'pdfs/',out_prefix,'_binned_vars','.pdf'))
ggarrange(ggarrange(na_variant,na_sample_count,na_recall_bin,na_precision_bin, ncol= 1, align ='hv',common.legend = TRUE, legend = "bottom"),ggarrange(vaf_mut_count,pur_sample_count,recall_bin,precision_bin, ncol = 1,common.legend = TRUE, legend = "bottom",align = "hv"),ncol = 2,align = "hv", widths = c(0.13,1.0),common.legend = TRUE, legend = "bottom")

    print(ggarrange(vaf_mut_count,pur_sample_count,recall_bin,precision_bin, ncol = 1,common.legend = TRUE, legend = "bottom",align = "hv"))

dev.off()
# ############################################

# ########## SAMPLE LEVEL PLOTS #############

sample_level_raw <- lapply(all_samples, function(sample){
  sample_ground <- ground[ground$Tumor_Sample_Barcode == sample,]
  sample_test <- test[test$Tumor_Sample_Barcode == sample ,]
  
  sample_level_stats <- calc_stats_by_variant_type(sample_ground,sample_test,'sample')
  sample_level_stats$Tumor_Sample_Barcode <- sample
  return(sample_level_stats)
})

sample_level_raw <- do.call(rbind,sample_level_raw)
cols.nums <- c(seq(2,8,1),seq(10,12,1))
sample_level_raw[cols.nums] <- lapply(sample_level_raw[cols.nums], as.numeric)
sample_level_raw <- sample_level_raw[,c('permission','type','Tumor_Sample_Barcode','statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]

write.table(sample_level_raw,paste0(directory,out_prefix,'_sample_overview_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')

statistics_graphs(sample_level_raw[sample_level_raw$permission == 'restrictive',],'type','boxplot',directory,out_prefix)



variable_parsing_and_graph_sample <- function(variable) {
  sample_level_df <- lapply(all_samples, function(sample){
    sample_ground <- ground[ground$Tumor_Sample_Barcode == sample,]
    sample_test <- test[test$Tumor_Sample_Barcode == sample ,]
    
    sample_level_stats <- parse_dataframe_on_var(sample_ground,sample_test,variable,'sample')
    if(nrow(sample_level_stats) != 0){
      sample_level_stats$Tumor_Sample_Barcode <- sample
    }
    return(sample_level_stats)
  })
  
  sample_level_df <- do.call(rbind,sample_level_df)
  cols.nums <- c(seq(2,8,1),seq(10,12,1))
  sample_level_df[cols.nums] <- lapply(sample_level_df[cols.nums], as.numeric)
  sample_level_df <- sample_level_df[,c('permission','type','Tumor_Sample_Barcode',variable,'statistic_name','value','lower','upper','total_var_count','n_samples','tps','fps','fns','truth_set_no_ev_not_detect','test_set_no_ev_not_detect')]

  write.table(sample_level_df,paste0(directory,out_prefix,'_',variable,'_sample_performance_measures.txt'),quote = FALSE,row.names = FALSE,sep = '\t')
  
  if(variable %nin% binned_variables){
    df1 <- sample_level_df %>% filter(permission == 'restrictive')
    statistics_graphs(df1,variable,'boxplot',directory,out_prefix)
  } 
  return(NULL)
}


sample_variables_to_parse <- c('oncogenic_tf','clonality','substitutions','is_non_syn_mut','t_var_freq_bin', additional_variables)

returning_null <- plyr::adply(sample_variables_to_parse, 1, variable_parsing_and_graph_sample, .parallel = T)

############################################

##### PR CURVE
pr_curve_df <- left_join(sample_level_raw[sample_level_raw$type == 'all',],ground %>% select(c(Tumor_Sample_Barcode, purity)) %>% distinct(), by = "Tumor_Sample_Barcode")

pr_curve_df <- pr_curve_df %>% filter(statistic_name %in% c('recall','precision')) %>% filter(permission == 'restrictive') %>%
  select(c(statistic_name,value,Tumor_Sample_Barcode,purity)) %>% distinct(statistic_name,value,Tumor_Sample_Barcode,purity) %>%
  pivot_wider(id_cols= c(Tumor_Sample_Barcode, purity),names_from = statistic_name, values_from =  value) 
  
pdf(paste0(directory,'pdfs/',out_prefix,'_precision_recall_plot','.pdf'))
  ggplot(pr_curve_df,aes(x = recall, y = precision, color = purity)) + geom_point()  + scale_color_gradient() +
    theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position = "right",legend.background=element_blank()) +
    ggtitle('Precision Recall Per Sample') + ylim(0,1) + xlim(0,1)
dev.off()
############################################

###### If Fillouts, check status of fillouts and wait for completion 
is_job_still_running <- function(job_name) {
  return (any(grepl("not found", system(paste0('bjobs -u all -J ', job_name, ' 2>&1'), intern=T))))
}

if (opt$fillouts){
  #####
  ##### Wait for job completion
  #####
  completed_jobs <- c('none')
  still_running_jobs <- c(queued_jobs$job_name_test[queued_jobs$job_name_test %nin% completed_jobs],queued_jobs$job_name_ground[queued_jobs$job_name_ground %nin% completed_jobs])
  while(1) {
    
    completed_jobs <- c(completed_jobs, sapply(still_running_jobs, function(job) {
      ifelse(is_job_still_running(job), return(job),return(NULL)) 
      }))
    still_running_jobs <- c(queued_jobs$job_name_test[queued_jobs$job_name_test %nin% completed_jobs],queued_jobs$job_name_ground[queued_jobs$job_name_ground %nin% completed_jobs])
    
    write(paste0("# of fillouts still running... ", length(still_running_jobs)),stderr())
    
    if (length(still_running_jobs) == 0) { write("All fillouts have completed!",stderr()); break }  
    
    
    Sys.sleep(120) ## Check for the jobs completion every couple of mins.
  }
  
  ground_dir <- paste0(fillout_results_dir,'ground/')
  test_dir <- paste0(fillout_results_dir,'test/')
  rscript_dir <- getwd()
  if(is.null(opt$bed_file)){
    write(paste0('bsub -J  collect_fillouts_results -e ',fillout_output_dir,' -n 2 -R rusage[mem=5] -We 0:59 "Rscript fillouts_restruct.R -r ', ground_dir, ' -d ', directory, ' -o ',out_prefix,' -p TRUE -e ', test_dir, ' -m ', fillout_combined_mafs, ' -j TRUE"' ),stderr())
    
    system(paste0('bsub -J  collect_fillouts_results -e ',fillout_output_dir,' -n 2 -R rusage[mem=5] -We 0:59 "Rscript fillouts_restruct.R -r ', ground_dir, ' -d ', directory, ' -o ',out_prefix,' -p TRUE -e ', test_dir, ' -m ', fillout_combined_mafs, ' -j TRUE"' ))
  } else{
    system(paste0('bsub -J  collect_fillouts_results -e ',fillout_output_dir,' -n 2 -R rusage[mem=5] -We 0:59 "Rscript fillouts_restruct.R -r ', ground_dir, ' -d ', directory, ' -o ',out_prefix,' -p TRUE -e ', test_dir, ' -m ', fillout_combined_mafs, ' -b ',opt$bed_file ,' -j TRUE"' ))
    
  }
  
  #
} else{
  if(!is.null(opt$bed_file)){
    ground <- rbind(ground,ground_off)
    test <- rbind(test,test_off)
  }
  if(opt$fillout_to_pr){
  
    
    test <- test %>% mutate(Detectable_in_Other_Run = ifelse(var_tag %in% ground$var_tag[ground$detectable] , TRUE, FALSE))
    ground <- ground %>% mutate(Detectable_in_Other_Run = ifelse(var_tag %in% test$var_tag[test$detectable], TRUE, FALSE))
 
    test <- test %>% mutate(Evidence_in_Other_Run = ifelse(var_tag %in% ground$var_tag[ground$evidence], TRUE, FALSE))
    ground <- ground %>% mutate(Evidence_in_Other_Run = ifelse(var_tag %in% test$var_tag[test$evidence], TRUE, FALSE))
    
    
    
    
  } else{
    test <- test %>% mutate(Called_in_Other_Run = ifelse(var_tag %in% ground$var_tag , TRUE, FALSE))
    ground <- ground %>% mutate(Called_in_Other_Run = ifelse(var_tag %in% test$var_tag , TRUE, FALSE))
  }
  ground <- as.data.frame(unnest(ground, substitutions))
  test <- as.data.frame(unnest(test, substitutions))
  
  
  write.table(ground,paste0(directory,out_prefix,'_ground_PR_labeled.maf'), row.names=FALSE,quote=FALSE, sep= '\t')
  write.table(test,paste0(directory,out_prefix,'_test_PR_labeled.maf'), row.names=FALSE,quote=FALSE, sep= '\t')
  
  
}

