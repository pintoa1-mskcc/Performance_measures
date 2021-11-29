suppressPackageStartupMessages({library(argparse) 
library(data.table)
library(dplyr)
library(stringr)})
#negate
'%nin%' = Negate('%in%')
opt = commandArgs(TRUE)

parser=ArgumentParser()

parser$add_argument('-r','--ground_directory', type='character', default = getwd(), help = 'Ground OR any fillouts directory. Must be ground directory if running performance measures')
parser$add_argument('-d','--directory',type = 'character',default = NULL, help ='Performance Measures / Final MAF Output Directory; default =[%ground_directory/top PR dir]')
parser$add_argument("-o", "--out_prefix" , type = 'character',default = NULL, help = 'Output prefix')
parser$add_argument('-p','--performance_measures',type = 'logical',default = FALSE, help = 'Do you wish to run performance measures on the fillouts result? Must provide a second fillouts directory and ')
parser$add_argument('-e','--test_directory', type='character', default = NULL, help = 'Test fillouts directory, must be provided if running performance measures')
parser$add_argument("-m", "--maf_dir" , type = 'character',default = NULL, help = 'MAF directory for MAFs used in the fillouts script. This is the -m flag in maf_fillouts.py')
parser$add_argument("-j", "--juno" , type = 'logical',default = FALSE)
parser$add_argument("-b","--bed_file", type = "character", default = NULL, help = "If you wish to use a bed file for targetted performance measures.")
parser$add_argument("-s", "--script",type = 'character',default = getwd(),help = 'If running performance measures, expects Rscript to be in working directory If not, please specify.')

opt=parser$parse_args()


opt$script <- paste0(opt$script,'/')
if(is.null(opt$directory)){
  opt$directory <- paste0(opt$ground_directory,'/')
}

if(opt$performance_measures) {
  if(is.null(opt$maf_dir) | is.null(opt$test_directory) ){
    stop('Something went wrong! A unified MAF directory AND a test directory is required if you wish to run performance measure.')
  } 
}



restructure_mafs <- function(x) {
  temp <- fread(x,data.table=FALSE)
  colnam <- colnames(temp)
  temp <- temp  %>% mutate(variant_loc_to_merge = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele1))
  
  ### We will calcaulate clonality for test var tags in case they are false positives, but true positives must match between ground and test var tags
  n <- temp[grepl("_N00|GCT02|uttcc",temp$Tumor_Sample_Barcode),] #2 special old samples with no N in the normal
  t <- temp[!grepl("_N00|GCT02|uttcc",temp$Tumor_Sample_Barcode),]
  colnames(n) <- gsub("^t_","n_",colnames(n))
  
  n$Matched_Norm_Sample_Barcode <- n$Tumor_Sample_Barcode 
  n$Tumor_Sample_Barcode <- NULL
  t$Matched_Norm_Sample_Barcode <- NULL
  n <- n %>% select(c( colnames(n)[grepl('^n',colnames(n))],'Matched_Norm_Sample_Barcode','variant_loc_to_merge'))
  t <- t[,!grepl('^n',colnames(t))]
  
  
  temp2 <- merge(t,n,by='variant_loc_to_merge')
  
  temp2 <- temp2  %>% mutate(var_tag = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele1,':',Tumor_Sample_Barcode), 
                             TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Sample_Barcode))
  temp2$variant_loc_to_merge <- NULL
  
  return(temp2)
}


performance_measures_expected_formatting <- function(fillout_maf,opt){
  fillout_maf$Tumor_Seq_Allele2 <- fillout_maf$Tumor_Seq_Allele1
  fillout_maf$t_var_freq <- fillout_maf$t_variant_frequency
  sample_to_find <- unique(fillout_maf$Tumor_Sample_Barcode)
  target_mafs <- list.files(opt$maf_dir)
  target <- target_mafs[grepl(unique(fillout_maf$Tumor_Sample_Barcode),target_mafs)]
  maf_with_more_info <- fread(paste0(opt$maf_dir,target), data.table = FALSE)

  t <- colnames(fillout_maf)[grepl('^t',colnames(fillout_maf))]
  n <- colnames(fillout_maf)[grepl('^n',colnames(fillout_maf))]
  maf_with_more_info <- maf_with_more_info[, colnames(maf_with_more_info) %nin% c(t,n)]
  fillout_maf <- fillout_maf[,c('var_tag',t,n)]
  fillout_maf$fillout_to_pr <- TRUE
  return(merge(fillout_maf,maf_with_more_info, by = 'var_tag',all.y = FALSE))
  
}

format_output_name <- function(opt,suffix) {
  if(is.null(opt$out_prefix)){
    return(paste0(opt$directory,suffix))
  }else {
   return(paste0(opt$directory,opt$out_prefix,'_',suffix))
  }
}

ground_files <- list.files(opt$ground_directory)
ground_files <- paste0(opt$ground_directory,ground_files)

if(opt$performance_measures) {
  opt$test_directory <- paste0(opt$test_directory,'/')
  opt$maf_dir <-  paste0(opt$maf_dir,'/')
  test_files <- list.files(opt$test_directory)
  test_files <- paste0(opt$test_directory,test_files)
  
  fillex_ground <- do.call(rbind,lapply(ground_files,function(file) {

    output <- restructure_mafs(file) 
    output1 <- performance_measures_expected_formatting(output,opt)
    return(output1)
  } ))
  
  fillex_test <- do.call(rbind,lapply(test_files,function(file) {
    output <- restructure_mafs(file) 
    output1 <- performance_measures_expected_formatting(output,opt)
    return(output1)
  } ))
  
  
  
  
  test_file_name <- format_output_name(opt,'test_mut_somatic_fillout.maf')

  ground_file_name <- format_output_name(opt,'ground_mut_somatic_fillout.maf')
  
  
  write.table(fillex_ground,ground_file_name,quote = FALSE, row.names = FALSE,sep = "\t")
  
  write.table(fillex_test,test_file_name,quote = FALSE, row.names = FALSE,sep = "\t")
  

  if(is.null(opt$bed_file)){
    write(paste0('bsub -e ', opt$directory, ' -n 2 -R "rusage[mem=8]" -W 0:59 "Rscript ',opt$script,'performance_measure_script.R -g ', ground_file_name,' -t ', test_file_name, ' -d ',opt$directory,' -o fillout_', opt$out_prefix,'"'),stderr())
    system(paste0('bsub -e ', opt$directory, ' -n 2 -R "rusage[mem=8]" -W 0:59 "Rscript ',opt$script,'performance_measure_script.R -g ', ground_file_name,' -t ', test_file_name, ' -d ',opt$directory,' -o fillout_', opt$out_prefix,'"'))
  } else {
    write(paste0('bsub -e ', opt$directory, ' -n 2 -R "rusage[mem=8]" -W 0:59 "Rscript ',opt$script,'performance_measure_script.R -g ', ground_file_name,' -t ', test_file_name, ' -d ',opt$directory,' -o fillout_', opt$out_prefix, ' -b ',opt$bed_file ,'"'),stderr())
    system(paste0('bsub -e ', opt$directory, ' -n 2 -R "rusage[mem=8]" -W 0:59 "Rscript ',opt$script,'performance_measure_script.R -g ', ground_file_name,' -t ', test_file_name, ' -d ',opt$directory,' -o fillout_', opt$out_prefix, ' -b ',opt$bed_file ,'"'))
      
  }
  
  
} else { 
  fillex<- do.call(rbind,lapply(ground_files,restructure_mafs))
  
  file_name <- format_output_name(opt$prefix,'mut_somatic_fillout.maf')
  
  write.table(fillex,file_name,quote = FALSE, row.names = FALSE,sep = "\t")
}

