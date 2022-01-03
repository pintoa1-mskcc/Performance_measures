suppressPackageStartupMessages({library(argparse) 
library(data.table)
library(dplyr)
library(stringr)})
#negate
'%nin%' = Negate('%in%')
opt = commandArgs(TRUE)

parser=ArgumentParser()

parser$add_argument('-r','--ground_directory', type='character', default = getwd(), help = 'Ground OR any fillouts directory. Must be ground directory if running performance measures')
parser$add_argument('-d','--directory',type = 'character',default = NULL, help ='Performance Measures / Final MAF Output Directory; default =[top PR dir]')
parser$add_argument("-o", "--out_prefix" , type = 'character',default = NULL, help = 'Output prefix')
parser$add_argument('-p','--performance_measures',type = 'logical',default = FALSE, help = 'Do you wish to run performance measures on the fillouts result? Must provide a second fillouts directory and ')
parser$add_argument('-e','--test_directory', type='character', default = NULL, help = 'Test fillouts directory, must be provided if running performance measures')
parser$add_argument("-m", "--maf_dir" , type = 'character',default = NULL, help = 'MAF directory for MAFs used in the fillouts script. This is the -m flag in maf_fillouts.py')
parser$add_argument("-j", "--juno" , type = 'logical',default = FALSE)
parser$add_argument("-b","--bed_file", type = "character", default = NULL, help = "If you wish to use a bed file for targetted performance measures.")
parser$add_argument("-s", "--script",type = 'character',default = getwd(),help = 'If running performance measures, expects Rscript to be in working directory If not, please specify directory.')
parser$add_argument('-c','--called_directory',type = "character", default=NULL, help = 'Specify the directory containing the performance measures for CALLED MAF results')
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


performance_measures_expected_formatting <- function(file,opt){
  fillout_maf <- fread(file,data.table=FALSE)
  fillout_maf <- fillout_maf  %>% mutate(var_tag = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele1,':',Tumor_Sample_Barcode), 
                             TAG = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Sample_Barcode))
  
  fillout_maf$Tumor_Seq_Allele2 <- fillout_maf$Tumor_Seq_Allele1
  fillout_maf$genotyped_variant_freq <- fillout_maf$t_variant_frequency
  fillout_maf$t_variant_frequency <- NULL
  sample_to_find <- unique(fillout_maf$Tumor_Sample_Barcode)
  target_mafs <- list.files(opt$maf_dir)
  target <- target_mafs[grepl(unique(fillout_maf$Tumor_Sample_Barcode),target_mafs)]
  maf_with_more_info <- fread(paste0(opt$maf_dir,target), data.table = FALSE)

  t <- colnames(fillout_maf)[grepl('^t',colnames(fillout_maf))]
  n <- colnames(fillout_maf)[grepl('^n',colnames(fillout_maf))]
  maf_with_more_info <- maf_with_more_info[, colnames(maf_with_more_info) %nin% c(t,n)]
  fillout_maf <- fillout_maf[,c('var_tag',t,n,'genotyped_variant_freq')]
  return(merge(fillout_maf,maf_with_more_info, by = 'var_tag',all.y = FALSE))
  
}

format_output_name <- function(opt,suffix) {
  if(is.null(opt$out_prefix)){
    return(paste0(opt$directory,'/fillout/',suffix))
  }else {
   return(paste0(opt$directory,'/fillout/',opt$out_prefix,'_',suffix))
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
    write(file,stderr())
    output <- performance_measures_expected_formatting(file,opt) 
    
    return(output)
  } ))
  
  fillex_test <- do.call(rbind,lapply(test_files,function(file) {
    output <- performance_measures_expected_formatting(file,opt) 
   
    return(output)
  } ))
  
  
  
  
  test_file_name <- format_output_name(opt,'genotyped_test.maf')

  ground_file_name <- format_output_name(opt,'genotyped_ground.maf')
  
  
  write.table(fillex_ground,ground_file_name,quote = FALSE, row.names = FALSE,sep = "\t")
  
  write.table(fillex_test,test_file_name,quote = FALSE, row.names = FALSE,sep = "\t")
  name_test <- basename(opt$test_directory)
  name_ground <- basename(opt$ground_directory)
  bsub_command <- paste0('bsub -e ', opt$directory,'logs/',opt$out_prefix, '_performance_measure_fillout.err -n 2 -R "rusage[mem=8]" -W 0:59 "Rscript ',opt$script,'performance_measure_script.R -g ', ground_file_name,' -t ', test_file_name, ' -d ',opt$directory,' -s ',name_test,' -n ',name_ground, ' -c ', opt$directory,' -m TRUE -p TRUE -o fillout_', opt$out_prefix)
  
  if(!is.null(opt$called_directory)){
    bsub_command <- paste0(bsub_command, ' -c ', opt$called_directory, ' -u ', opt$out_prefix)
  }
  
  
  if(is.null(opt$bed_file)){
    bsub_command <- paste0(bsub_command,'"')
  } else {
    bsub_command <- paste0(bsub_command, ' -b ',opt$bed_file ,'"')
  }
  
  system(bsub_command)
} 
