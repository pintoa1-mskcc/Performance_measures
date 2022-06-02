library(argparse) 
library(dplyr)
opt = commandArgs(TRUE)

parser=ArgumentParser()
parser$add_argument("-b", "--bam_dir", type='character', default=NULL,
                    help="Directory containing BAMS ")
parser$add_argument("-m",'--mapping',type='character',default=NULL,help='Tumor-Normal mapping. Expected header format: TUMOR_ID /t NORMAL_ID')
parser$add_argument("-o",'--output_file',type="character",default = NULL)
opt=parser$parse_args()

bam_dir <- opt$bam_dir
ccs_mapping <- read.table(opt$mapping,header = TRUE)

all_samples <- list.files(bam_dir)
bams <- c()
for (i in 1:length(all_samples)){
  list.files(paste0(bam_dir,all_samples[i]))
  bams <- c(bams,paste0(bam_dir,all_samples[i],'/',all_samples[i],'.bam'))
}

ccs_mapping$TUMOR_BAM <- sapply(ccs_mapping$TUMOR_ID, function(id) {bams[grep(id,bams)]})
ccs_mapping$NORMAL_BAM <- sapply(ccs_mapping$NORMAL_ID, function(id) {bams[grep(id,bams)]})
colnames(ccs_mapping) <- tolower(colnames(ccs_mapping))

write.table(ccs_mapping,opt$output_file,sep = '\t',quote = FALSE)
