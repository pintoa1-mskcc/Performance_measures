library(dplyr)

bam_dir <- '/Users/apinto1/juno/work/tempo/wes_repo/Results/v1.4.x/cohort_level/CCS_NFCXOZVP/bams/'
ccs_mapping <- read.table("../../fillout_testing/CCS_NFCXOZVP.cohort.txt",header = TRUE)

all_samples <- list.files(bam_dir)
bams <- c()
for (i in 1:length(all_samples)){
  list.files(paste0(bam_dir,all_samples[i]))
  bams <- c(bams,paste0(bam_dir,all_samples[i],'/',all_samples[i],'.bam'))
}

ccs_mapping$TUMOR_BAM <- sapply(ccs_mapping$TUMOR_ID, function(id) {bams[grep(id,bams)]})
ccs_mapping$NORMAL_BAM <- sapply(ccs_mapping$NORMAL_ID, function(id) {bams[grep(id,bams)]})
colnames(ccs_mapping) <- tolower(colnames(ccs_mapping))

write.table(ccs_mapping,'../develop/CCS_bam_mapping_for_fillouts.txt',sep = '\t',quote = FALSE)
