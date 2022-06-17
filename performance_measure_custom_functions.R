#########################################################################
########################### Functions ###################################
#negate
'%nin%' = Negate('%in%')


mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


# Function to break up dataframe based on given variable
### This function will always break down the variable based on Variable_Type then calculater statistics
#### For general information across a variable name  look at 'all' Variant type in returned table/graphs

parse_dataframe_on_var <- function(ground_df,test_df,variable_id,type_of_analysis){
  if(variable_id == 'substitutions'){
    ground_df <- ground_df[ground_df$Variant_Type == 'SNV',]
    test_df <- test_df[test_df$Variant_Type == 'SNV',]
  }
  
  levels_from_variable_id <- unique(c(ground_df[,variable_id],test_df[,variable_id]))
  levels_from_variable_id <- levels_from_variable_id[!is.na(levels_from_variable_id)]
  names(levels_from_variable_id) <- levels_from_variable_id
  
  variable_id_stats <- lapply(levels_from_variable_id, function(id) {
    targeted_ground <- ground_df %>% filter(get(variable_id) == id) 
    targeted_test <- test_df %>% filter(get(variable_id) == id) 
    if(variable_id == 'substitutions'){
        output <-f1_stats(targeted_ground,targeted_test,type_of_analysis)
        output[variable_id] <- id
        output['type'] <- 'SNV'
        output <- output %>% filter(tag_type == 'restrictive')
      
    } else {
      output <- calc_stats_by_variant_type(targeted_ground,targeted_test,type_of_analysis) 
      output[,variable_id] <- id
    } 
    
    return(output)
    
  })
  
  # Formatting
  variable_id_stats <- as.data.frame(do.call(rbind,variable_id_stats))
  return(variable_id_stats)
}

calc_stats_by_variant_type <- function(ground,test,type_of_analysis) {
  variant_types <- unique(c('all','indel',ground$Variant_Type,test$Variant_Type))
  output1 <- lapply(variant_types, function(type) {
    
    if(type == 'all'){
      targ_g <- ground
      targ_t <- test
      
    } else if(type == 'indel'){
      
      indels <- c('INS','DEL')
      targ_g <- ground %>% filter(Variant_Type %in% indels)
      targ_t <- test %>% filter(Variant_Type %in% indels)
      
      
    }else {
      
      targ_g <-ground %>% filter(Variant_Type == type)
      targ_t <- test %>% filter(Variant_Type == type)
      
    }
    stats_return <- f1_stats(targ_g,targ_t,type_of_analysis)
    stats_return$type <- type
    return(stats_return)
    
    
  })
  output1 <- as.data.frame(do.call(rbind,output1))
  cols.num <- c(seq(2,9,1),seq(11,13,1))
  output1[cols.num] <- lapply(output1[cols.num], as.numeric)
  return(output1)
}

# Calculate Statistics 
## Function will calculate recall, precision and f1 with error bars for dataframe
## Df must which contain var_tag variable

f1_stats <- function(ground_set,test_set,type_of_analysis){
  n_samples <- length(unique(c(ground_set$Tumor_Sample_Barcode,test_set$Tumor_Sample_Barcode)))
  if(any(c(colnames(ground_set),colnames(test_set)) == 'detectable')){
    run_types <- c('restrictive')
    
  } else {
    run_types <- c('restrictive','permissive')
    
  }
  stats <- do.call(rbind,lapply(run_types, function(tag_type){
    
    
    # If fillouts has been run, (detectable present), ensure opposing set is detectable at that location 
    if(any(c(colnames(ground_set),colnames(test_set)) == 'detectable')){
      test_set_must <- test_set[!test_set$evidence & test_set$detectable,]
      ground_set_must <- ground_set[!ground_set$evidence & ground_set$detectable,]
        tps <- length(ground_set$var_tag[ground_set$evidence & ground_set$var_tag %in% test_set$var_tag[test_set$evidence]])
        
        ## Remove variants from count if test_set$var_tag is NOT detectable (these variants cannot be used in analysis as FNs)
        fns <- ground_set$var_tag[ground_set$evidence & ground_set$var_tag %nin% test_set$var_tag[test_set$evidence]]
        test_set_no_ev_not_detect <- length(fns[which(fns %nin% test_set_must$var_tag)])
        fns <- length(fns[which(fns %in% test_set_must$var_tag)])
        
        ## Remove variants from count if ground_set$var_tag is NOT detectable (these variants cannot be used in analysis as FPs)
        fps <- test_set$var_tag[test_set$evidence & test_set$var_tag %in% ground_set$var_tag[!ground_set$evidence]]
        ground_set_no_ev_not_detect <- length(fps[which(fps %nin% ground_set_must$var_tag)])
        
        fps <- length(fps[which(fps %in% ground_set_must$var_tag)])
        vars_with_no_evidence_in_either_test_or_ground <- length(ground_set$var_tag[!ground_set$evidence & ground_set$var_tag %in%  test_set$var_tag[!test_set$evidence]])

      total_var_count <- tps + fps + fns + ground_set_no_ev_not_detect + test_set_no_ev_not_detect + vars_with_no_evidence_in_either_test_or_ground
      
    } else{
      if(tag_type == 'restrictive'){
        test_set_tags <- test_set$var_tag
        ground_set_tags <- ground_set$var_tag
      } else {
        test_set_tags <- test_set$TAG
        ground_set_tags <- ground_set$TAG
      }
      
      
      tps <- length(test_set_tags[test_set_tags %in% ground_set_tags])
      fps <- length(test_set_tags[test_set_tags %nin% ground_set_tags])
      fns <- length(ground_set_tags[ground_set_tags %nin% test_set_tags])
      test_set_no_ev_not_detect <- NA
      ground_set_no_ev_not_detect <- NA
      vars_with_no_evidence_in_either_test_or_ground <- NA
      total_var_count <- tps + fps + fns 
      
    }
    precision <- tps / (tps + fps)
    recall <- tps / (tps + fns)
    f1 <- (2 * precision * recall)/(precision + recall)
    
    bin_conf_recall <- binom.confint(tps,(tps + fns), conf.level = 0.95, method = 'wilson')
    bin_conf_precision <- binom.confint(tps,(tps + fps), conf.level = 0.95, method = 'wilson')
    # If cohort level bootstrap, if not return NA

      
    f1_confidence <- c((2 * bin_conf_recall$lower * bin_conf_precision$lower)/(bin_conf_precision$lower + bin_conf_recall$lower),(2 * bin_conf_recall$upper * bin_conf_precision$upper)/(bin_conf_precision$upper + bin_conf_recall$upper))
      
    
    repo <- c(tag_type,total_var_count,n_samples,tps,fps,fns,ground_set_no_ev_not_detect,test_set_no_ev_not_detect,vars_with_no_evidence_in_either_test_or_ground)
    stats_p <- c(repo,'precision',precision,bin_conf_precision$lower,bin_conf_precision$upper)
    stats_r <- c(repo,'recall', recall,bin_conf_recall$lower,bin_conf_recall$upper)
    stats_f <- c(repo,'f1_Score',f1,f1_confidence[1],f1_confidence[2])
    
    stats <- as.data.frame(rbind(stats_r,stats_p,stats_f))
    colnames(stats) <- c('tag_type','total_var_count','n_samples','tps','fps','fns','ground_set_no_ev_not_detect','test_set_no_ev_not_detect','vars_with_no_evidence_in_either_test_or_ground','statistic_name','value','lower','upper')
    col.nums <- c(seq(2,9,1),seq(11,13,1))
    
    stats[col.nums] <- lapply(stats[col.nums], as.numeric)
    
    stats
  }))
  
  return(stats)
}
statistics_graphs <- function(dataframe,variable_id,graph_type,dir,out,opt){

  if(variable_id == 'type'){
    dataframe <- dataframe %>% filter(type %nin% c('INS','DEL'))
  }
  
  general_theme <- theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=1,size = 12),legend.position = "none",legend.background=element_blank(),legend.title=element_blank())
  
  if(!opt$fillout_to_pr){
    col_scale <-  if(variable_id != 'substitutions')
    {
      scale_fill_jama()
    } else{
      scale_fill_manual(
        values = c("C>T" = "red","C>G" = "black","C>A" = "skyblue", "T>A" = "gray","T>C" = "lightgreen","T>G" ="pink"),
        labels = c("C>T", "C>G", "C>A","T>A","T>C","T>G")
      )
    }
    
    for (stat_name in unique(dataframe$statistic_name)){
      assign(paste0('base_',stat_name), (ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=value,fill = get(variable_id))) + col_scale+ 
                                                          scale_x_discrete(labels = unique(dataframe[variable_id]))+
                                                      labs(x=' ', y = stat_name) + general_theme + ylim(0,1)))
      if( graph_type == 'bar') {
        assign(paste0('base_',stat_name), (get(paste0('base_',stat_name)) +geom_bar(stat='identity',position=position_dodge(),width=0.75)  +
                                             geom_errorbar(aes(ymin=upper, ymax = lower),position = position_dodge(0.75), width =0.65)))
      }else {
        assign(paste0('base_',stat_name), (get(paste0('base_',stat_name)) +geom_boxplot()))
      }
  
      }
    
  
    base_tps <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=tps,fill = get(variable_id))) + col_scale +
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'True Positive Count') + general_theme + ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'tps'])+1)
    base_fps <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=fps,fill = get(variable_id))) + col_scale +
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'False Positive Count') + general_theme+ ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'fps'])+1)
    
    base_fns <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=fns,fill = get(variable_id))) + col_scale + 
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'False Negative Count') + general_theme + ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'fns'])+1)
  } else{
    
    col_scale <- scale_fill_jama()
    
    
    for (stat_name in unique(dataframe$statistic_name)){
      assign(paste0('base_',stat_name), (ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=value,fill = Genotyped)) + col_scale+ 
                                           scale_x_discrete(labels = unique(dataframe[variable_id]))+
                                           labs(x=' ', y = stat_name) + general_theme + ylim(0,1)))
      if( graph_type == 'bar') {
        assign(paste0('base_',stat_name), (get(paste0('base_',stat_name)) +geom_bar(stat='identity',position=position_dodge(),width=0.75)  +
                                             geom_errorbar(aes(ymin=upper, ymax = lower, group = Genotyped),position = position_dodge(0.75), width =0.65)))
      }else {
        assign(paste0('base_',stat_name), (get(paste0('base_',stat_name)) +geom_boxplot()))
      }
      
    }
    
    
    base_tps <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=tps,fill = Genotyped)) + col_scale +
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'True Positive Count') + general_theme + ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'tps'])+1)
    base_fps <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=fps,fill = Genotyped)) + col_scale +
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'False Positive Count') + general_theme+ ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'fps'])+1)
    
    base_fns <- ggplot(dataframe[dataframe$statistic_name == stat_name,],aes(x=get(variable_id),y=fns,fill = Genotyped)) + col_scale + 
      scale_x_discrete(labels = unique(dataframe[variable_id]))+
      labs(x=' ', y = 'False Negative Count') + general_theme + ylim(0,max(dataframe[dataframe$statistic_name == stat_name,'fns'])+1)
  }
  
  
  
  
  
  if (graph_type == 'bar') {
    base_tps <- base_tps + geom_bar(stat='identity',position=position_dodge(),width=0.75)  
    base_fps  <- base_fps + geom_bar(stat='identity',position=position_dodge(),width=0.75)  
    base_fns <- base_fns + geom_bar(stat='identity',position=position_dodge(),width=0.75)  
  }else {
    base_tps <- base_tps + geom_boxplot()
    
    base_fns <- base_fns + geom_boxplot()
    base_fps <- base_fps + geom_boxplot()
  }
  if(opt$multiqc){
    png(paste0(opt$mq_dir,variable_id,'_',graph_type,'_mqc.png'),type = 'cairo', bg = "white",width=1000)
    
    print(annotate_figure(ggarrange(base_recall,base_precision,base_f1_Score,base_tps,base_fps,base_fns, ncol=3,nrow=2,common.legend=TRUE),bottom = variable_id, top = paste0('Statistics for ',variable_id)))
    dev.off()
  }
  pdf(paste0(dir,'images/',out,'_',variable_id,'_',graph_type,'.pdf'),paper = 'a4r')  


    print(annotate_figure(ggarrange(base_recall,base_precision,base_f1_Score,base_tps,base_fps,base_fns, ncol=3,nrow=2,common.legend=TRUE),bottom = variable_id, top = paste0('Statistics for ',variable_id)))
  dev.off()
  
}

restruct_for_multiqc <- function(df,variable,level,directory){

  if(any(colnames(df) == 'Genotyped')){
    df[,variable] <- paste0(df[,variable],'/',df$Genotyped)
  }
  
  statistics_to_parse <- c('recall','precision','f1_Score','total_var_count','tps','fps','fns','ground_set_no_ev_not_detect','test_set_no_ev_not_detect','vars_with_no_evidence_in_either_test_or_ground')
  names(statistics_to_parse) <- statistics_to_parse
  
  if(level != 'sample'){
    tmp <- pivot_wider(df[,c(variable,'statistic_name','value')], names_from = 'statistic_name', values_from = "value")
    tmp2 <- unique(df[,c(variable,'total_var_count','n_samples','tps','fps','fns','ground_set_no_ev_not_detect','test_set_no_ev_not_detect','vars_with_no_evidence_in_either_test_or_ground')])
    tmp3 <- merge(tmp,tmp2)

    tmp3$ID <- tmp3[,variable]
    if(any(colnames(df) == 'Genotyped')){
      tmp3 <- tmp3 %>% separate(get(variable), c(variable, "Genotyped"), "/")

      tmp3 <- tmp3[,c('ID',variable,'Genotyped','n_samples',statistics_to_parse)]
      
    }else {
      tmp3 <- tmp3[,c('ID',variable,'n_samples',statistics_to_parse)]
    }
    
    suppressWarnings(cat("#plot_type: 'table' \n ",file=paste0(directory,variable,'_',level,'_mqc.tsv')))
    suppressWarnings(write.table(tmp3,paste0(directory,variable,'_',level,'_mqc.tsv'),sep= '\t',append=TRUE,row.names = FALSE))
  } else {
    df$ID <- paste(df$Tumor_Sample_Barcode,df[,variable], sep = '/')
    
    tmp <- pivot_wider(df[,c('ID','statistic_name','value')], names_from = 'statistic_name', values_from = "value")
    tmp2 <- unique(df[,c('ID','total_var_count','tps','fps','fns','ground_set_no_ev_not_detect','test_set_no_ev_not_detect','vars_with_no_evidence_in_either_test_or_ground')])
    tmp <- merge(tmp,tmp2)
    tmp3 <- tmp
   #  variables <- unique(tmp[,variable])
   #  names(variables) <- unique(tmp[,variable])
   # 
   #  tmp3 <-  lapply(variables,function(variable_id){
   #    quantile_reports <-  lapply(statistics_to_parse,function(stat) {
   #      
   #        return(quantile(tmp[tmp[,variable] == variable_id,stat], na.rm = TRUE))
   #      
   #    })
   #    quantile_reports <- as.data.frame(do.call(cbind, quantile_reports))
   #    quantile_reports[,'Quantile'] <- row.names(quantile_reports)
   #    quantile_reports[,variable] <- variable_id
   #    quantile_reports[,'ID'] <- paste0(variable_id,'_',quantile_reports$Quantile)
   #    return(quantile_reports)
   # }) 

    if(any(colnames(df) == 'Genotyped')){
      tmp3 <- tmp3 %>% separate(ID, c('Tumor_Sample_Barcode',variable, "Genotyped"), "/",remove = FALSE)
      tmp3 <- tmp3[,c('ID','Tumor_Sample_Barcode',variable,'Genotyped',statistics_to_parse)]
    }else {
      tmp3 <- tmp3 %>% separate(ID, c('Tumor_Sample_Barcode',variable), "/",remove = FALSE)
      tmp3 <- tmp3[,c('ID','Tumor_Sample_Barcode',variable,statistics_to_parse)]
    }
    if(nrow(tmp3) <= 500){
      suppressWarnings(cat("#plot_type: 'table' \n ",file=paste0(directory,variable,'_',level,'_mqc.tsv')))
      suppressWarnings(write.table(tmp3,paste0(directory,variable,'_',level,'_mqc.tsv'),sep= '\t',append=TRUE,row.names = FALSE))
    }
  }
  

}
