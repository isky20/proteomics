# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)
library(tidyr)
make.enrichment.2samples.prot <- function(file.name, total.sample,g1,g2,g1name, g2name,TESTNAME) {
  df<- read.csv(file.name) # insert only the data frame with all possible replication
  
  # df <- aggregate(df[], list(df$Gene), FUN = max)
  
  # dim(df)
  
  
  df_list <- list()
  
  for (i in colnames(df)[2:length(colnames(df))] ){
    protlist <- df$Gene[df[i] > 0]
    
    enrichment <- rbioapi::rba_string_enrichment(protlist, 9606, split_df = F)
    
    #clean dataframe
    categories_to_match <- c("Process", "Function","Component","KEGG","RCTM","WikiPathways")
    subset_df <- enrichment[enrichment$category %in% categories_to_match, ][, c("term", "description","category","number_of_genes")]
    colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_",i) 
    df_list[[i]] <- subset_df
  }
  
  final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description","category"), all = TRUE), df_list)
  final_df[is.na(final_df)] <- 0 
  
  numeric_cols <- sapply(final_df, is.numeric)
  wsum.col.of.subjs <- colSums(final_df[, numeric_cols])
  
  ##calculate DAV without normalized data and cut the "term_name", "source" from the dataframe
  
  psm <- final_df[,5:ncol(final_df)] 
  
  Data.psm <- as.data.frame(psm)
  
  # calculate the mean for each protein throgh the cases and control
  Data.psm$Means_pre <- rowMeans(Data.psm[, grepl(g1name, names(Data.psm))]) #--#
  Data.psm$Means_post <- rowMeans(Data.psm[, grepl(g2name, names(Data.psm))]) #--#
  
  # build function to collect DAV for each comparison  #--#
  
  custom_functioncount_g2_count_g1<- function(row) {
    if ((row['Means_pre'] + row['Means_post']) != 0) {
      return (((row['Means_pre'] - row['Means_post']) / (row['Means_pre'] + row['Means_post'])) / 0.5)
    } else {
      return (0)
    }
  }
  
  # Apply the custom function to Data.psm.norm and create a new column holding the DAV values in each protein  #--#
  Data.psm$DAvs_count_pre_count_post <- apply(Data.psm, 1, custom_functioncount_g2_count_g1)
  
  
  #calculate p-value LDA
  psm.norm.t <- t(psm[,1:total.sample]) #--#
  
  labels <- c(rep(g1name,g1),rep(g2name,g2)) #--#
  
  final_psm_norm_t <- data.frame(psm.norm.t, labels)
  
  #MANOVA test
  dependent_vars <- cbind(as.matrix(final_psm_norm_t[,1:(ncol(final_psm_norm_t)-1)]))
  
  #independent variables (conditions)
  independent_var <- as.factor(final_psm_norm_t$labels)
  
  #MANOVA 
  #for manova model insert dependent variables before tilde and independent ones after. After comma indicate the dataframe cointaning both labels and quantitative data
  manova_model <- manova(dependent_vars ~ independent_var, data = final_psm_norm_t)
  #this command create a list of lists. Each one cointains the gene name and 2 values. the first one correspond to the Pvalue obtained from JMP.
  summary <-summary.aov(manova_model)
  
  output_LDA<- data.frame()
  for(i in 1:length(summary)){
    output_LDA<- rbind(output_LDA,summary[[i]][["Pr(>F)"]][1])
  }
  
  output_LDA$description <- final_df$description
  output_LDA$term <- final_df$term
  
  
  colnames(output_LDA)[1] <- "Pvalue"
  
  selected_sig <- output_LDA[output_LDA$Pvalue < 0.05,]
  
  # marge between two dataframe 
  Data.psm$description <- final_df$description
  Data.psm$term <- final_df$term
  
  merged_df <- merge(Data.psm , selected_sig, by = c("description" ,"term"), all.x = TRUE)
  
  
  the_data <- merge(merged_df, final_df[,c("description" ,"category" ,"term")], by = c("description","term"))
  
  no.na <- na.omit(the_data)
  df <- subset(no.na, !duplicated(no.na))
  
  new_column_order <- c(names(df)[1], "description", names(df)[-c(1, which(names(df) == "description"))])
  
  # Reorder the columns in the dataframe
  df <- df[, new_column_order]
  
  
  write.csv(df, paste(g1name,g2name,TESTNAME,"sig_enrichmentFC.csv", sep = "_"), row.names = FALSE)

}
make.enrichment.2samples.prot("Balb_norm_AP.csv",8,4,4,"pre","post","AP")
make.enrichment.2samples.prot("Balb_norm_CAP.csv",10,5,5,"pre","post","CAP")
