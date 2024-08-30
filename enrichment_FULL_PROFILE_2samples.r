# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)
library(tidyr)


make.enrichment.2samples.prot <- function(file.name, output.file) {
  
  # Read input data
  df <- read.csv(file.name)
  # Perform enrichment analysis for each column and collect results
  df_list <- lapply(colnames(df), function(i) {
    protlist <- df[[i]]
    enrichment <- rbioapi::rba_string_enrichment(protlist, 9606, split_df = F)
    categories_to_match <- c("Process", "Function","Component","KEGG","RCTM","WikiPathways")
    subset_df <- enrichment[enrichment$category %in% categories_to_match, ][, c("term", "description","category","preferredNames","number_of_genes","fdr")]
    colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_",i)
    colnames(subset_df)[colnames(subset_df) == "preferredNames"] <- paste0("genes_",i)
    colnames(subset_df)[colnames(subset_df) == "fdr"] <- paste0("fdr_",i)
    subset_df
  })
  
  # Merge enrichment results
  final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description", "category"), all = TRUE), df_list)
  colnames(final_df)
  # Replace NA values with 0
  final_df[is.na(final_df)] <- 0
  final_df <- final_df %>%
    mutate(across(4, ~ sapply(., function(x) paste(x, collapse = ", ")))) %>%
    mutate(across(7, ~ sapply(., function(x) paste(x, collapse = ", "))))
  
  custom_functioncount_g2_count_g1 <- function(row) {
    count_C <- as.numeric(row[5])
    count_P <- as.numeric(row[8])
    
    if ((count_C + count_P) != 0) {
      return (((count_P - count_C) / (count_P + count_C)) / 0.5)
    } else {
      return (0)
    }
  }
  final_df$DAVs <- apply(final_df, 1, custom_functioncount_g2_count_g1)
  write.csv(final_df, file = output.file, row.names = FALSE)
}


make.enrichment.2samples.prot("Striatum.csv","Striatum_enrichment.csv")
