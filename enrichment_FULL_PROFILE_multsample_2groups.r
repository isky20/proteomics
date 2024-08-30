# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)
library(tidyr)


make.enrichment.2samples.prot <- function(file.name, total.sample, g1, g2, g1name, g2name, TESTNAME) {
  
  # Read input data
  df <- read.csv(file.name)
  
  # Perform string enrichment for each column
  df_list <- lapply(colnames(df)[-1], function(i) {
    protlist <- df$Gene[df[i] > 0]
    enrichment <- rbioapi::rba_string_enrichment(protlist, 9606, split_df = FALSE)
    categories_to_match <- c("Process", "Function", "Component", "KEGG", "RCTM", "WikiPathways")
    subset_df <- enrichment[enrichment$category %in% categories_to_match, c("term", "description", "category", "number_of_genes")]
    colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_", i)
    subset_df
  })
  
  # Merge enrichment results
  final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description", "category"), all = TRUE), df_list)
  final_df[is.na(final_df)] <- 0
  
  # Calculate DAV
  Data.psm <- final_df[, 5:ncol(final_df)]
  Data.psm$Means_pre <- rowMeans(Data.psm[, grepl(g1name, names(Data.psm))])
  Data.psm$Means_post <- rowMeans(Data.psm[, grepl(g2name, names(Data.psm))])
  
  custom_functioncount_g2_count_g1 <- function(row) {
    if ((row['Means_pre'] + row['Means_post']) != 0) {
      return (((row['Means_pre'] - row['Means_post']) / (row['Means_pre'] + row['Means_post'])) / 0.5)
    } else {
      return (0)
    }
  }
  
  Data.psm$DAVs_count_pre_count_post <- apply(Data.psm, 1, custom_functioncount_g2_count_g1)
  
  # Perform MANOVA
  psm.norm.t <- t(Data.psm[, 1:total.sample])
  labels <- c(rep(g1name, g1), rep(g2name, g2))
  final_psm_norm_t <- data.frame(psm.norm.t, labels)
  
  dependent_vars <- as.matrix(final_psm_norm_t[, 1:(ncol(final_psm_norm_t) - 1)])
  independent_var <- as.factor(final_psm_norm_t$labels)
  manova_model <- manova(dependent_vars ~ independent_var, data = final_psm_norm_t)
  summary <- summary.aov(manova_model)
  
  # Extract P-values
  output_LDA <- data.frame()
  for (i in 1:length(summary)) {
    output_LDA <- rbind(output_LDA, summary[[i]][["Pr(>F)"]][1])
  }
  output_LDA$description <- final_df$description
  output_LDA$term <- final_df$term
  colnames(output_LDA)[1] <- "Pvalue"
  
  # Select significant results based on P-value threshold
  selected_sig <- output_LDA[output_LDA$Pvalue < 0.05, ]
  
  # Merge significant results with Data.psm
  Data.psm$description <- final_df$description
  Data.psm$term <- final_df$term
  merged_df <- merge(Data.psm, selected_sig, by = c("description", "term"), all.x = TRUE)
  the_data <- merge(merged_df, final_df[, c("description", "category", "term")], by = c("description", "term"))
  no.na <- na.omit(the_data)
  df <- subset(no.na, !duplicated(no.na))
  
  # Reorder columns
  new_column_order <- c(names(df)[1], "description", names(df)[-c(1, which(names(df) == "description"))])
  df <- df[, new_column_order]
  
  # Write output to CSV
  write.csv(df, paste(g1name, g2name, TESTNAME, "sig_enrichmentFC.csv", sep = "_"), row.names = FALSE)
}

