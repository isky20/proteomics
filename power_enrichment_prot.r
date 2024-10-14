# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(openxlsx)

process_enrichment <- function(obj, spid, int, output_file) {
  # Load the CSV file into a dataframe
  df <- read.csv(obj)
  
  # Aggregate the data using the "Gene" column, taking the max of each group
  df <- aggregate(df[], list(df$Gene), FUN = max)
  
  # Drop the old "Gene" column and rename the first column back to "Gene"
  df["Gene"] <- c()
  colnames(df)[1] <- "Gene"
  
  # Initialize an empty list to store the enrichment results
  df_list <- list()
  
  # Loop through each column name except the first one ("Gene")
  for (i in colnames(df)[2:length(colnames(df))]) {
    # Find genes with values greater than 0
    protlist <- df$Gene[df[i] > 0]
    
    # Perform enrichment analysis using rbioapi
    enrichment <- rbioapi::rba_string_enrichment(protlist, spid, split_df = F)
    
    # Filter relevant categories and clean up the dataframe
    categories_to_match <- c("Process", "Function", "Component", "KEGG", "RCTM", "WikiPathways")
    subset_df <- enrichment[enrichment$category %in% categories_to_match, ][, c("term", "description", "category", "preferredNames", "number_of_genes")]
    colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_", i)
    colnames(subset_df)[colnames(subset_df) == "preferredNames"] <- paste0("genes_", i)
    
    # Add the cleaned dataframe to the list
    df_list[[i]] <- subset_df
  }
  
  # Merge all the dataframes from df_list based on term, description, and category
  final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description", "category"), all = TRUE), df_list)
  final_df[is.na(final_df)] <- 0  # Replace NA values with 0
  
  ## Process the genes and clean up the dataframe
  gene_columns <- grep("^(term|genes_)", names(final_df), value = TRUE)
  gene_df <- final_df[, gene_columns]
  
  # Function to find unique gene lists across multiple columns
  find_unique <- function(lists) {
    combined_list <- unlist(lists)
    combined_list <- combined_list[combined_list != 0]  # Remove 0
    unique_values <- paste(unique(combined_list), collapse = ",")
    return(unique_values)
  }
  
  # Apply the function to create a unique gene list for each row
  gene_df$geneid <- apply(gene_df[, grepl("genes_", names(gene_df))], 1, find_unique)
  
  # Keep only the term and geneid columns
  gene_df <- gene_df[, c(1, ncol(gene_df))]
  
  ## Process the counts
  counts_df <- final_df[, !names(final_df) %in% grep("genes_", names(final_df), value = TRUE)]
  numeric_cols <- sapply(counts_df, is.numeric)
  wsum.col.of.subjs <- colSums(counts_df[, numeric_cols])
  
  ## Calculate means for int variables
  psm <- counts_df[, 4:ncol(counts_df)]
  Data.psm <- as.data.frame(psm)
  
  mean_list <- list()
  for (i in c(int)) {
    mean <- rowMeans(Data.psm[, grepl(i, names(Data.psm))])
    mean_list[[i]] <- mean
  }
  
  mean_df <- as.data.frame(mean_list)
  Data.psm <- cbind(Data.psm, mean_df)
  
  ## Calculate DAvs_X_Y for each combination of the int variables
  combinations <- t(combn(colnames(mean_df), 2))
  sub_dataframes <- list()
  
  DAvs_X_Y <- function(row) {
    if ((row[1] + row[2]) != 0) {
      return (((row[1] - row[2]) / (row[1] + row[2])) / 0.5)
    } else {
      return (0)
    }
  }
  
  for (i in 1:nrow(combinations)) {
    col1 <- combinations[i, 1]
    col2 <- combinations[i, 2]
    h <- paste0("DAvs_", col1, "_", col2)
    
    sub_df <- mean_df[, c(col1, col2)]
    Da <- apply(sub_df, 1, DAvs_X_Y)
    sub_dataframes[[h]] <- Da
  }
  
  sub_dataframes_df <- as.data.frame(sub_dataframes)
  Data.psm <- cbind(Data.psm, sub_dataframes_df)
  
  ## Perform MANOVA
  psm.norm.t <- t(psm[, 1:ncol(psm)])
  
  counts <- c()
  for (subword in int) {
    count <- sum(grepl(subword, colnames(df)))
    counts <- c(counts, count)
  }
  
  labels <- rep(int, counts)
  final_psm_norm_t <- data.frame(psm.norm.t, labels)
  
  dependent_vars <- as.matrix(final_psm_norm_t[, 1:(ncol(final_psm_norm_t) - 1)])
  independent_var <- as.factor(final_psm_norm_t$labels)
  
  manova_model <- manova(dependent_vars ~ independent_var, data = final_psm_norm_t)
  summary <- summary.aov(manova_model)
  
  output_LDA <- data.frame()
  for (i in 1:length(summary)) {
    output_LDA <- rbind(output_LDA, summary[[i]][["Pr(>F)"]][1])
  }
  
  output_LDA$description <- final_df$description
  output_LDA$term <- final_df$term
  colnames(output_LDA)[1] <- "Pvalue"
  
  selected_sig <- output_LDA[output_LDA$Pvalue < 0.05,]
  
  # Merge selected significant terms with Data.psm
  Data.psm$description <- final_df$description
  Data.psm$term <- final_df$term
  merged_df <- merge(Data.psm, selected_sig, by = c("description", "term"), all.x = TRUE)
  
  # Final processing and merging with gene_df
  the_data <- merge(merged_df, final_df[, c("description", "category", "term")], by = c("description", "term"))
  no.na <- na.omit(the_data)
  df <- subset(no.na, !duplicated(no.na))
  
  new_column_order <- c("category", names(df)[-which(names(df) == "category")])
  thedf <- merge(df[, new_column_order], gene_df, by = "term")
  
  # Save the result as an Excel file
  write.xlsx(thedf, output_file, overwrite = TRUE)
  return(thedf)
  
}

obj <- "test3.csv"
spid <- 9606
int <- c("Control","ALlambda", "ATTR")
output_file <- "enrichment_analysis_results_test3.xlsx"

result <- process_enrichment(obj, spid, int, output_file)
