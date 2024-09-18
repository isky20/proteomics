# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)
library(tidyr)
library(openxlsx)

df<- read.csv("tot.csv") # insert only the data frame with all possible replication
head(df)


df_list <- list()

for (i in colnames(df)[2:length(colnames(df))] ){
  protlist <- df$Gene[df[i] > 0]
  enrichment <- rbioapi::rba_string_enrichment(protlist, 10090, split_df = F)
  categories_to_match <- c("Process", "Function","Component","KEGG","RCTM","WikiPathways")
  subset_df <- enrichment[enrichment$category %in% categories_to_match, ][, c("term", "description","category","preferredNames","number_of_genes","fdr")]
  colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_",i)
  colnames(subset_df)[colnames(subset_df) == "preferredNames"] <- paste0("genes_",i)
  colnames(subset_df)[colnames(subset_df) == "fdr"] <- paste0("fdr_",i)
  df_list[[i]] <- subset_df
}

final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description","category"), all = TRUE), df_list)
final_df[is.na(final_df)] <- 0



genes_y_counts <- grep("count.*WT", names(final_df), value = TRUE)
final_df$Means_WT <- rowMeans(final_df[, genes_y_counts])

genes_z_counts <- grep("count.*DEL", names(final_df), value = TRUE)
final_df$Means_DEL <- rowMeans(final_df[, genes_z_counts])


final_df <- final_df %>%
  mutate(across(starts_with("genes"), ~ sapply(., function(x) paste(x, collapse = ", "))))


genes_Y_columns <- grep("genes.*WT", names(final_df), value = TRUE)
genes_z_columns <- grep("genes.*DEL", names(final_df), value = TRUE)

# Create a unique string for each row
final_df$genesWT <- apply(final_df[genes_Y_columns], 1, function(row) paste(unique(row), collapse = ", "))
final_df$genesDEL <- apply(final_df[genes_z_columns], 1, function(row) paste(unique(row), collapse = ", "))




final_df <- final_df %>%
  select(-all_of(genes_Y_columns))%>%
  select(-all_of(genes_z_columns))

genes_fdr_columns <- grep("fdr", names(final_df), value = TRUE)

final_df <- final_df %>%
  select(-all_of(genes_fdr_columns))



custom_functionATTR_ALlambda<- function(row) {
  if ((row['Means_WT'] + row['Means_DEL']) != 0) {
    return (((row['Means_DEL'] - row['Means_WT']) / (row['Means_WT'] + row['Means_DEL'])) / 0.5)
  } else {
    return (0)
  }
}



final_df$DAvs_DEL_WT <- apply(final_df[, c("Means_WT", "Means_DEL")], 1, custom_functionATTR_ALlambda)

dim(final_df)

psm <- final_df[,4:15]
head(psm)
psm.norm.t <- t(psm) #--#

labels <- c(rep("WT",6),rep("DEL",6))

final_psm_norm_t <- data.frame(psm.norm.t, labels)

#MANOVA test
dependent_vars <- cbind(as.matrix(final_psm_norm_t[,1:(ncol(final_psm_norm_t))-1]))

#independent variables (conditions)
independent_var <- as.factor(final_psm_norm_t$labels)

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

merged_df <- merge(final_df , selected_sig, by = c("description" ,"term"), all.x = TRUE)





write.xlsx(merged_df, "sig_enrichment_tot.xlsx")
