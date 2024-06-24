# load the library			
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)
library(tidyr)

rba_options()
df<- read.csv("Gliofil_norm.csv") # insert only the data frame with all possible replication
dim(df)

#df <- aggregate(df[], list(df$Gene), FUN = max)

head(df)


df_list <- list()

for (i in colnames(df)[2:length(colnames(df))] ){
  protlist <- df$Gene[df[i] > 0]
  
  enrichment <- rbioapi::rba_string_enrichment(protlist, 10090, split_df = F)
  
  #clean dataframe
  categories_to_match <- c("Process", "Function","Component","KEGG","RCTM","WikiPathways")
  subset_df <- enrichment[enrichment$category %in% categories_to_match, ][, c("term", "description","category","preferredNames","number_of_genes")]
  colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_",i)
  colnames(subset_df)[colnames(subset_df) == "preferredNames"] <- paste0("genes_",i) 
  
  df_list[[i]] <- subset_df
}

final_df <- Reduce(function(x, y) 
  merge(x, y, by = c("term", "description","category"),
        all = TRUE), df_list)

final_df[is.na(final_df)] <- 0 


##clean gene_df

#clean the data frame and make gene list and clean it 
gene_columns <- grep("^(term|genes_)", names(final_df), value = TRUE)
gene_df <- final_df[, gene_columns]
# Function to find common elements between multiple lists

find_unique <- function(lists) {
  combined_list <- unlist(lists)
  # Remove 0 from the combined list
  combined_list <- combined_list[combined_list != 0]
  unique_values <- paste(unique(combined_list), collapse = ",")
  return(unique_values)
}
# Apply the function to find common elements
gene_df$geneid <- apply(gene_df[, grepl("genes_", names(gene_df))], 1, find_unique)

gene_df <- gene_df[, c(1,ncol(gene_df))]


## count list 
counts_df  <- final_df[, !names(final_df) %in% grep("genes_", names(final_df), value = TRUE)]

## check if it is numeric and sum the bd reads  
numeric_cols <- sapply(counts_df , is.numeric)
wsum.col.of.subjs <- colSums(counts_df [, numeric_cols])

##calculate DAV without normalized data and cut the "term_name", "source" from the dataframe

psm <- counts_df[,4:ncol(counts_df)] 

Data.psm <- as.data.frame(psm)

# calculate the mean for each protein throgh the cases and control
Data.psm$count_g1 <- rowMeans(Data.psm[, grepl("count_C", names(Data.psm))]) #--#
Data.psm$count_g2 <- rowMeans(Data.psm[, grepl("count_NT", names(Data.psm))]) #--#
Data.psm$count_g3 <- rowMeans(Data.psm[, grepl("count_T", names(Data.psm))]) #--#

# build function to collect DAV for each comparison  #--#

custom_functioncount_g2_count_g1<- function(row) {
  if ((row['count_g1'] + row['count_g2']) != 0) {
    return (((row['count_g1'] - row['count_g2']) / (row['count_g1'] + row['count_g2'])) / 0.5)
  } else {
    return (0)
  }
}

custom_functioncount_g2_count_g3<- function(row) {
  if ((row['count_g3'] + row['count_g2']) != 0) {
    return (((row['count_g3'] - row['count_g2']) / (row['count_g3'] + row['count_g2'])) / 0.5)
  } else {
    return (0)
  }
}

custom_function_count_g1_count_g3<- function(row) {
  if ((row['count_g1'] + row['count_g3']) != 0) {
    return (((row['count_g1'] - row['count_g3']) / (row['count_g1'] + row['count_g3'])) / 0.5)
  } else {
    return (0)
  }
}


# Apply the custom function to Data.psm.norm and create a new column holding the DAV values in each protein  #--#
Data.psm$DAvs_count_g2_count_g1 <- apply(Data.psm, 1, custom_functioncount_g2_count_g1)
Data.psm$DAvs_count_g2_count_g3 <- apply(Data.psm, 1, custom_functioncount_g2_count_g3)
Data.psm$DAvs_count_g1_count_g3 <- apply(Data.psm, 1, custom_function_count_g1_count_g3)

#calculate p-value LDA
dim(psm)
psm.norm.t <- t(psm[,1:18]) #--#

labels <- c(rep("count_g1",6),rep("count_g2",6),rep("count_g3",6)) #--#

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

new_column_order <- c("category",names(df)[-c(which(names(df) == "category"))])

# Reorder the columns in the dataframe
df <- df[, new_column_order]


thedf <- merge(df, gene_df, by = "term")


write.csv(thedf, "sig_enrichment.csv", row.names = FALSE)
