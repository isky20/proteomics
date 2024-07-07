library(ggplot2)
library(dplyr)

generate_plots <- function(df, file_name) {
  # Step 1: Select the first three columns
  first_three_columns <- df[, 1:3]
  
  # Step 2: Identify columns starting with 'DAV'
  dav_columns <- grep("^DAv", names(df), value = TRUE)
  
  # Step 3: Create sub-dataframes
  sub_dfs <- list()
  for (dav_col in dav_columns) {
    sub_df <- cbind(first_three_columns, df[, dav_col, drop = FALSE])
    filtered_data <- subset(sub_df, abs(sub_df[[dav_col]]) > 0.4)
    filtered_data <- filtered_data %>%
      mutate(MainCategory = case_when(
        category %in% c("WikiPathways", "RCTM", "KEGG") ~ "Pathways",
        category == "Function" ~ "Function",
        category == "Component" ~ "Component",
        category == "Process" ~ "Process",
        TRUE ~ "other"  # Default category if no match
      ))
    sub_dfs[[dav_col]] <- filtered_data
  }
  
  separated_dataframes <- list()
  column_to_split <- "MainCategory"
  
  # Loop through each dataframe
  for (name in names(sub_dfs)) {
    df <- sub_dfs[[name]]
    
    # Split dataframe based on column
    grouped <- split(df, df[[column_to_split]])
    
    # Store grouped dataframes in the separated_dataframes list
    separated_dataframes[[name]] <- grouped
  }
  
  for (name in names(separated_dataframes)) {
    df_list <- separated_dataframes[[name]]
    
    # Loop through each category in df_list
    for (category_name in names(df_list)) {
      df <- df_list[[category_name]]
      
      # Create a plot (example: bar plot of counts by MainCategory)
      plot_title <- paste("Plot for", name, "-", category_name)
      p <- ggplot(df, aes(x = reorder(description, df[[4]]), y = df[[4]], fill = category)) +
        geom_col(stat = "identity", width = 0.1) +  # Use geom_col() to create column chart
        scale_fill_brewer(palette = "Accent") +  # Use a Brewer color palette for fill colors
        labs(title = plot_title, x = "Description", y = name) +   # Labels for title, x-axis, and y-axis
        theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Adjusted angle and hjust
              axis.text.y = element_text(size = 10),  # Adjusted size for y-axis text
              panel.grid.major = element_line(color = "white"),
              panel.grid.minor = element_blank()) +
        coord_flip() +
        ylim(-2, 2)
      
      # Save the plot
      ggsave(filename = paste(name, "_", category_name, "_", file_name, ".svg", sep = ""), plot = p,
             width = 8, height = 6, units = "in", dpi = 600, bg = "white")
    }
  }
}

# Example usage:
df <- read.csv("sig_enrich_PSMs_Breast.csv")
generate_plots(df, "PSMs_breast")
