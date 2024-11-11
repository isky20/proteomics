import sys
import pandas as pd

def filter_and_save(input_file, output_file, D, pvalue_threshold):
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Filter the DataFrame based on columns starting with 'DAv' and pvalue sorting
    df_filtered = df[
        (df.filter(like='DAv').abs().ge(D).any(axis=1)) &
        (df['Pvalue'] <= pvalue_threshold)
    ].sort_values(by='Pvalue', ascending=True)
    
    # Get unique categories
    categories = df_filtered['category'].unique()
    
    # Sort each category by Pvalue and take the top 20 rows
    dfs_by_category = {
        cat: df_filtered[df_filtered['category'] == cat].sort_values(by='Pvalue', ascending=True).head(20)
        for cat in categories
    }
    
    # Combine all category DataFrames into one DataFrame
    df_combined = pd.concat(dfs_by_category.values(), axis=0)
    
    # Reset the index
    df_combined.reset_index(drop=True, inplace=True)
    
    # Save the combined DataFrame to the output CSV file
    df_combined.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 5:
        print("Usage: python3 code.py <input_file> <output_file> <D> <pvalue_threshold>")
        sys.exit(1)
    
    # Get input and output file paths and D and Pvalue values from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    D = float(sys.argv[3])
    pvalue_threshold = float(sys.argv[4])
    
    # Call the function to filter and save the data
    filter_and_save(input_file, output_file, D, pvalue_threshold)
