library(org.Hs.eg.db)
library(dplyr)

# Step 1: Load the data
count_matrix <- read.table("/Users/diyahafiz/Documents/Systems_Biology/Final_Project/GSE153873_summary_count.star.txt",
                           header = TRUE,
                           row.names = 1,
                           sep = "\t")

gene_info <- read.table("/Users/diyahafiz/Documents/Systems_Biology/Final_Project/mart_export.txt",
                        header = TRUE,
                        sep = "\t")

# Step 2: Clean up gene_info
names(gene_info) <- c("ENSEMBL", "Length")
gene_info <- gene_info[, c("ENSEMBL", "Length")]  # Keep only needed columns

# Step 3: Map gene symbols to Ensembl IDs
gene_symbols <- rownames(count_matrix)
symbol_to_ensembl <- mapIds(
  org.Hs.eg.db,
  keys = gene_symbols,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Step 4: Create clean mapping dataframe
mapping_df <- data.frame(
  SYMBOL = names(symbol_to_ensembl),
  ENSEMBL = symbol_to_ensembl,
  stringsAsFactors = FALSE
)

# Step 5: Merge data
# Add symbols to count matrix
count_matrix$SYMBOL <- rownames(count_matrix)

# Merge counts with mapping
merged_counts <- merge(mapping_df, count_matrix, by = "SYMBOL")

# Merge with gene lengths
final_data <- merge(merged_counts, gene_info, by = "ENSEMBL")
write.table(final_data,"Merged_IDS_starcount.csv")

# Step 6: Calculate TPM
# Get count columns
count_cols <- setdiff(colnames(final_data), c("SYMBOL", "ENSEMBL", "Length"))

# Calculate RPK (reads per kilobase)
rpk_matrix <- as.matrix(final_data[, count_cols]) / (final_data$Length / 1000)

# Calculate scaling factors
scaling_factors <- colSums(rpk_matrix) / 1e6

# Calculate TPM
tpm_matrix <- sweep(rpk_matrix, 2, scaling_factors, FUN = "/")
rownames(tpm_matrix) <- final_data$SYMBOL

# Step 7: Verify results
print("Data dimensions:")
print(paste("Original count matrix:", nrow(count_matrix), "genes"))
print(paste("After mapping and cleaning:", nrow(tpm_matrix), "genes"))
print(paste("Number of unique genes:", length(unique(rownames(tpm_matrix)))))
print(paste("Number of duplicated genes:", sum(duplicated(rownames(tpm_matrix)))))

# Step 8: Save results
# Save TPM matrix
write.table(tpm_matrix,
            "tpm_matrix_unique.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Create and save log-transformed version
log10_matrix <- log10(tpm_matrix + 1)
write.table(log10_matrix,
            "log10_tpm_matrix_unique.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Keep only the first occurrence of each SYMBOL
unique_data <- final_data[!duplicated(final_data$SYMBOL), ]

# Print dimensions to verify
print("Original dimensions:")
print(dim(final_data))
print("Dimensions after keeping unique symbols:")
print(dim(unique_data))
# Get count columns (exclude SYMBOL, ENSEMBL, Length)
count_cols <- setdiff(colnames(unique_data), c("SYMBOL", "ENSEMBL", "Length")
                      
# Step 1: Load the data from the text file
# Replace "your_file.txt" with the path to your file
data <- read.table("C:/Users/Gunja Gupta/Downloads/final_log10_tpm_matrix.txt", header = TRUE, row.names = 1)
                      
# Step 2: Convert row names into a column
data$RowNames <- rownames(data)
                      
# Step 3: (Optional) Reorder columns to place "RowNames" first
data <- data[, c("RowNames", setdiff(names(data), "RowNames"))]
                      
# Step 4: Save the modified data frame to a new text file
write.table(data, "modified_file.txt", sep = "\t", row.names = FALSE, quote = FALSE)
                      
# View the updated data frame
print(data)
                      
# Step 1: Read the text file
data <- read.table("modified_file.txt", header = TRUE)
                      
# Step 2: Change the name of the first column
colnames(data)[1] <- "Gene"  # Replace "NewColumnName" with your desired column name
                      
# Step 3: Save the modified data frame back to a text file (optional)
write.table(data, "modified_file.txt", sep = "\t", quote = FALSE)
                      
# View the updated data frame# View the updated data frameTRUE
print(data)