# ===================================================================
#   AI and Omics Research Internship - Module I / Assignment II
#   Differential Gene Expression (DGE) Analysis Script
#
#   Author: Sulaimon Araromi
#   Date:   August 31, 2025
# ===================================================================


# Task 1: Define the Gene Classification Function
# -------------------------------------------------------------------
# This function categorizes genes based on log2 fold change and adjusted p-value.

classify_gene <- function(logFC, padj, log2FC_threshold = 1, padj_threshold = 0.05) {
  # Handle cases where p-value or logFC is missing
  if (is.na(padj) || is.na(logFC)) {
    return("Not_Significant")
  }
  
  # First, check for statistical significance
  if (padj < padj_threshold) {
    # If significant, check the direction of expression change
    if (logFC > log2FC_threshold) {
      return("Upregulated")
    } else if (logFC < -log2FC_threshold) {
      return("Downregulated")
    } else {
      return("Not_Significant") # Significant p-value but logFC is below threshold
    }
  } else {
    return("Not_Significant") # Not statistically significant
  }
}


# Task 2: Define a Function to Process DGE Datasets
# -------------------------------------------------------------------
# This function reads a DGE file, applies the classification, and returns a summary.

process_dge_data <- function(file_path) {
  # Read the dataset
  cat("--> Reading data from:", file_path, "\n")
  data <- read.csv(file_path)
  
  # Handle missing padj values as per assignment instructions

  missing_padj_count <- sum(is.na(data$padj))
  if (missing_padj_count > 0) {
    cat("--> Found and replaced", missing_padj_count, "missing padj values with 1.\n")
    data$padj[is.na(data$padj)] <- 1
  }
  
  # Apply the classification function to every row of the dataset
  # mapply is used here to iterate the function over the logFC and padj columns
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  
  # Print a summary of the results for this file
  cat("--> Classification Summary:\n")
  print(table(data$status))
  
  return(data)
}


# Task 3: Process Both DGE Datasets Using the Function
# -------------------------------------------------------------------
# Ensure the .csv files are in your 'raw_data' sub-folder.

cat("====================================================\n")
cat("Processing DEGs_Data_1.csv\n")
cat("====================================================\n")

processed_data_1 <- process_dge_data("raw_data/DEGs_Data_1.csv")

cat("\n====================================================\n")
cat("Processing DEGs_Data_2.csv\n")
cat("====================================================\n")

processed_data_2 <- process_dge_data("raw_data/DEGs_Data_2.csv")


# Task 4: Save the Processed Dataframes to CSV files
# -------------------------------------------------------------------
# The processed data, including the new 'status' column, is saved.

write.csv(processed_data_1, "processed_DEGs_Data_1.csv", row.names = FALSE)
write.csv(processed_data_2, "processed_DEGs_Data_2.csv", row.names = FALSE)
cat("\n--> Processed datasets have been saved to 'processed_DEGs_Data_1.csv' and 'processed_DEGs_Data_2.csv'.\n")


# Task 5: Generate a Final Summary Report
# -------------------------------------------------------------------
# This creates a summary table of the analysis results.

summary_table <- data.frame(
  Dataset = c("DEGs_Data_1", "DEGs_Data_2"),
  Total_Genes = c(nrow(processed_data_1), nrow(processed_data_2)),
  Upregulated = c(sum(processed_data_1$status == "Upregulated"), sum(processed_data_2$status == "Upregulated")),
  Downregulated = c(sum(processed_data_1$status == "Downregulated"), sum(processed_data_2$status == "Downregulated")),
  Not_Significant = c(sum(processed_data_1$status == "Not_Significant"), sum(processed_data_2$status == "Not_Significant"))
)

cat("\n====================================================\n")
cat("Final Analysis Summary\n")
cat("====================================================\n")
print(summary_table)


# Task 6: Save the R Workspace for Submission
# -------------------------------------------------------------------
# This saves all key objects from the analysis into one file.

# List of objects to save for a clean workspace file
objects_to_save <- c("classify_gene", "process_dge_data", "processed_data_1", "processed_data_2", "summary_table")

# Save the specified objects to the .RData file
save(list = objects_to_save, file = "SulaimonAraromi_Class_2_Assignment.RData")

cat("\n--> Workspace successfully saved to 'SulaimonAraromi_Class_2_Assignment.RData'.\n")
cat("\nAssignment completed.\n")
