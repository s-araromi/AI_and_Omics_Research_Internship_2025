#!/usr/bin/env Rscript

# Title: Microarray Data Preprocessing and Normalization in R
# Author: Sulaimon Araromi
# Date: October 5, 2025
# Description: This script performs quality control before and after normalization,
#              normalizes microarray data using RMA, filters low-intensity probes,
#              and defines target groups based on phenotype information.

# --- 1. Load Required Libraries ---
# Ensure these packages are installed in your R environment.
# If not, install them using BiocManager::install() or install.packages().
# For example: BiocManager::install(c("affy", "Biobase"))
library(affy)
library(Biobase)

# --- 2. Set Up Data Paths ---
# Define the path to the directory containing the CEL files.
# This path is relative to the working directory inside the Docker container.
cel_path <- "/home/ruser/workdir/Downloads/E-GEOD-42252/CELs"

# List all CEL files in the specified directory.
cel_files <- list.files(cel_path, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)

# Check if CEL files were found and print their names.
if (length(cel_files) == 0) {
  stop("No CEL files found in the specified directory. Please check the path and file names.")
} else {
  cat("Found CEL files:\n")
  print(cel_files)
}

# --- 3. Read Microarray Data ---
# Read the CEL files into an AffyBatch object using the 'affy' package.
# This object will contain the raw intensity data.
raw_data <- ReadAffy(filenames = cel_files)

# Extract phenotype data. If actual phenotype data is available, it should be loaded here.
# For this exercise, we create dummy sample names from the CEL file names.
pData(raw_data)$sample_name <- sampleNames(raw_data)

# --- 4. Quality Control Before Normalization ---
cat("\n--- Performing Quality Control BEFORE Normalization ---\n")

# 4.1. Boxplot of raw data intensities
# Visualizes the distribution of raw probe intensities across all arrays.
# Helps identify arrays with unusual intensity distributions (potential outliers).
png("boxplot_raw.png")
boxplot(exprs(raw_data), main = "Raw Data Boxplot", ylab = "Log2 Intensity", col = "lightblue")
dev.off()
cat("Generated: boxplot_raw.png\n")

# 4.2. Histogram of raw data intensities
# Shows the density distribution of raw probe intensities.
png("hist_raw.png")
hist(exprs(raw_data), main = "Raw Data Histogram", xlab = "Log2 Intensity", col = "lightgreen")
dev.off()
cat("Generated: hist_raw.png\n")

# 4.3. Outlier detection before normalization
# A simple method to flag potential outliers based on median intensity.
# Arrays with median intensity outside 1.5*IQR of all medians are flagged.
median_intensities_raw <- apply(exprs(raw_data), 2, median)
IQR_raw <- IQR(median_intensities_raw)
Q1_raw <- quantile(median_intensities_raw, 0.25)
Q3_raw <- quantile(median_intensities_raw, 0.75)
lower_bound_raw <- Q1_raw - 1.5 * IQR_raw
upper_bound_raw <- Q3_raw + 1.5 * IQR_raw
outliers_before_norm <- sum(median_intensities_raw < lower_bound_raw | median_intensities_raw > upper_bound_raw)
cat(paste0("Number of arrays flagged as outliers BEFORE normalization: ", outliers_before_norm, "\n"))

# --- 5. Normalization (RMA) ---
cat("\n--- Performing RMA Normalization ---\n")
# Robust Multi-array Average (RMA) is a widely used method for background correction,
# normalization, and summarization of Affymetrix GeneChip data.
normalized_data <- rma(raw_data)
cat("RMA normalization complete.\n")

# --- 6. Quality Control After Normalization ---
cat("\n--- Performing Quality Control AFTER Normalization ---\n")

# 6.1. Boxplot of normalized data intensities
# Should show more consistent distributions across arrays after normalization.
png("boxplot_normalized.png")
boxplot(exprs(normalized_data), main = "Normalized Data Boxplot", ylab = "Log2 Intensity", col = "lightblue")
dev.off()
cat("Generated: boxplot_normalized.png\n")

# 6.2. Histogram of normalized data intensities
# Should show more similar density distributions.
png("hist_normalized.png")
hist(exprs(normalized_data), main = "Normalized Data Histogram", xlab = "Log2 Intensity", col = "lightgreen")
dev.off()
cat("Generated: hist_normalized.png\n")

# 6.3. Outlier detection after normalization
# Ideally, the number of outliers should decrease or be zero after effective normalization.
median_intensities_norm <- apply(exprs(normalized_data), 2, median)
IQR_norm <- IQR(median_intensities_norm)
Q1_norm <- quantile(median_intensities_norm, 0.25)
Q3_norm <- quantile(median_intensities_norm, 0.75)
lower_bound_norm <- Q1_norm - 1.5 * IQR_norm
upper_bound_norm <- Q3_norm + 1.5 * IQR_norm
outliers_after_norm <- sum(median_intensities_norm < lower_bound_norm | median_intensities_norm > upper_bound_norm)
cat(paste0("Number of arrays flagged as outliers AFTER normalization: ", outliers_after_norm, "\n"))

# Print summary of normalized data (ExpressionSet object)
print(normalized_data)

# --- 7. Filtering to Remove Low-Intensity Probes ---
cat("\n--- Applying Filtering to Remove Low-Intensity Probes ---\n")

# Get expression values from normalized data
expr_values <- exprs(normalized_data)

# Calculate the maximum intensity for each probe across all samples.
# Probes with consistently low expression across all samples are often considered noise.
max_intensities <- apply(expr_values, 1, max)

# Define a threshold for low intensity. Here, we use the 20th percentile of max intensities.
# Probes whose maximum intensity across all samples falls below this threshold are removed.
intensity_threshold <- quantile(max_intensities, 0.20)

# Identify probes to keep (those with max intensity above the threshold).
probes_to_keep <- max_intensities > intensity_threshold

# Filter the normalized data to retain only the selected probes.
filtered_data <- normalized_data[probes_to_keep, ]

# Count and report the number of transcripts before and after filtering.
num_transcripts_before_filter <- nrow(expr_values)
num_transcripts_after_filter <- nrow(exprs(filtered_data))

cat(paste0("Number of transcripts BEFORE filtering: ", num_transcripts_before_filter, "\n"))
cat(paste0("Number of transcripts AFTER filtering: ", num_transcripts_after_filter, "\n"))
cat("Filtering complete.\n")

# --- 8. Define Target Groups and Re-label ---
cat("\n--- Defining Target Groups and Re-labeling ---\n")

# Extract sample names from the normalized data (or filtered_data if subsequent analysis uses it).
# For this step, we use the normalized_data as the phenotype information is associated with it.
sample_names <- sampleNames(normalized_data)

# Initialize a vector to store group labels for each sample.
group_labels <- character(length(sample_names))

# Assign labels based on patterns in sample names.
# Based on the dataset E-GEOD-42252, 'GC' samples are Gastric Cancer and 'PC' samples are Pancreatic Cancer.
# This is a simplified approach; in a real scenario, phenotype data would be loaded from a separate file.
group_labels[grep("GC", sample_names)] <- "Gastric_Cancer"
group_labels[grep("PC", sample_names)] <- "Pancreatic_Cancer"

# Add the newly defined group labels to the phenoData slot of the normalized_data object.
# This integrates the phenotype information directly with the expression data.
pData(normalized_data)$Group <- group_labels

# Print the updated phenoData to verify the new group assignments.
cat("Updated phenoData with target groups:\n")
print(pData(normalized_data))
cat("Target groups defined and re-labeled.\n")

cat("\nScript execution complete. Please review the generated plots and output messages.\n")



# --- 9. Principal Component Analysis (PCA) ---
cat("\n--- Performing Principal Component Analysis (PCA) ---\n")

# PCA is used to reduce the dimensionality of the data and visualize sample relationships.
# It helps identify batch effects, outliers, and clustering of samples based on their expression profiles.

# Ensure that the 'stats' package is loaded for prcomp (usually loaded by default).
# We will also use 'ggplot2' for enhanced visualization.
library(ggplot2)

# Perform PCA on the transpose of the expression matrix from the normalized data.
# We transpose because prcomp expects variables (genes/probes) as columns and observations (samples) as rows.
pca_result <- prcomp(t(exprs(normalized_data)), scale. = TRUE)

# Prepare data for ggplot2
pca_data <- as.data.frame(pca_result$x)

# Add phenotype information to the PCA data for coloring the plot
pca_data$Group <- pData(normalized_data)$Group

# Plot PCA results
png("pca_plot_normalized.png", width = 800, height = 600)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_data)), vjust = -1, hjust = 0.5) +
  labs(title = "PCA Plot of Normalized Microarray Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 2), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 2), "% variance)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
cat("Generated: pca_plot_normalized.png\n")
cat("PCA analysis complete. Please review the generated plot.\n")


