# Example script demonstrating CNV hierarchical clustering analysis
# This script shows how to use the cnv_clustering_analysis.R functions

# Load the clustering analysis functions
source("cnv_clustering_analysis.R")

# Example 1: Basic clustering analysis using actual CNV regions
cat("Running Example 1: Basic clustering analysis\n")
cat("============================================\n")

if (file.exists("cyto_result.txt")) {
  results1 <- cnv_clustering_analysis(
    input_file = "cyto_result.txt",
    output_dir = "example_output",
    use_bins = FALSE,  # Use actual CNV regions
    distance_method = "binary",
    clustering_method = "complete"
  )
  
  cat("\nExample 1 completed. Check 'example_output' directory for results.\n\n")
} else {
  cat("cyto_result.txt not found. Please ensure the example data file exists.\n\n")
}

# Example 2: Clustering analysis using genomic bins
cat("Running Example 2: Clustering with genomic bins\n")
cat("===============================================\n")

if (file.exists("cyto_result.txt")) {
  results2 <- cnv_clustering_analysis(
    input_file = "cyto_result.txt",
    output_dir = "example_output_bins", 
    use_bins = TRUE,     # Use genomic bins
    bin_size = 25000000, # 25MB bins
    distance_method = "binary",
    clustering_method = "ward.D2"  # Different clustering method
  )
  
  cat("\nExample 2 completed. Check 'example_output_bins' directory for results.\n\n")
} else {
  cat("cyto_result.txt not found. Please ensure the example data file exists.\n\n")
}

# Example 3: Analysis with custom data
cat("Running Example 3: Creating sample data and analysis\n")
cat("====================================================\n")

# Create sample CNV data
sample_data <- data.frame(
  Sample.ID = c("Sample1", "Sample1", "Sample2", "Sample2", "Sample3", "Sample3"),
  Chr = c("chr1", "chr2", "chr1", "chr3", "chr2", "chr4"),
  Start = c(1000000, 5000000, 2000000, 1000000, 6000000, 500000),
  End = c(3000000, 8000000, 4000000, 2000000, 9000000, 1500000),
  Type = c("Gain", "Loss", "Loss", "Gain", "Gain", "Loss"),
  Cells.Present = c("unknown", "unknown", "unknown", "unknown", "unknown", "unknown")
)

# Save sample data
write.table(sample_data, "sample_cnv_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Run analysis on sample data
results3 <- cnv_clustering_analysis(
  input_file = "sample_cnv_data.txt",
  output_dir = "example_output_custom",
  use_bins = FALSE,
  distance_method = "euclidean",  # Different distance method
  clustering_method = "average"   # Different clustering method
)

cat("\nExample 3 completed. Check 'example_output_custom' directory for results.\n")

cat("\n=== All Examples Completed ===\n")
cat("Generated outputs:\n")
cat("- example_output/: Basic analysis results\n")
cat("- example_output_bins/: Analysis with genomic bins\n") 
cat("- example_output_custom/: Analysis with custom sample data\n")
cat("- sample_cnv_data.txt: Sample CNV data file\n")