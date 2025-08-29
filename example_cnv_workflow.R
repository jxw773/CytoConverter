#' Example Usage of CytoConverter CNV Analysis Workflow
#' 
#' This script demonstrates how to use the integrated CNV analysis workflow
#' with various types of input data.

# Load required libraries
library(stringr)
library(dplyr)

# Source the integrated workflow
source("cnv_analysis_workflow.R")

cat("=========================================================\n")
cat("  CYTOCONVERTER CNV ANALYSIS WORKFLOW EXAMPLES\n")
cat("=========================================================\n\n")

# Example 1: Single karyotype string
cat("EXAMPLE 1: Single karyotype analysis\n")
cat("=====================================\n")

single_karyotype <- "47,XX,+8"
cat(sprintf("Input karyotype: %s\n\n", single_karyotype))

result1 <- cyto_cnv_analysis(
  input_data = single_karyotype,
  build = "GRCh38",
  perform_pathways = TRUE,
  pathway_databases = c("KEGG", "Reactome"),
  save_results = FALSE
)

# Display summary
display_cnv_summary(result1, show_pathways = TRUE, max_pathways = 3)

cat("\n\n")

# Example 2: Multiple samples with mock data
cat("EXAMPLE 2: Multiple sample analysis\n")
cat("===================================\n")

# Create sample data that simulates CytoConverter output
multi_sample_data <- data.frame(
  "Sample ID" = c("Patient_001", "Patient_001", "Patient_002", "Patient_003", "Patient_003"),
  "Chr" = c("chr17", "chr8", "chr7", "chr3", "chr9"),
  "Start" = c(7000000, 127000000, 55000000, 179000000, 21900000),
  "End" = c(8000000, 128000000, 56000000, 180000000, 22000000),
  "Type" = c("Loss", "Gain", "Gain", "Gain", "Loss"),
  "Cells Present" = c("45/50", "30/50", "unknown", "40/50", "25/50"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input CNV data:\n")
print(multi_sample_data)
cat("\n")

result2 <- cyto_cnv_analysis(
  input_data = multi_sample_data,
  build = "GRCh38",
  perform_pathways = TRUE,
  pathway_databases = c("KEGG", "GO_BP"),
  save_results = FALSE
)

# Display summary
display_cnv_summary(result2, show_pathways = TRUE, max_pathways = 5)

cat("\n\n")

# Example 3: Analysis with file output
cat("EXAMPLE 3: Analysis with file output\n")
cat("====================================\n")

# Use the same data but save results
result3 <- cyto_cnv_analysis(
  input_data = multi_sample_data,
  build = "GRCh38",
  perform_pathways = TRUE,
  pathway_databases = c("KEGG", "Reactome", "GO_BP"),
  save_results = TRUE,
  output_prefix = "example_analysis"
)

# List created files
created_files <- list.files(pattern = "example_analysis_.*\\.txt$")
if (length(created_files) > 0) {
  cat("Files created:\n")
  for (file in created_files) {
    cat(sprintf("  - %s\n", file))
  }
} else {
  cat("No output files were created (this is expected in the demo with mock data)\n")
}

cat("\n\n")

# Example 4: Focus on specific CNV types
cat("EXAMPLE 4: Analysis by CNV type\n")
cat("===============================\n")

# Extract genes by CNV type
if (nrow(result2$annotated_cnv) > 0) {
  gains_only <- result2$gene_lists$gained_genes
  losses_only <- result2$gene_lists$lost_genes
  
  cat(sprintf("Genes in GAINS (%d): %s\n", length(gains_only), paste(gains_only, collapse = ", ")))
  cat(sprintf("Genes in LOSSES (%d): %s\n", length(losses_only), paste(losses_only, collapse = ", ")))
  
  # Analyze pathways for gains only
  if (length(gains_only) > 0) {
    cat("\nPathway analysis for GAINED genes only:\n")
    source("pathway_analysis.R")
    gains_pathways <- perform_pathway_enrichment(gains_only, database = "KEGG", p_cutoff = 0.1)
    
    if (nrow(gains_pathways) > 0) {
      visualize_pathway_results(gains_pathways, top_n = 3)
    } else {
      cat("No significant pathways found for gained genes.\n")
    }
  }
  
  # Analyze pathways for losses only
  if (length(losses_only) > 0) {
    cat("\nPathway analysis for LOST genes only:\n")
    losses_pathways <- perform_pathway_enrichment(losses_only, database = "KEGG", p_cutoff = 0.1)
    
    if (nrow(losses_pathways) > 0) {
      visualize_pathway_results(losses_pathways, top_n = 3)
    } else {
      cat("No significant pathways found for lost genes.\n")
    }
  }
}

cat("\n\n")

# Example 5: Different genome builds
cat("EXAMPLE 5: Different genome builds\n")
cat("==================================\n")

builds_to_test <- c("GRCh38", "hg19")

for (build in builds_to_test) {
  cat(sprintf("Testing with genome build: %s\n", build))
  
  # Simple analysis
  small_data <- data.frame(
    "Sample ID" = "Test_Sample",
    "Chr" = "chr17",
    "Start" = 7000000,
    "End" = 8000000,
    "Type" = "Loss",
    "Cells Present" = "unknown",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  build_result <- cyto_cnv_analysis(
    input_data = small_data,
    build = build,
    perform_pathways = FALSE,  # Skip pathways for speed
    save_results = FALSE
  )
  
  cat(sprintf("  Genes found: %s\n", paste(build_result$gene_lists$all_genes, collapse = ", ")))
  cat("\n")
}

cat("=========================================================\n")
cat("              ALL EXAMPLES COMPLETED\n")
cat("=========================================================\n")

# Clean up example files
example_files <- list.files(pattern = "example_analysis_.*\\.txt$")
if (length(example_files) > 0) {
  file.remove(example_files)
  cat("Cleaned up example output files.\n")
}

cat("\nWorkflow demonstration complete!\n")
cat("Use cyto_cnv_analysis() for your own data analysis.\n")