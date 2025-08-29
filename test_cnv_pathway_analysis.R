#' Test and Demonstration Script for CNV Annotation and Pathway Analysis
#' 
#' This script demonstrates the functionality of the CNV annotation and pathway analysis
#' functions integrated with CytoConverter.

# Load required libraries
library(stringr)
library(dplyr)

# Source the new functions
source("cnv_annotation.R")
source("pathway_analysis.R")

cat("=======================================================\n")
cat("  CNV ANNOTATION AND PATHWAY ANALYSIS DEMONSTRATION\n")
cat("=======================================================\n\n")

# Create sample CNV data (simulating CytoConverter output)
cat("1. Creating sample CNV data...\n")
sample_cnv_data <- data.frame(
  "Sample ID" = c("Patient_A", "Patient_A", "Patient_B", "Patient_C", "Patient_C"),
  "Chr" = c("chr17", "chr7", "chr8", "chr3", "chr9"),
  "Start" = c(7000000, 55000000, 127000000, 179000000, 21900000),
  "End" = c(8000000, 56000000, 128000000, 180000000, 22000000),
  "Type" = c("Loss", "Gain", "Gain", "Gain", "Loss"),
  "Cells Present" = c("45 of 50", "30 of 50", "unknown", "40 of 50", "25 of 50"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(sample_cnv_data)
cat("\n")

# Test CNV annotation
cat("2. Annotating CNV regions with genes...\n")
annotated_cnv <- annotate_cnv_with_genes(sample_cnv_data, build = "GRCh38")

cat("Annotated CNV results:\n")
print(annotated_cnv[, c("Sample.ID", "Chr", "Type", "Gene_Symbol", "Gene_Description", "Overlap_Percentage")])
cat("\n")

# Extract genes for pathway analysis
cat("3. Extracting genes for pathway analysis...\n")

# Extract all genes
all_genes <- extract_genes_from_cnv(annotated_cnv, cnv_type = "All")
cat("All affected genes:", paste(all_genes, collapse = ", "), "\n")

# Extract only gained genes
gained_genes <- extract_genes_from_cnv(annotated_cnv, cnv_type = "Gain")
cat("Gained genes:", paste(gained_genes, collapse = ", "), "\n")

# Extract only lost genes
lost_genes <- extract_genes_from_cnv(annotated_cnv, cnv_type = "Loss")
cat("Lost genes:", paste(lost_genes, collapse = ", "), "\n\n")

# Perform pathway enrichment analysis
cat("4. Performing pathway enrichment analysis...\n\n")

if (length(all_genes) > 0) {
  # Analyze KEGG pathways
  cat("KEGG Pathway Analysis:\n")
  cat("====================\n")
  kegg_results <- perform_pathway_enrichment(all_genes, database = "KEGG", p_cutoff = 0.1)
  
  if (nrow(kegg_results) > 0) {
    visualize_pathway_results(kegg_results, top_n = 5)
  } else {
    cat("No significant KEGG pathways found.\n\n")
  }
  
  # Analyze Reactome pathways
  cat("Reactome Pathway Analysis:\n")
  cat("=========================\n")
  reactome_results <- perform_pathway_enrichment(all_genes, database = "Reactome", p_cutoff = 0.1)
  
  if (nrow(reactome_results) > 0) {
    visualize_pathway_results(reactome_results, top_n = 3)
  } else {
    cat("No significant Reactome pathways found.\n\n")
  }
  
  # Comprehensive summary across databases
  cat("5. Comprehensive pathway analysis summary...\n")
  summary_results <- pathway_enrichment_summary(all_genes, 
                                                databases = c("KEGG", "Reactome", "GO_BP"), 
                                                p_cutoff = 0.1)
  
  cat("\nSummary of results:\n")
  for (db_name in names(summary_results)) {
    result <- summary_results[[db_name]]
    cat(sprintf("- %s: %d significant pathways\n", db_name, nrow(result)))
  }
  
} else {
  cat("No genes found for pathway analysis.\n")
}

cat("\n=======================================================\n")
cat("  DEMONSTRATION COMPLETE\n")
cat("=======================================================\n")

# Clean up
rm(sample_cnv_data, annotated_cnv, all_genes, gained_genes, lost_genes)
if (exists("kegg_results")) rm(kegg_results)
if (exists("reactome_results")) rm(reactome_results)
if (exists("summary_results")) rm(summary_results)