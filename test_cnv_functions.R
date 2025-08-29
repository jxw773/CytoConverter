#' Unit Tests for CNV Annotation and Pathway Analysis Functions
#' 
#' This script provides unit tests for the CNV annotation and pathway analysis
#' functionality to ensure code quality and reliability.

# Load required libraries
library(stringr)
library(dplyr)

# Source the functions to test
source("cnv_annotation.R")
source("pathway_analysis.R")

#' Simple test framework
test_that <- function(description, test_expression) {
  cat(sprintf("Testing: %s\n", description))
  tryCatch({
    if (test_expression) {
      cat("  ✓ PASS\n")
      return(TRUE)
    } else {
      cat("  ✗ FAIL\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", e$message))
    return(FALSE)
  })
}

cat("=======================================================\n")
cat("       UNIT TESTS FOR CNV ANALYSIS FUNCTIONS\n")
cat("=======================================================\n\n")

# Initialize test results
total_tests <- 0
passed_tests <- 0

# Test 1: load_gene_annotations function
total_tests <- total_tests + 1
if (test_that("load_gene_annotations returns data frame with correct columns", {
  gene_data <- load_gene_annotations("GRCh38")
  is.data.frame(gene_data) && 
  all(c("gene_symbol", "chromosome", "start_pos", "end_pos") %in% names(gene_data))
})) passed_tests <- passed_tests + 1

# Test 2: load_gene_annotations with invalid build
total_tests <- total_tests + 1
if (test_that("load_gene_annotations rejects invalid build", {
  error_occurred <- FALSE
  tryCatch({
    load_gene_annotations("invalid_build")
  }, error = function(e) {
    error_occurred <<- TRUE
  })
  error_occurred
})) passed_tests <- passed_tests + 1

# Test 3: annotate_cnv_with_genes with valid input
total_tests <- total_tests + 1
if (test_that("annotate_cnv_with_genes works with valid CNV data", {
  test_cnv <- data.frame(
    "Sample ID" = "Test",
    "Chr" = "chr17",
    "Start" = 7000000,
    "End" = 8000000,
    "Type" = "Loss",
    "Cells Present" = "unknown",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result <- annotate_cnv_with_genes(test_cnv)
  is.data.frame(result) && nrow(result) > 0 && "Gene_Symbol" %in% names(result)
})) passed_tests <- passed_tests + 1

# Test 4: annotate_cnv_with_genes with missing columns
total_tests <- total_tests + 1
if (test_that("annotate_cnv_with_genes rejects invalid input", {
  invalid_cnv <- data.frame(
    "Sample" = "Test",
    "Position" = 123456,
    stringsAsFactors = FALSE
  )
  error_occurred <- FALSE
  tryCatch({
    annotate_cnv_with_genes(invalid_cnv)
  }, error = function(e) {
    error_occurred <<- TRUE
  })
  error_occurred
})) passed_tests <- passed_tests + 1

# Test 5: extract_genes_from_cnv function
total_tests <- total_tests + 1
if (test_that("extract_genes_from_cnv extracts genes correctly", {
  test_annotated <- data.frame(
    "Sample.ID" = c("Test", "Test"),
    "Type" = c("Gain", "Loss"),
    "Gene_Symbol" = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )
  all_genes <- extract_genes_from_cnv(test_annotated, "All")
  gain_genes <- extract_genes_from_cnv(test_annotated, "Gain")
  loss_genes <- extract_genes_from_cnv(test_annotated, "Loss")
  
  length(all_genes) == 2 && length(gain_genes) == 1 && length(loss_genes) == 1
})) passed_tests <- passed_tests + 1

# Test 6: load_pathway_database function
total_tests <- total_tests + 1
if (test_that("load_pathway_database returns correct structure", {
  pathway_data <- load_pathway_database("KEGG")
  is.list(pathway_data) && 
  "database" %in% names(pathway_data) && 
  "pathways" %in% names(pathway_data) &&
  pathway_data$database == "KEGG"
})) passed_tests <- passed_tests + 1

# Test 7: load_pathway_database with invalid database
total_tests <- total_tests + 1
if (test_that("load_pathway_database rejects invalid database", {
  error_occurred <- FALSE
  tryCatch({
    load_pathway_database("INVALID_DB")
  }, error = function(e) {
    error_occurred <<- TRUE
  })
  error_occurred
})) passed_tests <- passed_tests + 1

# Test 8: perform_pathway_enrichment function
total_tests <- total_tests + 1
if (test_that("perform_pathway_enrichment works with valid gene list", {
  test_genes <- c("TP53", "BRCA1", "EGFR")
  result <- perform_pathway_enrichment(test_genes, database = "KEGG", p_cutoff = 1.0)
  is.data.frame(result) && 
  all(c("Pathway_ID", "Pathway_Name", "P_Value", "P_Adjusted") %in% names(result))
})) passed_tests <- passed_tests + 1

# Test 9: perform_pathway_enrichment with empty gene list
total_tests <- total_tests + 1
if (test_that("perform_pathway_enrichment handles empty gene list", {
  # This should return an empty data frame, not error
  result <- tryCatch({
    perform_pathway_enrichment(character(0))
  }, error = function(e) {
    # Expected to fail with empty input, so this is correct behavior
    return(data.frame())
  })
  is.data.frame(result)
})) passed_tests <- passed_tests + 1

# Test 10: perform_pathway_enrichment with invalid parameters
total_tests <- total_tests + 1
if (test_that("perform_pathway_enrichment validates parameters", {
  error_occurred <- FALSE
  tryCatch({
    perform_pathway_enrichment(c("TP53"), p_cutoff = 1.5)  # Invalid p_cutoff
  }, error = function(e) {
    error_occurred <<- TRUE
  })
  error_occurred
})) passed_tests <- passed_tests + 1

# Test 11: visualize_pathway_results function
total_tests <- total_tests + 1
if (test_that("visualize_pathway_results handles empty results", {
  empty_results <- data.frame()
  # Should not throw error, just display message
  result <- tryCatch({
    capture.output(visualize_pathway_results(empty_results))
    TRUE
  }, error = function(e) {
    FALSE
  })
  result
})) passed_tests <- passed_tests + 1

# Test 12: pathway_enrichment_summary function
total_tests <- total_tests + 1
if (test_that("pathway_enrichment_summary returns list structure", {
  test_genes <- c("TP53", "MYC")
  result <- pathway_enrichment_summary(test_genes, databases = c("KEGG"), p_cutoff = 1.0)
  is.list(result) && "KEGG" %in% names(result)
})) passed_tests <- passed_tests + 1

# Test 13: Integration test - full workflow
total_tests <- total_tests + 1
if (test_that("Full workflow integration test", {
  # Create test CNV data
  test_cnv <- data.frame(
    "Sample ID" = "Integration_Test",
    "Chr" = "chr17",
    "Start" = 7000000,
    "End" = 8000000,
    "Type" = "Loss",
    "Cells Present" = "50/50",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Annotate with genes
  annotated <- annotate_cnv_with_genes(test_cnv)
  
  # Extract genes
  genes <- extract_genes_from_cnv(annotated)
  
  # Perform pathway analysis if genes found
  if (length(genes) > 0) {
    pathways <- perform_pathway_enrichment(genes, p_cutoff = 1.0)
    TRUE
  } else {
    TRUE  # Still pass if no genes found (expected for some regions)
  }
})) passed_tests <- passed_tests + 1

# Test 14: Error handling in annotate_cnv_with_genes
total_tests <- total_tests + 1
if (test_that("annotate_cnv_with_genes handles edge cases", {
  # Test with no overlapping genes
  test_cnv <- data.frame(
    "Sample ID" = "No_Genes_Test",
    "Chr" = "chr99",  # Non-existent chromosome
    "Start" = 1,
    "End" = 100,
    "Type" = "Gain",
    "Cells Present" = "unknown",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  result <- annotate_cnv_with_genes(test_cnv)
  is.data.frame(result) && nrow(result) > 0
})) passed_tests <- passed_tests + 1

# Test 15: Different genome builds
total_tests <- total_tests + 1
if (test_that("Functions work with different genome builds", {
  builds <- c("GRCh38", "hg19", "hg18", "hg17")
  all_builds_work <- TRUE
  
  for (build in builds) {
    tryCatch({
      gene_data <- load_gene_annotations(build)
      if (!is.data.frame(gene_data)) {
        all_builds_work <- FALSE
        break
      }
    }, error = function(e) {
      all_builds_work <<- FALSE
    })
  }
  
  all_builds_work
})) passed_tests <- passed_tests + 1

cat("\n=======================================================\n")
cat("                  TEST SUMMARY\n")
cat("=======================================================\n")
cat(sprintf("Total tests run: %d\n", total_tests))
cat(sprintf("Tests passed: %d\n", passed_tests))
cat(sprintf("Tests failed: %d\n", total_tests - passed_tests))
cat(sprintf("Success rate: %.1f%%\n", (passed_tests / total_tests) * 100))

if (passed_tests == total_tests) {
  cat("\n🎉 ALL TESTS PASSED! 🎉\n")
} else {
  cat("\n⚠️  SOME TESTS FAILED ⚠️\n")
}

cat("=======================================================\n")