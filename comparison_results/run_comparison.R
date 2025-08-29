#!/usr/bin/env Rscript

# Master comparison script
# Run this after installing R dependencies and extracting master version

source("../compare_cytoconverter.R")

cat("=== RUNNING CYTOCONVERTER COMPARISON ===\n\n")

# Step 1: Prepare input data
input_data <- prepare_input_data("../case_result.txt")

if (is.null(input_data)) {
    stop("Could not prepare input data")
}

cat("Input data prepared:", nrow(input_data), "samples\n\n")

# Step 2: Run current branch analysis
cat("=== RUNNING CURRENT BRANCH ANALYSIS ===\n")
setwd("current_branch")

# Load current branch CytoConverter
tryCatch({
    source("../../modules/cytoconverter.R")
    cat("Successfully loaded current branch CytoConverter\n")
    
    # Run analysis
    result_current <- run_cytoconverter_analysis(input_data, "current", "current_branch")
    
}, error = function(e) {
    cat("Error with current branch:", e$message, "\n")
})

setwd("..")

# Step 3: Run master branch analysis  
cat("\n=== RUNNING MASTER BRANCH ANALYSIS ===\n")
setwd("master_branch")

# Load master branch CytoConverter
tryCatch({
    source("../master_version/cytoscript_master.R")
    cat("Successfully loaded master branch CytoConverter\n")
    
    # Run analysis
    result_master <- run_cytoconverter_analysis(input_data, "master", "master_branch")
    
}, error = function(e) {
    cat("Error with master branch:", e$message, "\n")
})

setwd("..")

# Step 4: Compare results
cat("\n=== COMPARING RESULTS ===\n")
compare_results("current_branch/current_results.txt", 
                "master_branch/master_results.txt", 
                "current", "master")

cat("\n=== COMPARISON COMPLETE ===\n")
cat("Results saved in comparison_results/\n")
