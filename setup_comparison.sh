#!/bin/bash

# Setup script for CytoConverter branch comparison
# This script sets up the environment to compare CytoConverter results between branches

echo "=== CYTOCONVERTER BRANCH COMPARISON SETUP ==="
echo "Setting up comparison between current branch and master (base commit)"

# Create directories for organized output
mkdir -p comparison_results
mkdir -p comparison_results/current_branch
mkdir -p comparison_results/master_branch

echo "Created output directories in comparison_results/"

# Get current branch information
current_branch=$(git branch --show-current)
current_commit=$(git rev-parse HEAD)
base_commit="89d9d18"

echo "Current branch: $current_branch"
echo "Current commit: $current_commit"
echo "Base commit (master equivalent): $base_commit"

# Save current state information
cat > comparison_results/branch_info.txt << EOF
CytoConverter Branch Comparison Setup
Date: $(date)
Current Branch: $current_branch
Current Commit: $current_commit
Base Commit: $base_commit

Files to compare:
- Input: case_result.txt ($(wc -l < case_result.txt) lines)
- Current branch CytoConverter: modules/cytoconverter.R
- Master branch CytoConverter: To be extracted from base commit

Comparison Plan:
1. Run CytoConverter with case_result.txt on current branch
2. Extract master version and run CytoConverter with case_result.txt  
3. Compare outputs and identify differences
EOF

echo "Saved branch information to comparison_results/branch_info.txt"

# Create a script to extract master version
cat > comparison_results/extract_master.sh << 'EOF'
#!/bin/bash

# Extract master version for comparison
echo "Extracting master version from base commit..."

# Create temporary directory for master version
mkdir -p master_version
cd master_version

# Extract files from base commit
git show 89d9d18:cytoscript_orig.R > cytoscript_master.R
git show 89d9d18:case_result.txt > case_result_master.txt
git show 89d9d18:README.md > README_master.md

# Copy build files (these should be the same)
cp -r ../Builds .

echo "Master version extracted to master_version/"
echo "Key files:"
echo "  - cytoscript_master.R (main CytoConverter function)"
echo "  - case_result_master.txt (input data)"
echo "  - Builds/ (reference genome builds)"

cd ..
EOF

chmod +x comparison_results/extract_master.sh

echo "Created extraction script: comparison_results/extract_master.sh"

# Create R script template for running the actual comparison
cat > comparison_results/run_comparison.R << 'EOF'
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
EOF

chmod +x comparison_results/run_comparison.R

echo "Created comparison runner: comparison_results/run_comparison.R"

# Create README for the comparison
cat > comparison_results/README.md << 'EOF'
# CytoConverter Branch Comparison

This directory contains the setup and results for comparing CytoConverter functionality between the current branch and master branch using `case_result.txt` as input.

## Files

- `branch_info.txt` - Information about branches and commits being compared
- `extract_master.sh` - Script to extract master version files
- `run_comparison.R` - R script to run the actual comparison
- `current_branch/` - Results from current branch CytoConverter
- `master_branch/` - Results from master branch CytoConverter  
- `master_version/` - Extracted files from master branch

## Usage

1. **Install R dependencies** (requires internet connection):
   ```bash
   sudo R -e "install.packages(c('stringr', 'stringi', 'DescTools', 'dplyr'), repos='https://cran.r-project.org')"
   ```

2. **Extract master version**:
   ```bash
   ./extract_master.sh
   ```

3. **Run comparison**:
   ```bash
   Rscript run_comparison.R
   ```

## Expected Output

- `current_results.txt` - CytoConverter results from current branch
- `master_results.txt` - CytoConverter results from master branch
- `comparison_current_vs_master.txt` - Detailed difference report
- Console output showing summary statistics and key differences

## Input Data

The comparison uses `case_result.txt` which contains 2116 rows of karyotype data with columns:
- `Refno` - Reference number
- `KaryShort` - Karyotype string (e.g., "47,XX,+8", "45,X,-Y")

The script converts this to the format expected by CytoConverter (Sample ID, Karyotype).
EOF

echo "Created README: comparison_results/README.md"

# Display summary
echo ""
echo "=== SETUP COMPLETE ==="
echo "Comparison framework is ready. Next steps:"
echo ""
echo "1. Install R dependencies (requires internet):"
echo "   sudo R -e \"install.packages(c('stringr', 'stringi', 'DescTools', 'dplyr'))\""
echo ""
echo "2. Extract master version:"
echo "   cd comparison_results && ./extract_master.sh"
echo ""
echo "3. Run the comparison:"
echo "   cd comparison_results && Rscript run_comparison.R"
echo ""
echo "4. Review results in comparison_results/ directory"
echo ""
echo "All setup files are in comparison_results/"