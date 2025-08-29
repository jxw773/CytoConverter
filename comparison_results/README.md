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
