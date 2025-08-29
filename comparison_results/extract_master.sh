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
