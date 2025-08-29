#!/usr/bin/env Rscript

# CytoConverter Comparison Script
# This script compares CytoConverter results between the current branch and master branch
# using case_result.txt as input

# Function to prepare input data from case_result.txt
prepare_input_data <- function(file_path = "case_result.txt") {
    cat("Reading input data from", file_path, "...\n")
    
    # Read the case_result.txt file
    case_data <- read.delim(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # Remove quotes from KaryShort column
    case_data$KaryShort <- gsub('"', '', case_data$KaryShort)
    
    # Create input matrix for CytoConverter (sample_name, karyotype)
    input_matrix <- matrix(c(
        paste("Sample", case_data$Refno, sep="_"),
        case_data$KaryShort
    ), ncol=2)
    
    colnames(input_matrix) <- c("SampleID", "Karyotype")
    
    cat("Prepared", nrow(input_matrix), "samples for analysis\n")
    return(input_matrix)
}

# Function to run CytoConverter and save results
run_cytoconverter_analysis <- function(input_data, output_prefix, branch_name) {
    cat("Running CytoConverter analysis for", branch_name, "...\n")
    
    # Check if CytoConverter function is available
    if (!exists("CytoConverter")) {
        cat("Error: CytoConverter function not found. Please source the appropriate script first.\n")
        return(NULL)
    }
    
    # Run CytoConverter
    tryCatch({
        result <- CytoConverter(input_data, build="GRCh38", constitutional=TRUE)
        
        # Save results
        results_file <- paste0(output_prefix, "_results.txt")
        errors_file <- paste0(output_prefix, "_errors.txt")
        
        # Write results table
        if (!is.null(result$Result)) {
            write.table(result$Result, file=results_file, sep='\t', row.names=FALSE, quote=FALSE)
            cat("Results saved to", results_file, "\n")
        }
        
        # Write error log
        if (!is.null(result$Error_log)) {
            write.table(result$Error_log, file=errors_file, sep='\t', row.names=FALSE, quote=FALSE)
            cat("Error log saved to", errors_file, "\n")
        }
        
        return(result)
        
    }, error = function(e) {
        cat("Error running CytoConverter:", e$message, "\n")
        return(NULL)
    })
}

# Function to compare two CytoConverter result files
compare_results <- function(results1_file, results2_file, branch1_name, branch2_name) {
    cat("\n=== COMPARING CYTOCONVERTER RESULTS ===\n")
    cat("Branch 1:", branch1_name, "->", results1_file, "\n")
    cat("Branch 2:", branch2_name, "->", results2_file, "\n")
    
    if (!file.exists(results1_file) || !file.exists(results2_file)) {
        cat("Error: One or both result files do not exist\n")
        return(NULL)
    }
    
    # Read both result files
    results1 <- read.delim(results1_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    results2 <- read.delim(results2_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    
    cat("\nSummary Statistics:\n")
    cat(sprintf("%-20s %10s %10s\n", "Metric", branch1_name, branch2_name))
    cat(sprintf("%-20s %10d %10d\n", "Total rows:", nrow(results1), nrow(results2)))
    cat(sprintf("%-20s %10d %10d\n", "Unique samples:", length(unique(results1$Sample.ID)), length(unique(results2$Sample.ID))))
    
    if ("Type" %in% colnames(results1) && "Type" %in% colnames(results2)) {
        gains1 <- sum(results1$Type == "Gain", na.rm=TRUE)
        gains2 <- sum(results2$Type == "Gain", na.rm=TRUE)
        losses1 <- sum(results1$Type == "Loss", na.rm=TRUE) 
        losses2 <- sum(results2$Type == "Loss", na.rm=TRUE)
        
        cat(sprintf("%-20s %10d %10d\n", "Gains:", gains1, gains2))
        cat(sprintf("%-20s %10d %10d\n", "Losses:", losses1, losses2))
    }
    
    # Find differences
    cat("\n=== DETAILED DIFFERENCES ===\n")
    
    # Create comparison keys
    if (all(c("Sample.ID", "Chr", "Start", "End", "Type") %in% colnames(results1)) && 
        all(c("Sample.ID", "Chr", "Start", "End", "Type") %in% colnames(results2))) {
        
        key1 <- paste(results1$Sample.ID, results1$Chr, results1$Start, results1$End, results1$Type, sep="|")
        key2 <- paste(results2$Sample.ID, results2$Chr, results2$Start, results2$End, results2$Type, sep="|")
        
        only_in_1 <- setdiff(key1, key2)
        only_in_2 <- setdiff(key2, key1)
        
        cat("Records only in", branch1_name, ":", length(only_in_1), "\n")
        cat("Records only in", branch2_name, ":", length(only_in_2), "\n")
        
        if (length(only_in_1) > 0) {
            cat("\nFirst 10 records only in", branch1_name, ":\n")
            indices1 <- which(key1 %in% only_in_1)[1:min(10, length(only_in_1))]
            print(results1[indices1, ])
        }
        
        if (length(only_in_2) > 0) {
            cat("\nFirst 10 records only in", branch2_name, ":\n")
            indices2 <- which(key2 %in% only_in_2)[1:min(10, length(only_in_2))]
            print(results2[indices2, ])
        }
        
        # Save detailed comparison
        comparison_file <- paste0("comparison_", branch1_name, "_vs_", branch2_name, ".txt")
        
        sink(comparison_file)
        cat("=== CYTOCONVERTER COMPARISON REPORT ===\n")
        cat("Date:", date(), "\n")
        cat("Branch 1:", branch1_name, "\n")
        cat("Branch 2:", branch2_name, "\n\n")
        
        cat("Summary:\n")
        cat("Total records in", branch1_name, ":", nrow(results1), "\n")
        cat("Total records in", branch2_name, ":", nrow(results2), "\n")
        cat("Records only in", branch1_name, ":", length(only_in_1), "\n")
        cat("Records only in", branch2_name, ":", length(only_in_2), "\n")
        
        if (length(only_in_1) > 0) {
            cat("\n=== RECORDS ONLY IN", toupper(branch1_name), "===\n")
            write.table(results1[key1 %in% only_in_1, ], sep='\t', row.names=FALSE, quote=FALSE)
        }
        
        if (length(only_in_2) > 0) {
            cat("\n=== RECORDS ONLY IN", toupper(branch2_name), "===\n")
            write.table(results2[key2 %in% only_in_2, ], sep='\t', row.names=FALSE, quote=FALSE)
        }
        
        sink()
        
        cat("\nDetailed comparison saved to:", comparison_file, "\n")
    }
}

# Main execution function
main <- function() {
    cat("=== CYTOCONVERTER BRANCH COMPARISON ===\n")
    cat("This script compares CytoConverter results between branches using case_result.txt\n\n")
    
    # Prepare input data
    input_data <- prepare_input_data("case_result.txt")
    
    if (is.null(input_data)) {
        cat("Error: Could not prepare input data\n")
        return()
    }
    
    # Display sample of input data
    cat("\nSample of input data:\n")
    print(head(input_data, 5))
    
    cat("\n=== INSTRUCTIONS FOR COMPARISON ===\n")
    cat("To complete the comparison, follow these steps:\n\n")
    
    cat("1. CURRENT BRANCH ANALYSIS:\n")
    cat("   - Source the current CytoConverter function\n")
    cat("   - Run: source('modules/cytoconverter.R') # or appropriate file\n")
    cat("   - Run: result_current <- run_cytoconverter_analysis(input_data, 'current_branch', 'current')\n\n")
    
    cat("2. MASTER BRANCH ANALYSIS:\n")
    cat("   - Check out master branch or base commit: git checkout 89d9d18\n")
    cat("   - Source the master version CytoConverter function\n")
    cat("   - Run: result_master <- run_cytoconverter_analysis(input_data, 'master_branch', 'master')\n\n")
    
    cat("3. COMPARE RESULTS:\n")
    cat("   - Run: compare_results('current_branch_results.txt', 'master_branch_results.txt', 'current', 'master')\n\n")
    
    cat("Note: This script framework is ready to run once R dependencies are installed.\n")
    cat("Required packages: stringr, stringi, DescTools, dplyr\n")
}

# Execute main function
if (!interactive()) {
    main()
}