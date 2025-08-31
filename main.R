#!/usr/bin/env Rscript

#' CytoConverter Main Entry Point
#' 
#' @description
#' This script provides the main command-line interface for CytoConverter.
#' It handles argument parsing, input/output operations, and calls the core
#' CytoConverter functionality through the modular implementation.
#' 
#' @usage
#' ./main.R --input input.txt --output results.txt [--log errors.txt]
#' 
#' @details
#' The script uses optparse for command-line argument handling and supports:
#' - Input file specification (or stdin if not provided)
#' - Output file specification (required)
#' - Optional log file for errors and warnings
#' 
#' Processing uses constitutional=FALSE and guess=TRUE by default for
#' somatic cytogenetic analysis with ambiguous region guessing enabled.

options(scipen = 999)  # Disable scientific notation for coordinate output


# Load modules
# Use the modular implementation for better maintainability
mod_cytoconverter <- modules::use('modules/cytoconverter.R')


# Parse arguments
# Set up command-line argument parsing with optparse for robust CLI interface

option_list = list(
    optparse::make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        help="Input file containing sample names and karyotypes (tab-delimited). If not specified, reads from stdin.",
        metavar="character"
    ),
    optparse::make_option(
        c("-o", "--output"),
        type="character",
        default=NULL,
        help="Output file for CytoConverter results (required). Contains genomic coordinates of aberrations.",
        metavar="character"
    ),
    optparse::make_option(
        c("-l", "--log"),
        type="character",
        default=NULL,
        help="Log file for errors and warnings (optional). Captures parsing issues and problematic karyotypes.",
        metavar="character"
    )
)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$output)) {
    optparse::print_help(opt_parser)
    stop("Missing required argument: --output -> output file name")
}



# Call CytoConverter
# Process input data and generate results using core CytoConverter functionality

input_m <- if (is.null(opt$input)) {
    # Read from stdin if no input specified - enables piping from other tools
    # Expected format: sample_name<TAB>karyotype (one per line)
    as.matrix(as.data.frame(read.delim(file=file("stdin"), header=F, sep="\t")))
} else {
    # Read from specified input file with same tab-delimited format
    as.matrix(as.data.frame(read.delim(file=opt$input, header=F, sep="\t")))
}

# Run CytoConverter with default settings optimized for somatic analysis
# constitutional=F: Optimized for somatic (cancer) karyotypes rather than germline
# guess=T: Enable guessing for ambiguous cytogenetic regions (improves parsing success)
# These defaults work well for most cancer cytogenetics datasets
result <- mod_cytoconverter$CytoConverter(input_m, constitutional=F, guess=T)

# Write output files
# Results are structured as: [[1]] = main results table, [[2]] = error/warning log

# Main results table (required output)
# Format: Sample, Chr, Start, End, Type, Percentage
write.table(result[[1]], file=opt$output, quote=FALSE, sep='\t', col.names=F, row.names=F)

# Error/warning log (optional output for troubleshooting)
# Contains information about parsing failures, ambiguous karyotypes, etc.
if (!is.null(opt$log)) {
    write.table(result[[2]], file=opt$log, quote=FALSE, sep='\t', col.names=F, row.names=F)
}

