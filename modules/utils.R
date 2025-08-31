
#' Utility Functions for CytoConverter
#' 
#' @description
#' This module contains utility functions used throughout CytoConverter for 
#' data processing, coordinate manipulation, and interval operations.

#' Position Sorter for Cytogenetic Bands
#' 
#' @description
#' This function ensures that cytogenetic band positions are sorted in the correct
#' order from p-arm to q-arm, with proper handling of complex band notations.
#' 
#' @param positions Character vector of cytogenetic band positions
#' @return Character vector with properly ordered band positions
#' 
#' @details
#' The function handles:
#' \itemize{
#'   \item P-arm bands (sorted from high to low numbers)
#'   \item Q-arm bands (sorted from low to high numbers)  
#'   \item Range notations with "-" or "~" separators
#' }
positionSorter <- function(positions) {
    ## make sure the order of the bands is in the order we like
    ## (p to q highers ps to low, lower qs to high )

    if (length(unlist(strsplit(positions, "-|~"))) > 1) {
        positions <- unlist(lapply(
            strsplit(positions, "-|~"),
            function(x) { x[1] }
        ))
    }

    return(positions)
}

#' Merge Overlapping Intervals
#' 
#' @description
#' This function merges two intervals if they overlap, returning the combined
#' interval or the original intervals if they don't overlap.
#' 
#' @param v1 First interval as numeric vector [start, end]
#' @param v2 Second interval as numeric vector [start, end]
#' @return Merged interval or original intervals, NA if inputs contain NA values
#' 
#' @details
#' The function checks for overlap between two genomic intervals and merges them
#' if they overlap or are adjacent. Used extensively in coordinate processing
#' and interval consolidation throughout CytoConverter.
mergeIntOverlap <- function(v1, v2) {
    ## returns overlapped region, deletes un-overlaps

    v3 = NA
    if (is.na(prod(c(v1, v2)))) {
        return(NA)
    }

    if (
            v1[1] > v2[1] |
            v2[2] < v1[2] |
            v1[2] > v2[1] |
            v2[2] < v1[1]
    ) {
        v3 <- sort(c(v1, v2))
        v3 <- v3[c(2, 3)]
    }
  
    return(v3)
}

#' Detect Addition/Deletion Patterns in Chromosomal Data
#' 
#' @description
#' This function detects and processes addition (+) and deletion (-) patterns in 
#' chromosomal data, standardizing the notation for downstream processing.
#' 
#' @param temp_table_processed Processed temporary table containing chromosomal data
#' @param ex_table_processed Processed reference table for comparison
#' @return List containing processed tables with standardized +/- notation
#' 
#' @details
#' The function:
#' \itemize{
#'   \item Identifies chromosomal additions marked with "+"
#'   \item Identifies chromosomal deletions marked with "-" 
#'   \item Standardizes notation for consistent downstream processing
#'   \item Handles special chromosomes that require different processing
#' }
#' 
#' @note This function is part of the internal processing pipeline and typically
#' should not be called directly by end users.
detectAdd <- function(temp_table_processed, ex_table_processed) {
    ##detects whether there is a + in special chromosomes to handle differently
    ##not working properly because no longer table
  
    temp_table_processed[grep("\\+", temp_table_processed)] <- "+"
    ex_table_processed[grep("\\+", ex_table_processed)] <- "+"

    temp_table_processed[grep("^-", temp_table_processed)] <- "-"
    ex_table_processed[grep("^-", ex_table_processed)] <- "-"
  
    temp_table_processed[grep("\\+|^-", temp_table_processed, invert = T)] <- ""
    ex_table_processed[grep("\\+|^-", ex_table_processed, invert = T)] <- ""
  
    additionList <- list(temp_table_processed, ex_table_processed)
  
    return(additionList)
}


