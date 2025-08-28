
#' Utility Functions for CytoConverter
#' 
#' @description
#' This module contains utility functions used throughout CytoConverter for
#' data processing, sorting, and coordinate manipulation.

#' Position Sorter Function
#' 
#' @description
#' Ensures the order of cytogenetic bands follows the expected sequence
#' from p-arm (short arm) to q-arm (long arm), with higher p values to lower,
#' and lower q values to higher.
#' 
#' @param positions Character vector of cytogenetic band positions
#' @return Character vector of sorted positions
#' 
#' @details
#' The function handles position ranges indicated by "-" or "~" separators
#' by taking the first position in each range for sorting purposes.
#' 
#' @examples
#' \dontrun{
#' positions <- c("p23", "q11", "p21-p22")
#' sorted_pos <- positionSorter(positions)
#' }
#' 
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
#' Returns the overlapped region between two genomic intervals,
#' removing non-overlapping portions.
#' 
#' @param v1 Numeric vector of length 2 representing first interval [start, end]
#' @param v2 Numeric vector of length 2 representing second interval [start, end]
#' @return Numeric vector of length 2 representing overlapped region, or NA if no overlap
#' 
#' @details
#' The function calculates the intersection of two genomic intervals.
#' If the intervals do not overlap, it returns the sorted middle values
#' to represent the gap between them.
#' 
#' @examples
#' \dontrun{
#' interval1 <- c(1000, 2000)
#' interval2 <- c(1500, 2500)
#' overlap <- mergeIntOverlap(interval1, interval2)
#' # Returns c(1500, 2000)
#' }
#' 
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

#' Detect Addition Patterns
#' 
#' @description
#' Detects whether there are "+" or "-" symbols in chromosome tables
#' to handle gain and loss patterns appropriately.
#' 
#' @param temp_table_processed Processed temporary table
#' @param ex_table_processed Processed exclusion table
#' @return List containing processed tables with standardized addition/deletion markers
#' 
#' @details
#' This function standardizes gain (+) and loss (-) indicators across
#' chromosome tables by:
#' \itemize{
#'   \item Converting any "+" patterns to standard "+"
#'   \item Converting any leading "-" patterns to standard "-"
#'   \item Clearing entries that don't match gain/loss patterns
#' }
#' 
#' @examples
#' \dontrun{
#' temp_table <- c("+21", "-5", "normal")
#' ex_table <- c("+X", "-Y", "other")
#' result <- detectAdd(temp_table, ex_table)
#' }
#' 
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


