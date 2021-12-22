
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


