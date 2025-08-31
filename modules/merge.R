
#' Merge Functions for CytoConverter
#' 
#' @description
#' This module contains functions for merging overlapping chromosomal intervals 
#' and handling complex gain/loss regions. These functions maintain non-overlapping
#' data structures while preserving the hierarchical nature of nested aberrations.
#' 
#' @details
#' The merge system uses a hash-based data structure to efficiently track:
#' \itemize{
#'   \item Non-overlapping genomic intervals
#'   \item Gain and loss regions within each interval
#'   \item Plus-loss events (+Loss) for specialized processing
#'   \item Order preservation for original aberration sequence
#' }

#' Insert Section into Merge Data Structure
#' 
#' @description
#' This function inserts chromosomal regions into the merge data structure h_, which
#' tracks all non-overlapping gain/loss sections within genomic regions. It handles
#' the complex logic of splitting existing intervals when new regions overlap.
#' 
#' @param h_ Hash data structure containing existing genomic intervals with the format:
#'   Start coordinate (key) -> hash containing:
#'     \itemize{
#'       \item End: end coordinate for this section
#'       \item Gain: list of end coordinates of "Gain" regions
#'       \item Loss: list of end coordinates of "Loss" regions  
#'       \item "+Loss": list of end coordinates of "+Loss" regions
#'     }
#' @param start Numeric start coordinate of the new region to insert
#' @param end Numeric end coordinate of the new region to insert
#' @param type Character string specifying the type of aberration ("Gain", "Loss", or "+Loss")
#' 
#' @return Modified hash data structure with the new region inserted and existing
#'   intervals split as necessary to maintain non-overlapping structure
#' 
#' @details
#' The function performs three main operations:
#' \enumerate{
#'   \item **Split at start position**: If the start coordinate falls within an existing
#'     interval, split that interval at the start position
#'   \item **Split at end position**: If the end coordinate falls within an existing  
#'     interval, split that interval at the end position
#'   \item **Update overlapping intervals**: For all intervals completely contained
#'     within the new region, add the aberration type to their respective lists
#' }
#' 
#' This approach ensures that overlapping gains and losses are properly tracked
#' without creating conflicting interval boundaries.
#' 
#' @examples
#' \dontrun{
#' # Create new hash structure
#' h <- hash::hash()
#' h["100"] <- hash::hash(End=200, Gain=list(), Loss=list(), "+Loss"=list())
#' 
#' # Insert a gain from 150-250
#' h <- insertSection(h, 150, 250, "Gain")
#' # This splits the 100-200 interval and creates appropriate gain annotations
#' }
insertSection <- function(h_, start, end, type) {
    # This function inserts regions into the data structure h_, which
    # keeps track of all non-overlapping gain/loss sections within regions
    
    # h_:
    #   Start coord:
    #      End: end coordinate for this section
    #      Gain: list of end coords of "Gain" regions, ordered by appearance in original list
    #      Loss: list of end coords of "Loss" regions, ordered by appearance in original list
    #   next start coord:
    #    ...
    #   next start coord:
    #    ...
    
    # For all existing sections in h_, split depending on new start and end
    
    ### Check Start Pos ###
    ## EXAMPLE:
    ## h_     |-----------|       |---------|
    ## new   start->|----------------|<-end
    ## BECOMES:
    ## h_     |-----|-----|       |---------|
    ## new   start->|----------------|<-end
    
    for (h_key in hash::keys(h_)) {
        if (start > as.numeric(h_key) && start < h_[[h_key]][['End']]) {
            # add a new section starting at "start"
            h_[start] <- hash::hash(
                End = h_[[h_key]][['End']],
                Gain = h_[[h_key]][['Gain']],
                Loss = h_[[h_key]][['Loss']],
                "+Loss" = h_[[h_key]][['+Loss']]
            )
            # modify the old section to end at "start"
            h_[[h_key]][['End']] <- start
        
            break

        }

    }
    
    ### Check End Pos ###
    ## EXAMPLE:
    ## h_     |-----|-----|       |---------|
    ## new   start->|----------------|<-end
    ## BECOMES:
    ## h_     |-----|-----|       |--|------|
    ## new   start->|----------------|<-end
    
    for (h_key in hash::keys(h_)) {
        if (end > as.numeric(h_key) & end < h_[[h_key]][['End']]) {
            # add a new section starting at "end"
            h_[end] <- hash::hash(
                End = h_[[h_key]][['End']],
                Gain = h_[[h_key]][['Gain']],
                Loss = h_[[h_key]][['Loss']],
                "+Loss" = h_[[h_key]][['+Loss']]
            )
            # modify the old section to end at "end"
            h_[[h_key]][['End']] <- end
        
            break

        }

    }
    
    ### Now add new section by splitting it accordingly
    start_vals <- sort(as.numeric(hash::keys(h_)))
    prev_end <- start
    for (start_val in start_vals) {
        ## Check if existing section is completely within new region
        ## Update existing section accordingly by adding to gain/loss
        ## EXAMPLE: Adjust gain/loss of ** sections
        ## h_     |-----|*****|       |**|------|
        ## new   start->|----------------|<-end
      
        if (
            start_val >= start
            && h_[[as.character(start_val)]][['End']] <= end
        ) {
            # keep track of original section order
            h_[[as.character(start_val)]][[type]][[
                length(h_[[as.character(start_val)]][[type]]) + 1
            ]] <- end
        
            ## Add a new section from prev_end
            ## This fills in the gaps
            ## EXAMPLE:
            ## h_     |-----|-----|       |--|------|
            ## new   start->|----------------|<-end
            ## BECOMES:
            ## h_     |-----|-----|-------|--|------|
            ## new   start->|----------------|<-end
        
            if (prev_end < start_val) {
                h_[prev_end] <- hash::hash(
                    End = start_val,
                    Gain = list(),
                    Loss = list(),
                    "+Loss" = list()
                )
                h_[[as.character(prev_end)]][[type]] <- list(end)
          
            }
        
            # set new previous end
            prev_end <- h_[[as.character(start_val)]][['End']]

        }

    }
    
    ## Check if no new sections have been added, this indicates no overlaps
    ## Create a new section
    ## EXAMPLE:
    ## h_       |-----|             |----|
    ## new         start->|-------|<-end
    ## BECOMES:
    ## h_       |-----|   |-------| |----|
    ## new         start->|-------|<-end
    
    if (prev_end == start) {
        h_[start] <- hash::hash(
            End = end,
            Gain = list(),
            Loss = list(),
            "+Loss" = list()
        )
        h_[[as.character(start)]][[type]] <- list(end)
      
    } else {
        ## Check if any remaining end section is outside of existing sections
        ## Create a new end section accordingly
        ## EXAMPLE:
        ## h_       |-----| |----|
        ## new       start->|-------|<-end
        ## BECOMES:
        ## h_       |-----| |----|--|
        ## new       start->|-------|<-end
      
        if (prev_end < end) {
            h_[prev_end] <- hash::hash(
                End = end,
                Gain = list(),
                Loss = list(),
                "+Loss" = list()
            )
            h_[[as.character(prev_end)]][[type]] <- list(end)
        
        }

    }
    
    return(h_)

} # insertSection

#' Delete Intersecting Gain/Loss Regions
#' 
#' @description
#' This function removes overlapping gain and loss regions from the merge data structure,
#' implementing the biological logic that simultaneous gains and losses in the same
#' genomic region cancel each other out. It prioritizes +Loss events over regular Loss events.
#' 
#' @param h_ Hash data structure containing genomic intervals with gain/loss annotations
#' 
#' @return Modified hash data structure with overlapping gain/loss regions removed
#' 
#' @details
#' The function implements a two-step deletion process:
#' \enumerate{
#'   \item **+Loss cancellation**: Remove equal numbers of Gain and +Loss events first,
#'     as +Loss events take priority in cytogenetic interpretation
#'   \item **Standard cancellation**: Remove equal numbers of remaining Gain and Loss events
#' }
#' 
#' Deletion rules:
#' \itemize{
#'   \item If gains == (+Loss + Loss), delete the entire section (complete cancellation)
#'   \item Otherwise, remove pairs of overlapping events while preserving order
#'   \item Maintains original event order by removing from the beginning of lists
#' }
#' 
#' This ensures that complex karyotypes with multiple overlapping aberrations
#' are resolved according to cytogenetic conventions.
deleteIntersections <- function(h_) {

    start_vals <- sort(as.numeric(hash::keys(h_)))

    for (start_val in start_vals) {
        ## If both gain and loss are > 0, then there is overlap
        ## If min_val is > 1, then there is duplicate overlap
        ## Only delete one gain for each loss or one loss for each gain
        gain <- length(h_[[as.character(start_val)]][['Gain']])
        loss <- length(h_[[as.character(start_val)]][['Loss']])
        ploss <- length(h_[[as.character(start_val)]][['+Loss']])

        if (gain == (ploss + loss)) {
            # Delete section
            hash::delete(start_val, h_)

        } else {
            # Delete ploss first
            min_val <- min(gain, ploss)
        
            # Delete leading elements from Gain and Loss lists
            if (min_val > 0) {
                for (i in 1:min_val) {
                    h_[[as.character(start_val)]][['Gain']][[1]] <- NULL
                    h_[[as.character(start_val)]][['+Loss']][[1]] <- NULL

                }

            }
        
            ## If both gain and loss are > 0, then there is overlap
            ## If min_val is > 1, then there is duplicate overlap
            ## Only delete one gain for each loss or one loss for each gain
            gain <- length(h_[[as.character(start_val)]][['Gain']])
            loss <- length(h_[[as.character(start_val)]][['Loss']])

            # Delete loss afterwards
            min_val <- min(gain, loss)

            # Delete leading elements from Gain and pLoss lists
            if (min_val > 0) {
                for (i in 1:min_val) {
                    h_[[as.character(start_val)]][['Gain']][[1]] <- NULL
                    h_[[as.character(start_val)]][['Loss']][[1]] <- NULL

                }

            }

        }

    }

    return(h_)

}

#' Get Contiguous Section Extension
#' 
#' @description
#' This recursive function crawls the hash data structure to build contiguous sections
#' by extending genomic intervals that share the same aberration type and connect
#' at their boundaries. It ensures that adjacent intervals of the same type are
#' merged into single, longer intervals.
#' 
#' @param h__ Hash data structure containing genomic intervals organized by chromosome
#' @param section Named vector representing the current genomic section being extended,
#'   containing: Chr, Start, End, Type
#' @param orig_end Numeric value of the original end coordinate used for matching
#'   adjacent sections
#' 
#' @return Extended section with updated End coordinate if contiguous sections
#'   of the same type are found, otherwise the original section unchanged
#' 
#' @details
#' The function works by:
#' \enumerate{
#'   \item Checking if there's an adjacent section starting where current section ends
#'   \item Verifying the adjacent section has the same aberration type (Gain/Loss/+Loss)
#'   \item Extending the current section to encompass the adjacent section
#'   \item Recursively continuing to find further extensions
#'   \item Cleaning up empty sections after merging
#' }
#' 
#' This recursive approach ensures that all contiguous regions of the same type
#' are merged into single intervals, simplifying downstream analysis.
getContiguousSection <- function(h__, section, orig_end) {

    # This is a recursive function that crawls the hash to build
    # contiguous sections
      
    sect_end <- section[['End']]
    if (hash::has.key(as.character(sect_end), h__)) {
        # extend section only if it matches the original ending of the previous section
        # and also matches Type (e.g., Gain/Loss)
        if (orig_end %in% h__[[sect_end]][[section[['Type']]]]) {
            # delete item from list
            orig_end_index <- match(orig_end, h__[[sect_end]][[section[['Type']]]])
            h__[[sect_end]][[section[['Type']]]][[orig_end_index]] <- NULL
          
            # extend section
            section[['End']] <- h__[[sect_end]][['End']]
          
            if (
                length(h__[[sect_end]][['Gain']]) == 0
                && length(h__[[sect_end]][['Loss']]) == 0
                && length(h__[[sect_end]][['+Loss']]) == 0
            ) {
                # delete entire section if no more gain or loss
                hash::delete(sect_end, h__)

            }
          
            # continue searching recursively
            section <- getContiguousSection(h__, section, orig_end)
          
        }
        
    }

    return(section)

}
  
#' Merge Adjacent Genomic Sections
#' 
#' @description
#' This function processes the hash data structure to create a final table of merged
#' genomic intervals. It systematically processes each chromosome and merges adjacent
#' sections of the same aberration type into contiguous regions, producing the final
#' output format for CytoConverter results.
#' 
#' @param h_ Hash data structure containing genomic intervals organized by chromosome,
#'   with each chromosome containing sections with Gain, Loss, and +Loss annotations
#' 
#' @return Data frame with columns Chr, Start, End, Type containing the final merged
#'   genomic intervals ready for output
#' 
#' @details
#' The function implements a systematic merging process:
#' \enumerate{
#'   \item **Initialize output**: Creates empty data frame with standard column structure
#'   \item **Process by chromosome**: Iterates through each chromosome in the hash
#'   \item **Build contiguous sections**: For each starting position, identifies the
#'     aberration type and extends it using getContiguousSection()
#'   \item **Clean up data structure**: Removes processed sections from the hash
#'   \item **Accumulate results**: Adds completed intervals to the output table
#' }
#' 
#' Priority order for aberration types when multiple types exist at the same position:
#' Gain → Loss → +Loss
#' 
#' This ensures consistent output formatting and proper merging of overlapping regions.
mergeAdjacentSections <- function(h_) {

   
    coord_table <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(coord_table) <- c('Chr', 'Start', 'End', 'Type')
    
    for (chr in hash::keys(h_)) {
        # while there are still values in the section hash
        start_vals <- sort(as.numeric(hash::keys(h_[[chr]])))
        while (length(start_vals)) {
            # initialize the contiguous section
            start <- start_vals[1]
            ch_start <- as.character(start)
            end <- as.numeric(h_[[chr]][[ch_start]][['End']])
            orig_end <- 0
        
            type <- ''
            if (length(h_[[chr]][[ch_start]][['Gain']]) > 0) {
                type <- 'Gain'

            } else if (length(h_[[chr]][[ch_start]][['Loss']]) > 0) {
                type <- 'Loss'

            } else if (length(h_[[chr]][[ch_start]][['+Loss']]) > 0) {
                type <- '+Loss'

            } else {
                # This should never happen
                print('Error, should have deleted this key')

            }
        
            # get original end value and delete from list
            orig_end <- h_[[chr]][[ch_start]][[type]][[1]]
            h_[[chr]][[ch_start]][[type]][[1]] <- NULL
        
            # delete section if no more gain or ploss
            if (
                length(h_[[chr]][[ch_start]][['Gain']]) == 0
                && length(h_[[chr]][[ch_start]][['+Loss']]) == 0
                && length(h_[[chr]][[ch_start]][['Loss']]) == 0
            ) {
                hash::delete(ch_start, h_[[chr]])

            }
        
            # first part of contiguous section
            section <- c(
                Chr = chr,
                Start = start,
                End = end,
                Type = type
            )

            # build the contiguous section recursively
            section <- getContiguousSection(h_[[chr]], section, orig_end)
        
            # add contiguous section to the table
            # find existing sections that fit at the end of the new section
            end_index <- which(
                as.numeric(coord_table[['Start']]) == end
                & coord_table[['Type']] == type
            )

            # or existing sections that fit at the beginning of the new section
            start_index <- which(
                as.numeric(coord_table[['End']]) == start
                & coord_table[['Type']] == type
            )

            if (length(start_index) > 0) {
                if (length(end_index) > 0) {
                    # modify existing row to combine with new section and
                    # existing end section
                    coord_table[start_index[1], ][['End']] <- coord_table[end_index[1], ][['End']]

                    # delete existing end section
                    coord_table <- coord_table[-end_index[1], ]

                } else {
                    # modify existing row to combine with new section at the end
                    coord_table[start_index[1], ][['End']] <- section[['End']]

                }

            } else {
                if (length(end_index) > 0) {
                    # modify existing row to combine with new section at the beginning
                    coord_table[end_index[1], ][['Start']] <- section[['Start']]

                } else {
                    # add new section to table
                    coord_table <- rbind(coord_table, as.data.frame(t(section)))

                }

            }
        
            # update list of start values for loop's stopping condition
            start_vals <- sort(as.numeric(hash::keys(h_[[chr]])))
            
        }
    }

    return(coord_table)

}

#' Merge Gain/Loss Table Using Hash-Based Algorithm
#' 
#' @description
#' This is the main merging function that processes a table of genomic intervals
#' containing gains, losses, and other aberrations. It uses a hash-based algorithm
#' to efficiently merge overlapping intervals while properly handling gain/loss
#' cancellations and preserving non-gain/loss aberrations.
#' 
#' @param M Data frame containing genomic intervals with columns:
#'   \itemize{
#'     \item Chr - Chromosome identifier
#'     \item Start - Start coordinate (numeric)
#'     \item End - End coordinate (numeric)
#'     \item Type - Aberration type ("Gain", "Loss", "+Loss", or other)
#'   }
#' @param keep_extras Logical flag indicating whether to preserve non-gain/loss
#'   aberrations in the output (default: FALSE)
#' 
#' @return Data frame with merged genomic intervals, containing only net gains
#'   and losses after cancellation, plus any extra aberrations if keep_extras=TRUE
#' 
#' @details
#' The merging algorithm:
#' \enumerate{
#'   \item **Separate data**: Extract non-gain/loss aberrations for later inclusion
#'   \item **Initialize hash structure**: Create chromosome-specific hash maps
#'   \item **Populate intervals**: Insert all gain/loss intervals using insertSection()
#'   \item **Cancel overlaps**: Remove intersecting gain/loss pairs using deleteIntersections()
#'   \item **Merge adjacent**: Combine contiguous intervals of the same type
#'   \item **Combine results**: Add back extra aberrations if requested
#' }
#' 
#' This approach efficiently handles complex karyotypes with multiple overlapping
#' aberrations while maintaining biological accuracy.
mergeTable <- function(M, keep_extras = F) {

    # Store non gains and losses for intermediate steps
    M_temp <- M[grep("Gain|Loss", M[, 4], invert = T), ]
  
    # Create an empty hash map for each chromosome
    chrs <- unique(M[, 1])
  
    h <- hash::hash()
    for (chr in chrs) {
        h[chr] <- hash::hash()
    }
  
    # Populate the hash map with sections
    for (row in 1:nrow(M)) {
        if (
            M[row, 'Type'] == 'Gain'
            || M[row, 'Type'] == 'Loss'
            || M[row, 'Type'] == '+Loss'
        ) {
            h[M[row, 'Chr']] <- insertSection(
                h[[M[row, 'Chr']]],
                as.numeric(M[row, 'Start']),
                as.numeric(M[row, 'End']),
                M[row, 'Type']
            )

        }
        
    }
  
    # Delete intersections
    for (key in hash::keys(h)) {
        h[key] <- deleteIntersections(h[[key]])

    }
  
    # Merge adjacent sections and create final table
    final_coord_table <- mergeAdjacentSections(h)
  
    # If want other info
    if (keep_extras) {
        final_coord_table <- rbind(final_coord_table, M_temp)

    }
  
    return(final_coord_table)

}


# Section to cancel out loss gains
#   First code that takes two intervals (start, loss), where gain is first, that
#   returns a 3-column data frame with the intersection missing.
#   Can have 0, 1, or two rowsAlso has a column for type.

# This function takes a table M (a data frame) with columns:
#   chromosome, start, stop, and type (Gain or Loss)
#   and does proper merging

# First we need a function that takes two intervals (start, loss) that
# returns NA if they don't overlap and merges them if they do.
# v1 and v2 are each two-vectors giving the interval

#' Merge Two Overlapping Genomic Intervals
#' 
#' @description
#' This function handles the merging of two genomic intervals, typically representing
#' opposing aberrations (gain vs loss). It calculates the non-overlapping portions
#' after cancellation and returns the remaining segments along with information
#' about whether the original intervals were modified.
#' 
#' @param v1 Numeric vector of length 3: c(start, end, type) for first interval
#' @param v2 Numeric vector of length 3: c(start, end, type) for second interval
#' 
#' @return List containing:
#'   \item{intervals}{Data frame with remaining non-overlapping intervals after merging}
#'   \item{modified}{Logical indicating whether intervals were changed (FALSE if no overlap)}
#' 
#' @details
#' The function handles various overlap scenarios:
#' \itemize{
#'   \item **Complete containment**: One interval completely contains another
#'   \item **Partial overlap**: Intervals overlap at their boundaries  
#'   \item **No overlap**: Intervals are completely separate
#'   \item **Exact match**: Intervals have identical coordinates (complete cancellation)
#' }
#' 
#' Order independence: The function automatically orders intervals by start position
#' to ensure consistent results regardless of input order.
#' 
#' Used primarily for gain/loss cancellation in cytogenetic analysis.
mergeDel <- function(v1, v2) {

    labs <- c(v1[3], v2[3])
    v1 <- as.vector(as.integer(as.numeric(as.character(v1[1:2]))))
    v2 <- as.vector(as.integer(as.numeric(as.character(v2[1:2]))))
    out <- as.data.frame(matrix(nrow = 0, ncol = 3))
    colnames(out) <- c("Start", "End", "Type")
    out <- cbind(
        as.integer(as.numeric(out[, 1])),
        as.integer(as.numeric(out[, 2])),
        out[, 3]
    )

    if (v1[1] > v2[1]) {
        labs <- labs[2:1]
        tmp <- v1
        v1 <- v2
        v2 <- tmp

    }
  
    if (v1[2] >= v2[2]) {
        out <- data.frame(
            Start = c(v1[1], v2[2]),
            End = c(v2[1], v1[2]),
            Type = rep(labs[1], 2)
        )

    } else {
        if (v1[2] >= v2[1]) {
            out <- data.frame(
                Start = c(v1[1], v1[2]),
                End = c(v2[1], v2[2]),
                Type = labs
            )

        } else {
            out <- data.frame(
                Start = c(v1[1], v2[1]),
                End = c(v1[2], v2[2]),
                Type = labs
            )

        }

    }
  
    bads <- which(out[, 1] >= out[, 2])
    if (length(bads) > 0) {
        out <- out[-bads, ]

    }

    origMat <- data.frame(
        Start = c(v1[1], v2[1]),
        End = c(v1[2], v2[2]),
        Type = labs
    )

    if (nrow(out) == nrow(origMat) && identical(out, origMat)) {
        out = list(out, TRUE)

    } else {
        out = list(out, FALSE)

    }

  return(out)

}

#' Merge Gain and Loss Matrices with Iterative Cancellation
#' 
#' @description
#' This function performs comprehensive merging of gain and loss intervals using
#' an iterative algorithm. It systematically compares each gain against all losses
#' and vice versa, removing overlapping regions to produce final non-conflicting
#' interval sets. This is an alternative to the hash-based mergeTable approach.
#' 
#' @param G Data frame containing gain intervals with columns: Start, End, Type
#' @param L Data frame containing loss intervals with columns: Start, End, Type
#' 
#' @return List containing:
#'   \item{gains}{Data frame with remaining gain intervals after cancellation}
#'   \item{losses}{Data frame with remaining loss intervals after cancellation}
#' 
#' @details
#' The iterative merging process:
#' \enumerate{
#'   \item **Process gains against losses**: For each gain interval, check against
#'     all loss intervals for overlaps, using mergeDel() to handle cancellations
#'   \item **Process losses against original gains**: Use the original gain set
#'     to process remaining losses, ensuring bidirectional cancellation
#'   \item **Handle special cases**: Processes intervals with "del" and "add" prefixes
#'     differently to preserve structural aberration information
#'   \item **Iterate until convergence**: Continue until no more modifications occur
#' }
#' 
#' This approach ensures complete cancellation of overlapping gain/loss pairs
#' while preserving the biological interpretation of complex karyotypes.
#' 
#' @note This function provides an alternative algorithm to mergeTable() for
#' scenarios requiring more explicit control over the merging process.
mergeDelmat <- function(G, L) {

    i <- 1
    ##G[,1:2]<-apply(G[,1:2],2,as.numeric)
    ##L[,1:2]<-apply(L[,1:2],2,as.numeric)
    origL <- L

    while (nrow(L) > 0 & i <= nrow(G)) {
        newL <- matrix(ncol = 3, nrow = 0)
        j <- 1
        modified <- TRUE
        #print(list(i, modified))
        while (j <= nrow(L) & modified) {
            nxt <- mergeDel(G[i, ], L[j, ])
            modified <- nxt[[2]]
            nxt <- nxt[[1]]
            w <- union(
                intersect(grep("del", nxt[, 3]), grep("^del", nxt[, 3], invert = T)),
                intersect(grep("add", nxt[, 3]), grep("^add", nxt[, 3], invert = T))
            )

            if (length(w) > 0) {
                newL <- rbind(
                    newL,
                    nxt[w, ],
                    if (nrow(L) > 1 && !is.vector(L[-j, ]) && nrow(L[-j, ]) >= j && !modified) {
                        L[-j, ][j:nrow(L[-j, ]), ]
                    }
                )

            } else if (nrow(nxt) > 0) {
                newL <- rbind(
                    newL,
                    if (nrow(L) > 1 && !is.vector(L[-j, ]) && nrow(L[-j, ]) >= j && !modified) {
                        L[-j, ][j:nrow(L[-j, ]), ]
                    }
                )

            } else if (!is.vector(L[-j, ])) {
                newL <- rbind(newL, apply(L[-j, ], 2, as.character))

            } else {
                newL <- rbind(newL, sapply(L[-j, ], as.character))

            }
            #print(list(j, modified, newL))
            j <- j + 1
        }
        L <- newL
        i <- i + 1
    }
  
    i <- 1
    while (nrow(G) > 0 & i <= nrow(origL)) {
        newG <- matrix(ncol = 3, nrow = 0)
        j <- 1
        modified <- TRUE
        #print(list(i, G))
    
        while (j <= nrow(G) & modified) {
            nxt <- mergeDel(G[j, ], origL[i, ])
            modified <- nxt[[2]]
            nxt <- nxt[[1]]
            w <- setdiff(
                grep(
                    "((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",
                    nxt[, 3]
                ),
                union(
                    intersect(
                        grep("del", nxt[, 3]), grep("^del", nxt[, 3], invert = T)
                    ),
                    intersect(
                        grep("add", nxt[, 3]), grep("^add", nxt[, 3], invert = T)
                    )
                )
            )

            if (length(w) > 0) {
                newG <- rbind(
                    newG,
                    nxt[w, ],
                    if (nrow(G) > 1 && !is.vector(G[-j, ]) && nrow(G[-j, ]) >= j && !modified) {
                        G[-j, ][j:nrow(G[-j, ]), ]
                    }
                )

            } else if (nrow(nxt) > 0) {
                newG <- rbind(
                    newG,
                    if (nrow(G) > 1 && !is.vector(G[-j, ]) && nrow(G[-j, ]) >= j && !modified) {
                        G[-j, ][j:nrow(G[-j, ]), ]
                    }
                )

            } else if (!is.vector(G[-j, ])) {
                newG <- rbind(
                    newG,
                    apply(G[-j, ], 2, as.character)
                )

            } else {
                newG <- rbind(
                    newG,
                    sapply(G[-j, ], as.character)
                )

            }
      
            # print(list(j, modified, newG))
            j <- j + 1

        } # while
        G <- newG
        i <- i + 1

    } # while
  
    G[, 1] <- as.numeric(as.character(G[, 1]))
    G[, 2] <- as.numeric(as.character(G[, 2]))
    L[, 1] <- as.numeric(as.character(L[, 1]))
    L[, 2] <- as.numeric(as.character(L[, 2]))
  
    out <- data.frame(
        Start = c(G[, 1], L[, 1]),
        End = c(G[, 2], L[, 2]),
        Type = c(G[, 3], L[, 3])
    )

    return(out)

}

#' Big Deletion Merge for Multiple Chromosomes
#' 
#' @description
#' This function applies the mergeDelmat algorithm across multiple chromosomes
#' in a single operation. It processes each chromosome separately while handling
#' special structural aberrations like translocations, isodicentric, Robertsonian,
#' tricentric, and dicentric chromosomes.
#' 
#' @param M Data frame containing genomic intervals with columns:
#'   Chr, Start, End, Type (where Type may include gains, losses, and structural aberrations)
#' 
#' @return Data frame with merged intervals where gain/loss cancellations have been
#'   applied chromosome by chromosome
#' 
#' @details
#' The function handles complex structural aberrations by:
#' \itemize{
#'   \item Identifying structural aberrations using regex patterns for t(, idic(, rob(, trc(, dic(
#'   \item Separating gains/structural aberrations from deletions/additions
#'   \item Applying mergeDelmat to opposing interval types within each chromosome
#'   \item Preserving structural aberration information during merging
#' }
#' 
#' Used primarily for complex karyotypes containing both simple gains/losses
#' and structural rearrangements.
bigDelMerge <- function(M) {

    ##ask what does this do
    out <- M[-(1:nrow(M)), ]
  
    chrs <- unique(M[, 1])
  
    for (i in 1:length(chrs)) {
        w <- which(M[, 1] == chrs[i])
        Msub <- M[w, ]
    
        startDel <- grep(
            "((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\())",
            Msub[, 4]
        )[1]
        
        wl <- union(
            intersect(
                grep("del", Msub[, 4]),
                grep("^del", Msub[, 4], invert = T)
            ),
            intersect(
                grep("add", Msub[, 4]),
                grep("^add", Msub[, 4], invert = T)
            )
        )
    
    if (!is.na(startDel)) {
        startDel <- 1:startDel
        wl <- setdiff(
            union(
                intersect(
                    grep("del", Msub[, 4]),
                    grep("^del", Msub[, 4], invert = T)
                ),
                intersect(
                    grep("add", Msub[, 4]),
                    grep("^add", Msub[, 4], invert = T)
                )
            ),
            startDel
        )

    }

    wg <- setdiff(
        grep(
            "((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",
            Msub[, 4]
        ),
        wl
    )
    
    if (length(wg) == 0 | length(wl) == 0) {
        out <- rbind(out, Msub)

    } else {
        nxt <- mergeDelmat(Msub[wg, -1], Msub[wl, -1])
        if (nrow(nxt) > 0) {
            nxt <- cbind(chrs[i], nxt)
            colnames(nxt)[1] <- "Chr"
        }
        out <- rbind(out, nxt)

    }

  }
  
  return(out)

}

#' Merge Deletions with Complex Structural Aberrations
#' 
#' @description
#' This function handles the merging of deletions and additions in the context of
#' complex structural aberrations. It categorizes different types of aberrations
#' and applies appropriate merging strategies while preserving the biological
#' significance of structural rearrangements.
#' 
#' @param M Data frame containing genomic intervals with structural aberrations
#' @param Mainchr Character vector indicating the main chromosomes involved in the analysis
#' 
#' @return Data frame with processed genomic intervals where deletions have been
#'   appropriately merged with structural aberrations
#' 
#' @details
#' The function implements a sophisticated categorization system:
#' \enumerate{
#'   \item **Structural aberrations**: Translocations, dicentric, Robertsonian, etc.
#'   \item **Deletions**: Intervals marked with "del" but not starting with "del"
#'   \item **Additions**: Intervals marked with "add" but not starting with "add"
#'   \item **Gains**: Standard gain intervals
#' }
#' 
#' The merging process:
#' \itemize{
#'   \item First applies bigDelMerge for complex structural interactions
#'   \item Then performs final coordinate adjustments
#'   \item Preserves non-conflicting aberrations in the final output
#' }
#' 
#' This function is particularly important for complex constitutional and somatic
#' karyotypes involving multiple types of chromosomal rearrangements.
mergeDeletions <- function(M, Mainchr) {

    OldM <- M
    startDel <- grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))", M[, 4])[1]
  
    wdel <- union(
        intersect(grep("del", M[, 4]), grep("^del", M[, 4], invert = T)),
        intersect(grep("add", M[, 4]), grep("^add", M[, 4], invert = T))
    )
  
    if (!is.na(startDel)) {
        startDel <- 1:startDel
        wdel <- setdiff(
            union(
                intersect(grep("del", M[, 4]), grep("^del", M[, 4], invert = T)),
                intersect(grep("add", M[, 4]), grep("^add", M[, 4], invert = T))
            ),
            startDel
        )
    
    }
  
    resttosee <- setdiff(
        grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))", M[, 4]),
        wdel
    )
    rest <- (1:nrow(M))[-1 * c(wdel, resttosee)]
  
    M <- M[c(wdel, resttosee), ]
    if (length(resttosee) > 0 & length(wdel) > 0) {
        M <- bigDelMerge(M)

    }
  
    startDel <- grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))", M[, 4])[1]
  
    wdel <- union(
        intersect(grep("del", M[, 4]), grep("^del", M[, 4], invert = T)),
        intersect(grep("add", M[, 4]), grep("^add", M[, 4], invert = T))
    )
  
    if (!is.na(startDel)) {
        startDel <- 1:startDel
        wdel <- setdiff(
            union(
                intersect(grep("del", M[, 4]), grep("^del", M[, 4], invert = T)),
                intersect(grep("add", M[, 4]), grep("^add", M[, 4], invert = T))
            ),
            startDel
        )
    
    }

    wg <- setdiff(
        grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))", M[, 4]),
        wdel
    )
    Mg <- M[wg, 1:3]
    Ml <- M[wdel, 1:3]
  
    if (length(wg) > 1) {
        M[wg, 1:3] <- Mg

    }
  
    if (length(wdel) > 1) {
        M[wdel, 1:3] <- Ml

    }
  
    # deletions that are to be removed
  
    w <- which(!is.na(M[, 2]))
    if (length(rest) > 0) {
        out <- rbind(M[w, ], OldM[rest, ])
        w <- c(1:2)

    } else {
        out <- M[w, ]

    }
  
    if (length(w) == 1) {
        out <- t(out)

    }
  
    return(out)

}




