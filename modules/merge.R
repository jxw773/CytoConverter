
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




