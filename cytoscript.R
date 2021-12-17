#' CytoConverter Main Function
#'
#' @description
#' This function accepts string or matrix input. Strings must be karyotypes and tables must have
#' sample name in the first column and karyotype in the second column. The function outputs a 
#' table. Assumes order of abberations is from top to bottom, assumes all output goes top to bottom
#' (eg no q20p23).
#' 
#' @param in_data
#' @param build
#' @param constitutional
#' @param guess
#' @param guess_q
#' @param forMtn
#' @param orOption
#' @param sexstimate
#' @param allow_Shorthand

CytoConverter <- function(
        in_data,
        build = "GRCh38",
        constitutional = T,
        guess = F,
        guess_q = F,
        forMtn = T,
        orOption = T,
        sexstimate = F,
        allow_Shorthand = F
) {
    
    # Load the cyto reference table, cyto_ref_table

    if (build == "GRCh38") {
        cyto_ref_table <- sapply(
            as.data.frame(read.delim(
                "Builds/cytoBand_GRCh38.txt", header = FALSE
            )),
            as.character
        )
      
    } else if (build == "hg19") {
        cyto_ref_table <- sapply(
            as.data.frame(read.delim(
                "Builds/cytoBand_hg19.txt", header = FALSE
            )),
            as.character
        )
      
    } else if (build == "hg18") {
        cyto_ref_table <- sapply(
            as.data.frame(read.delim(
                "Builds/cytoBand_hg18.txt", header = FALSE
            )),
            as.character
        )
      
    } else if (build == "hg17") {
        cyto_ref_table <- sapply(
            as.data.frame(read.delim(
                "Builds/cytoBand_hg17.txt", header = FALSE
            )),
            as.character
        )
      
    } else if (is.null(build)) {
        cyto_ref_table <- sapply(
            as.data.frame(read.delim(
                "Builds/cytoBand_GRCh38.txt", header = FALSE
            )),
            as.character
        )
      
    } else {
        return("Error : build incorrectly specified")

    }
   
    # ref_table stores the end coordinate for each chromosome

    ref_table <- as.data.frame(
        cyto_ref_table[
            sapply(
                unique(cyto_ref_table[, 1]),
                function(x) {
                    grep(x, cyto_ref_table[, 1])[
                        length(grep(paste(x, "$", sep = ""), cyto_ref_table[, 1]))
                    ]
                }
            ), 
        ][, c(1, 3)]
    )
    ref_table <- apply(ref_table, 2, as.character)

    # Table with desired output
    Final_table <- matrix(ncol = 5, nrow = 0)
    colnames(Final_table) <- c("Sample ID", "Chr", "Start", "End", "Type")
    
    # Convert any single string into a table
    if (is.vector(in_data)) {
        in_data <- t(matrix(c("sample", in_data)))
    }
    
    # Dump table of stuff containing unprocessed reads
    # Write this later, get all fish recorded
    Dump_table <- matrix(ncol = 3, nrow = 0)

    # Double check that this does not delete later data potentially
    fish_table <- in_data[grep("ish.*$", in_data[, 2]), ]
    if (is.vector(fish_table) && length(fish_table) > 0) {
        Dump_table <- rbind(Dump_table, c(fish_table, "Warning in fish reading"))
        # Now take out ish readings
        in_data[, 2] <- gsub("ish.*$", "", as.character(in_data[, 2]))
      
    } else {
        Dump_table <- rbind(Dump_table, cbind(fish_table, "Warning in fish reading"))
        # Now take out ish readings
        in_data[, 2] <- gsub("ish.*$", "", as.character(in_data[, 2]))

    }
    
    # Set everything to lowercase, then set x and y to uppercase
    in_data[, 2] <- tolower(in_data[, 2])
    in_data[, 2] <- chartr("x", "X", in_data[, 2])
    in_data[, 2] <- chartr("y", "Y", in_data[, 2])
    
    Con_data = matrix(nrow = 0, ncol = 2)

    # Split cell lines
    for (i in 1:nrow(in_data)) {
        # If activated, allow for shorthand outside of clones
        # and outside of translocations. E.g., del(1), or if its a clone

        if (
            !is.na(in_data[i, 2])
            && nchar(as.character(in_data[i, 2])) > 0
            && (allow_Shorthand == T || grepl("/", in_data[i, 2]))
        ) {

            # Split string
            data_split <- unlist(strsplit(in_data[i, 2], ","))
            data_split <- unlist(strsplit(
                gsub("/", "/-;-;opj", data_split), "/"
            ))
            mut_index <- grep("\\(.*\\)(?!\\()", data_split, perl = T)

            if (length(mut_index) > 0) {
                mut_list <- sapply(
                    data_split[mut_index],
                    function(x) {
                        gsub("\\[|\\]", "",
                            gsub("\\)", "\\\\)",
                                gsub("\\(", "\\\\(",
                                    gsub("\\?", "\\\\?",
                                        gsub("\\+", "",
                                            substr(x, 2, nchar(x))
                                        )
                                    )
                                )
                            )
                        )
                    }
                )

                index_match <- sapply(
                    mut_list,
                    function(x) {
                        grep(x, data_split)[1]
                    }
                )

                index_match <- index_match[
                    intersect(
                        which(!is.na(index_match)),
                        which(index_match > 0)
                    )
                ]

                if (
                    is.numeric(mut_index)
                    && is.numeric(index_match)
                    && length(index_match) > 0
                ) {
                    additional_to_add <- substring(
                        names(index_match), 0, regexpr("\\+", names(index_match))
                    )
                    data_split[mut_index] <- paste(
                        additional_to_add , data_split[index_match], sep = ""
                    )

                }

            }
        
            in_data[i, 2] <- gsub(
                ",-;-;opj",
                "/",
                paste(data_split, collapse = ",")
            )
        
        } # if
      
        if (
            !is.na(in_data[i, 2])
            & nchar(as.character(in_data[i, 2])) > 0
        ) {
            c2 <- strsplit(in_data[i, 2], split = "/")[[1]]
            c1 <- paste(in_data[i, 1], 1:length(c2), sep = "_")
            Con_data <- rbind(Con_data, cbind(c1, c2))

        } else {
            Con_data <- rbind(
                Con_data, cbind(paste(in_data[i, 1], "_1"), in_data[i, 2])
            )

        }
      
        # Make idems here
        if (any(grepl("idem|sl|sdl", Con_data[, 2]))) {
            idem_index <- grep("idem|sl|sdl", Con_data[, 2])
            ##### !! CHECK THIS !! #####
            temp_data <- Con_data[(idem_index[1] - 1):idem_index[length(idem_index)], ]

            if (!is.vector(temp_data) && nrow(temp_data) > 1) {

                for (j in 2:nrow(temp_data)) {
                    # Make sure temp data has more than 1 row
                    if (grepl("idem|sl", temp_data[j, 2])) {
                        prev <- unlist(strsplit(temp_data[1, 2], "\\["))[1]

                    } else {
                        # Implement this so it can handel two sl1 in sucession and sdl1 sdl2
                        prev <- unlist(strsplit(temp_data[j - 1, 2], "\\["))[1]

                    }
                    prev <- unlist(strsplit(prev, ","))
                    prev <- prev[2:length(prev)]
                    sexchromprev <- prev[grep("^[XY]+", prev)[1]]
                    autochromprev <- prev[grep("^[XY]+", prev, invert = T)]
            
                    clonecount = 1
                    if (grepl("idemx|idemX|slx|slX|sdl*x|sdl*X", temp_data[j, 2])) {
                        # Check if there is a X3 etc value, if there is pick that up, store,
                        # add to addtot, make it process through twice later
              
                        clonecount = unlist(strsplit(temp_data[j, 2], "idemx|idemX|slx|slX"))[2]
                        if (grepl("-|~", clonecount)) {
                            clonecount <- unlist(strsplit(clonecount, "~|-"))[1]
                        }
                        clonecount = as.numeric(
                            unlist(strsplit(as.character(clonecount), "\\[|,"))[1]
                        )
                    }

                    cur <- unlist(strsplit(temp_data[j, 2], ","))
                    if (!is.na(sexchromprev)) {
                        if (grepl("^[XY]+", cur[2])) {
                            cur[2] <- paste(
                                paste(
                                    rep(sexchromprev, clonecount),
                                    collapse = ''
                                ),
                                cur[2],
                                sep = ''
                            )

                        } else if (grepl("idem|sl|sdl", cur[2])) {
                            cur[2] <- paste(rep(sexchromprev, clonecount), collapse = '')

                        } else {
                            cur <- c(
                                cur[1],
                                paste(rep(sexchromprev, clonecount), collapse = ''),
                                cur[2:length(cur)]
                            )
                
                        }

                    }

                    if (!is.null(autochromprev)) {
                        if (length(cur) > 2) {
                            cur <- c(
                                cur[1:2],
                                rep(autochromprev, clonecount),
                                cur[3:length(cur)]
                            )

                        } else {
                            cur <- c(cur[1:2], rep(autochromprev, clonecount))

                        }

                    }

                    cur <- cur[grep("idem|sl|sdl", cur, invert = T)]
                    temp_data[j, 2] <- paste(c(cur, "ids"), sep = '', collapse = ',')

                } # for (j in 2:nrow(temp_data)) {

                Con_data[
                    (idem_index[1] - 1):idem_index[length(idem_index)],
                ] <- temp_data

            } else {
                # Put it in the error table
                Dump_table <- rbind(
                    Dump_table,
                    c(temp_data, "no proceeding cell line to refer to")
                )

            }

        } # Make idems here
    
    } # for (i in 1:nrow(in_data)) { # Split cell lines

    rownames(Con_data) <- 1:nrow(Con_data)
    
    Con_data[, 2] <- gsub(" ", "", Con_data[, 2])
    Con_data[, 2] <- gsub("crYp", "", Con_data[, 2])
    
    # Taking all reads with ? and ~ as well into dump table
    # Adding marker and add material as well
    unsure_table <- Con_data[grep("\\?|\\~|inc|mar|add", Con_data[, 2]), ]
    Dump_table <- if (is.vector(unsure_table)) {
        rbind(
            Dump_table,
            c(
                unsure_table,
                "Warning in ?,~,marker, unknown additional material or incomplete karyotype detected"
            )
        )

    } else {
        rbind(
            Dump_table,
            cbind(
                unsure_table,
                "Warning in ?,~, marker, unknown additional material or incomplete karyotype detected"
            )
        )

    }
    
    # Main program qdx qd
    # Rewrite for new data per row
    transloctable <- data.frame()
    
    for (i in 1:nrow(Con_data)) {
        Cyto_sample <- unlist(strsplit(as.character(Con_data[i, 2]), split = ",|\\["))

        # If X or Y in first col, move it
        if (grepl("X|Y", Cyto_sample[1])) {
            Cyto_sample[1] <- gsub(
                substr(
                    Cyto_sample[1],
                    regexec("X|Y", Cyto_sample[1])[[1]][1],
                    nchar(Cyto_sample[1])
                ),
                paste(
                    ",",
                    substr(
                        Cyto_sample[1],
                        regexec("X|Y", Cyto_sample[1])[[1]][1],
                        nchar(Cyto_sample[1])
                    ),
                    sep = ""
                ),
                Cyto_sample[1]
            )
            Cyto_sample <- unlist(
                strsplit(as.character(Cyto_sample), split = ",|\\[")
            )

        }

        # Take away extra accidental commas
        if (length(which(str_length(Cyto_sample) == 0)) > 0) {
            Cyto_sample <- Cyto_sample[-1 * which(str_length(Cyto_sample) == 0)]
        }
      
        # If idem or sl, cancel out any - details, add in any shorthands between code
        if (any(grepl("ids|idem|sl|sd", Cyto_sample)) || allow_Shorthand) {
            # index of anything in a clonal evolution step with no )(
            # think about what to do if assymmetric (2 defined, one mystery, one defined, 2 mystery)
            # takes first occurance if more than 2 of the same category -- f
            # fix it to skip
            mut_index <- grep("\\(.*\\)(?!\\()", Cyto_sample, perl = T)
            if (length(mut_index) > 0) {
                mut_list <- sapply(
                    Cyto_sample[mut_index],
                    function(x) {
                        gsub("\\)", "\\\\)",
                            gsub("\\(", "\\\\(",
                                gsub("\\?", "\\\\?",
                                    gsub("\\+", "",
                                        substr(x, 2, nchar(x))
                                    )
                                )
                            )
                        )
                    }
                )

                index_match <- sapply(
                    mut_list,
                    function(x) {
                        grep(x, Cyto_sample)[1]
                    }
                )

                index_match <- index_match[
                    intersect(which(!is.na(index_match)), which(index_match > 0))
                ]

                if (
                    is.numeric(mut_index)
                    && is.numeric(index_match)
                    && length(index_match) > 0
                ) {
                    additional_to_add <- substring(
                        names(index_match),
                        0,
                        regexpr("\\+", names(index_match))
                    )

                    Cyto_sample[mut_index] <- paste(
                        additional_to_add , Cyto_sample[index_match], sep = ""
                    )

                }

            }

            # Index of anything thats - in a clonal evolution step
            mut_gone_index <- grep("-[[:alpha:]]+\\(", Cyto_sample)
            if (length(mut_gone_index) > 0) {
                mut_list <- sapply(
                    Cyto_sample[mut_gone_index],
                    function(x) {
                        gsub("\\)", "\\\\)",
                            gsub("\\(", "\\\\(",
                                gsub("\\?", "\\\\?",
                                    substr(x, 2, nchar(x))
                                )
                            )
                        )
                    }
                )

                index_cancel <- sapply(
                    mut_list,
                    function(x) {
                        grep(x, Cyto_sample)[1]
                    }
                )

                index_cancel <- index_cancel[
                    intersect(which(!is.na(index_cancel)), which(index_cancel > 0))
                ]

                if (
                    is.numeric(mut_gone_index)
                    && is.numeric(index_cancel)
                    && length(index_cancel) > 0
                ) {
                    Cyto_sample <- Cyto_sample[-1 * c(mut_gone_index, index_cancel)]

                }

            }

        } # if (any(grepl("ids|idem|sl|sd", Cyto_sample)) || allow_Shorthand) {
      
        # Check to make sure the input is roughly a karyotype before processing
        if (
            grepl(
                "[[:digit:]]+((~|-)[[:digit:]]+)*(<[[:digit:]]+n>)*,",
                paste(Cyto_sample, collapse = ',', sep = '')
            )
            && grepl("[[:digit:]]", Cyto_sample[1])
        ) {
            tottable <- tryCatch({
                mod_rowparser$rowparse(
                    cyto_ref_table,
                    ref_table,
                    Cyto_sample,
                    Con_data,
                    i,
                    guess,
                    guess_q,
                    orOption,
                    constitutional,
                    transloctable,
                    Dump_table,
                    sexstimate,
                    forMtn
                )
            }, error = function(e) {
                return(gsub("\n", " ", paste(e, "in", i, "sample")))
            }, finally = {
                print(paste("Parsed sample: ", Con_data[i, 1], Con_data[i, 2]))
            })

            if (is.character(tottable) & length(tottable) == 1) {
                Dump_table <- rbind(Dump_table, c(Con_data[i,], tottable))
                transloctable <- data.frame()

            } else if (!is.list(tottable)) {
                Dump_table <- tottable
                transloctable <- data.frame()

            } else {
                sorted_sample_table <- tottable[[1]]
                Dump_table <- tottable[[2]]

                # Master table with names
                temp2 <- cbind(
                    rep(as.character(Con_data[i, 1]), nrow(sorted_sample_table)),
                    sorted_sample_table
                )
          
                # Sort temp2 to remove gains and losses that are cooresponding
                    
                colnames(temp2) <- colnames(Final_table)
                Final_table <- rbind(Final_table, temp2)

                # Correct data formats from factors into chr int int chr
                if (nrow(Final_table) == 1) {
                    Final_table <- as.data.frame(Final_table)
                    Final_table[, 1] <- as.character(Final_table[, 1])
                    Final_table[, 2] <- as.character(Final_table[, 2])
                    Final_table[, 3] <- as.integer(as.numeric(as.character(Final_table[, 3])))
                    Final_table[, 4] <- as.integer(as.numeric(as.character(Final_table[, 4])))
                    Final_table[, 5] <- as.character(Final_table[, 5])
                }
          
                # Make list to store translocations and insertions
                # If clone line, then, carry over transloctable to next one
                if (i < nrow(Con_data)) {
                    curname <- strsplit(Con_data[i, 1], "_")[[1]]
                    nextname <- strsplit(Con_data[i + 1, 1], "_")[[1]]
                    if (
                        curname[-length(curname)] == nextname[-length(nextname)]
                        & as.numeric(curname[length(curname)]) + 1
                            == as.numeric(nextname[length(nextname)])
                    ) {
                        # Keep previous transloctable if sample names the same and cell line
                        # is in correct order
                        transloctable <- tottable[[3]]

                    } else {
                        transloctable <- data.frame()

                    }

                }

            }

        } else {
            if (grepl("[[:digit:]]", Cyto_sample[1])) {
                Dump_table <- rbind(
                    Dump_table,
                    c(Con_data[i,], "Warning in karyotype number not specified")
                )
          
            } else {
                Dump_table <- rbind(
                    Dump_table,
                    c(Con_data[i,], "Warning in karyotype is incorrect")
                )

            }

        }

    } # for (i in 1:nrow(Con_data)) {
    
    # Remove duplicates and redundancy (same region loss/gain within group do this first,
    # then duplicates)
    
    # counting function for % mutated
    out <- Con_data
    
    col3 <- rep("unknown", nrow(out))
    
    IDs <- as.character(out[, 1])
    
    IDs <- lapply(
        strsplit(IDs, split = "_"),
        function(x) {
            x[-length(x)]
        }
    )
    
    IDs <- unlist(
        lapply(
            IDs,
            function(x) {
                paste(x, collapse = "_")
            }
        )
    )

    IDs <- unique(IDs)
    
    for (i in 1:length(IDs)) {
        # Get all ids with same base name
        # gsub away all regular expressions
        w <- grep(
            gsub("\\+", "\\\\+",
                gsub("\\*", "\\\\*",
                    gsub("\\!", "\\\\!",
                        paste("^", IDs[i], "_", sep = "")
                    )
                )
            ),
            out[, 1]
        )

        nxt <- unlist(
            lapply(
                strsplit(
                    as.character(out[w, 2]),
                    split = "\\["
                ),
                function(x) {
                    x[2]
                }
            )
        )
      
        # Deal with cp more gracefully
        is_cp <- grepl("cp", nxt)
        nxt <- as.numeric(gsub("\\](,ids)*|cp", "", nxt))

        # If its cp just skip it--get more sophisticated processing later, that one unknonwn
        # other stuff can be calculated
        # maybe modifify later to have total, something sometihng etc
        # can make more complicated, so that cp stuff is unknown but added up to second sample,
        # it creates a percentage for that

        if (any(is.na(nxt))) {
            col3[w] <- "unknown"

        } else {
            if (any(is_cp)) {
                col3[w] <- paste("1-", nxt, " of ", sum(nxt), sep = '')

            } else {
                col3[w] <- paste(nxt, "of", sum(nxt))

            }

        }

    }

    # Something is going wrong when its not just a normal cell and theres no more
    # (says unknown instead of 1)
    # turn ids back into over all list
    IDs <- as.character(Con_data[, 1])
    
    # Compare to final table
    if (is.vector(Final_table)) {
        colvect <- vector(length = 1)
      
    } else {
        colvect <- vector(length = length(Final_table[, 1]))

    }
    
    for (i in 1:length(IDs)) {
        if (is.vector(Final_table)) {
            colvect[grep(IDs[i], Final_table[1])] <- col3[i]
        
        } else {
            colvect[grep(IDs[i], Final_table[, 1])] <- col3[i]

        }

    }
    
    Final_Final_table <- cbind(Final_table, colvect)
    colnames(Final_Final_table) <- c(
        "Sample ID", "Chr", "Start", "End", "Type", "Cells Present"
    )

    if (is.vector(Final_Final_table)) {
        Final_Final_table <- Final_Final_table[grep("Loss|Gain", Final_Final_table[5]),]
      
    } else {
      Final_Final_table <- Final_Final_table[grep("Loss|Gain", Final_Final_table[, 5]),]

    }

    Final_Final_table <- na.omit(Final_Final_table)
    
    if (is.vector(Final_Final_table)) {
        Final_Final_table <- t(Final_Final_table)
        colnames(Final_Final_table) <- c(
            "Sample ID", "Chr", "Start", "End", "Type", "Percent Present"
        )

    }

    colnames(Dump_table) <- c("Sample ID", "Karyotype", "Error Message")
    
    result <- list(Final_Final_table, Dump_table)
    names(result) <- list("Result", "Error_log")
    
    return(result)

}




