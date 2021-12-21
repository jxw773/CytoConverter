#' CytoConverter Row Parser
#' 
#' @description
#' This function parses each row of a karyotype table.

mod_utils <- modules::use('modules/utils.R')
mod_merge <- modules::use('modules/merge.R')
mod_colparser <- modules::use('modules/colparser.R')

rowparse <- function(
        cyto_ref_table,
        ref_table,
        Cyto_sample,
        Con_data,
        transloctable,
        Dump_table,
        constitutional,
        guess,
        guess_q,
        forMtn,
        orOption,
        sexstimate
) {

    sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sample_table) <- c("Chr", "Start", "End", "Type")

    sorted_sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sorted_sample_table) <- c("Chr", "Start", "End", "Type")
  
    # Temporary table for storage for mutations
    temp_table <- matrix(
        byrow = TRUE,
        nrow = 1,
        ncol = 4
    )
    colnames(temp_table) <- c("Chr", "Start", "End", "Type")
  
  
    normX <- 0   # normal number of X chormosomes
    normY <- 0   # normal number of Y chromosomes
    xcount <- 0  # counts number of x in 2nd slot
    ycount <- 0  # counts number of y in 2nd slot
    xadd <- 0    # counts if +X occurs
    yadd <- 0    # counts if +Y occurs
    xmod <- 0    # counts modifications that arent whole chromosome add/del for X
    ymod <- 0    # counts modifications that arent whole chromosome add/del for Y
    xdel <- 0    # counts if -X occures as constitutional
    ydel <- 0    # counts if -Y occures as consitutional
  
    xconstitutional <- 0  # shift counts for xc indications (kind of a correction factor)
    yconstitutional <- 0  # shift counts for yc indications (kind of a correction factor)
  
    idealx <- 0  # estimate of what the x value should be
    idealy <- 0  # estimate of what the y value should be
  
    addtot <- 0  # counts total "new chromosomes"
    deltot <- 0  # counts total complete chrom deletions
    modtot <- 0  # for idems only, counts modification chromosomes
  
    n <- 1       # ploidy count
    ploidy <- 2  # ploidy non additive ##default 2 for diploid
  
    startcol <- 2



    if (
        length(Cyto_sample) > 0
        & grepl(
            "[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
            Con_data[2]
        )
    ) {

        # Initiate clonal evolution counter if ids is present
        if (grepl("^ids$", Cyto_sample[length(Cyto_sample)])) {
            # dont count sex chromosomes here
            # this has 46 chromosomes, as they get deleted or modified (without \\+), 
            # knocks em out of this tracker, stuff that remains is either gained 
            # or lost (likely gained)
            clone_chrom_tracker <- rep(1:22, 2)
        }
    
        if (
            grepl(
                "[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
                Con_data[2]
            )
        ) {
            # set normal count for XY chromosomes
            # think about what 46,X,+Y would mean
            if (grepl("(c$)|(c\\?$)", Cyto_sample[2])) {
                xconstitutional <- stringr::str_count(Cyto_sample[2], "X")
                yconstitutional <- stringr::str_count(Cyto_sample[2], "Y")
        
            }
      
            if (
                any(grepl("Y", Cyto_sample))
                & (sum(grepl("Y", Cyto_sample), na.rm = TRUE) 
                    > sum(grepl("\\+Y", Cyto_sample), na.rm = TRUE))
            ) {
                normX <- 1
                normY <- 1

            } else if (sexstimate == F && any(grepl("\\?", Cyto_sample[2]))) {
                # if ? is a chromosome and sexstimate is off, dont make guesses
                # on constitutional change
                normX <- 2
                normY <- 0

            } else {
                normX <- 2

            }
      
            # count number of XY in 2nd slot
            # make sure 2nd slot is sex chromosomes, change this code later
            if (!grepl("^[XY]+", Cyto_sample[2])) {
                # make sure reg ecpression means what you want it to

            } else {
                xcount <- stringr::str_count(Cyto_sample[2], "X")
                xcount <- xcount + stringr::str_count(Cyto_sample[2], "\\?")
                ycount <- stringr::str_count(Cyto_sample[2], "Y")
                # x,y calculation

            }

        }
    
        # check for ploidy levels and digit is over 2
        if (grepl("<[[:digit:]]n>", Cyto_sample[1])) {
            # extract range of stuff before n
            # JP: Redundant assignments
            n <- as.numeric(strsplit(Cyto_sample[1], "<|n>")[[1]][2])
            n <- ploidy
            n <- n - 2

            if (n > 0) {
                # JP: This block is never executed, n always = 0

                for (p in 1:n) {
                    temp_table <- data.frame(
                        ref_table[, 1],
                        rep(0, nrow(ref_table)),
                        ref_table[, 2],
                        rep("Gain", nrow(ref_table))
                    )

                    temp_table <- temp_table[grep("chrY|chrX|chrM", temp_table[, 1], invert = T),]
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
          
                    # not accoutning for sex chromosomes (will acount for later)
                    addtot <- addtot + 22
                }

            }
            # n= # change ploidy count for x and y calculations

        }
    
    
        for (j in startcol:length(Cyto_sample)) {
            # temporary table for storage for mutations
            temp_table <- matrix(
                byrow = TRUE,
                nrow = 1,
                ncol = 4
            )
      
            # for if deletions get completely eliminated in t() type abberations
            deletions = F
      
            #########
            # if guess is true, try to process ? marks
            # think about how this can affect counting  + and \\?
            if (guess_q == T) {
                Cyto_sample[j] <- gsub("\\?", "", Cyto_sample[j])
            }
      
            # if or, delete second option based on boolean
            if (orOption == T) {
                Cyto_sample[j] <- gsub("or.*$", "", Cyto_sample[j])
            }
      
            # addition and deletions of entire chromosomes can be handeled easily
            # and separate from the rest
            if (
                grepl(
                    "mar|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*r\\(*[[:digit:]]*\\)*$|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*neo[[:digit:]]*$",
                    Cyto_sample[j]
                )
                || (constitutional == F && grepl("(c$)|(c\\?$)", Cyto_sample[j]))
            ) {
                if (grepl("\\+", Cyto_sample[j])) {
                    tem = 1
                    # figure out how to do this
                    if (grepl("-|~", Cyto_sample[j])) {
                        if (grepl("mar", Cyto_sample[j])) {
                            Cyto_sample[j] <- paste(
                                unlist(strsplit(Cyto_sample[j], "~|-"))[1],
                                "mar",
                                sep = ""
                            )

                        }

                        if (
                            grepl(
                                "^\\+*([[:digit:]]((~|-)[[:digit:]])*)*r\\(*[[:digit:]]*\\)*$",
                                Cyto_sample[j]
                            )
                        ) {
                            Cyto_sample[j] <- paste(
                                unlist(strsplit(Cyto_sample[j], "~|-"))[1],
                                "r",
                                sep = ""
                            )

                        }
            
                        if (grepl("neo", Cyto_sample[j])) {
                            Cyto_sample[j] <- paste(
                                unlist(strsplit(Cyto_sample[j], "~|-"))[1],
                                "neo",
                                sep = ""
                            )

                        }
            
                    }
          
                    if (grepl("\\+[[:digit:]]", Cyto_sample[j])) {
                        if (constitutional == F & grepl("\\+[[:digit:]]+c\\?*$", Cyto_sample[j])) {
                            # just one addition
                            tem <- 1

                        } else {
                            if (!grepl("\\+[[:digit:]]+c\\?*$", Cyto_sample[j])) {
                                tem <- as.numeric(
                                    strsplit(
                                        strsplit(
                                            Cyto_sample[j],
                                            "mar|r\\(*[[:digit:]]*\\)*$|neo|c$|c\\?$"
                                        )[[1]][1],
                                        "\\+"
                                    )[[1]][2]
                                )

                            }

                        }

                    }
              
                    # only if tem is a number
                    if (is.numeric(tem) || is.integer(tem) || is.double(tem)) {
                        addtot <- addtot + tem

                    } else {
                        # output an error
                        Dump_table <- rbind(
                            Dump_table,
                            c(
                                Con_data,
                                "Error in markers and other ambiguous objects not accounted for"
                            )
                        )
                
                    }

                }

            } else if (
                grepl("^\\+[[:digit:]]+c*$", Cyto_sample[j])
                | grepl("\\+X", Cyto_sample[j])
                | grepl("\\+Y", Cyto_sample[j])
            ) {
                cytoName <- gsub("c", "", substring(Cyto_sample[j], first = 2))
                chr_name <- ref_table[
                    grep(paste("chr", as.character(cytoName), "$", sep = ""), ref_table),
                ]

                temp_table[1, 1] <- chr_name[1]
                temp_table[1, 2] <- "0"
                temp_table[1, 3] <- chr_name[2]
                temp_table[1, 4] <- "Gain"

                # if whole additions occur
                if (grepl("X", chr_name[1])) {
                    xadd <- xadd + 1
                }
            
                if (grepl("Y", chr_name[1])) {
                    yadd <- yadd + 1
                }

                addtot <- addtot + 1

            } else if (
                grepl("^-[[:digit:]]+c*$", Cyto_sample[j])
                | grepl("-X", Cyto_sample[j])
                | grepl("-Y", Cyto_sample[j])
            ) {
                cytoName <- gsub("c", "", substring(Cyto_sample[j], first = 2))
                chr_name <- ref_table[
                    grep(paste("chr", as.character(cytoName), "$", sep = ""), ref_table),
                ]
            
                # exclude x and y deletions until the end
                if (
                    !grepl("Y", chr_name[1])
                    & !grepl("X", chr_name[1])
                ) {
                    temp_table[1, 1] <- chr_name[1]
                    temp_table[1, 2] <- "0"
                    temp_table[1, 3] <- chr_name[2]
                    temp_table[1, 4] <- "Loss"

                }
            
                # if deletions occur in X or Y, up count for the respective mutation
                if (grepl("X", chr_name[1])) {
                    xdel <- xdel + 1

                }
            
                if (grepl("Y", chr_name[1])) {
                    ydel <- ydel + 1

                }
            
                deltot <- deltot + 1

            } else {
                # for all other cases call this first
                # check for x modifications in parser

                inc_table <- tryCatch({
                    mod_colparser$colparse(
                        cyto_ref_table,
                        ref_table,
                        j,
                        xmod,
                        ymod,
                        transloctable,
                        addtot,
                        Cyto_sample,
                        guess_q,
                        constitutional,
                        forMtn
                    )
                }, error = function(e) {
                    return(gsub("\n", " ", paste(e, "in", j, "field")))
                }, finally = {
                    # print(paste("  Parsed field: ", Cyto_sample[j]))
                })

                if (is.character(inc_table)) {
                    Dump_table <- rbind(Dump_table, c(Con_data, inc_table))

                } else if (is.null(inc_table)) {
                  
                } else if (length(inc_table) == 1 && is.na(inc_table)) {
                    Dump_table <- rbind(
                        Dump_table,
                        c(
                            Con_data,
                            "Error in more than one band associated with a chromosome in a translocation"
                        )
                    )
                  
                } else {
                    temp_table <- inc_table[[1]]
                    original_temp_table <- inc_table[[1]]
                    ex_table <- inc_table[[2]]
                    original_ex_table <- inc_table[[2]]
                    xmod <- inc_table[[3]]
                    ymod <- inc_table[[4]]
                    Mainchr <- inc_table[[5]]
                    multi <- inc_table[[6]]
                    transloctable <- inc_table[[7]]
                    addtot <- inc_table[[8]]
                  
                    # print(multi)
                    # for future implementation: switch to searching at ex table if temp table is empty,
                    # default is temp_table
                    pointerConditional <- temp_table
                  
                    if (nrow(temp_table) > 0) {
                        # if its a translocated vector, fix now
                        if (nrow(temp_table) == 4 & ncol(temp_table) == 1) {
                            temp_table <- t(temp_table)
                            colnames(temp_table) <- c("Chr", "Start", "End", "Type")

                        }
                    
                        # special cases (if, else if, else if ,else if, else )
                        # ider, i , ins, dic etc
                        # else, grep stuff on 4th col to detirmine addition or deletions, dont include
                        #   translations add etc
                        # deal with + values differently
                    
                        # if it goes through special cases, flip mod = true, excoord does not combine at end
                        mod <- FALSE
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        colnames(temp_table) <- c("Chr", "Start", "End", "Type")
                    
                        if (length(ex_table) > 0) {
                            colnames(ex_table) <- c("Chr", "Start", "End", "Type")
                            ex_table[, 4] <- as.character(ex_table[, 4])
                            ex_table[, 2:3] <- apply(
                                ex_table[, 2:3],
                                2,
                                function(x) {
                                    as.numeric(as.character(x))
                                }
                            )
                      
                        }

                        temp_table[, 2:3] <- apply(
                            temp_table[, 2:3],
                            2,
                            function(x) {
                                as.numeric(as.character(x))
                            }
                        )
                    
                        # see if its an addition
                        additionTable <- mod_utils$detectAdd(
                            temp_table[, 4],
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4]
                            } else {
                                NULL
                            }
                        )

                        # if its not a del # think about this one for all cases
                    
                        # need to improve counter, maybe put it way upstream in each j so if it finds +,
                        #   adds, then take it out of everything else
                        # probbbly doesnt take into account x3 (multi -2)
                        if (grepl("\\+", Cyto_sample[j])) {
                            addtot <- (1 * multi) + addtot

                        } else if (multi > 2) {
                            addtot <- addtot + multi - 2

                        }

                        # temp_table[grep("+\(.\)",ex_table[,4]),4]<-"Gain"
                        # ex_table[,4]<-"Gain"
                        # rbind(temp_table,ex_table)
                        # incorporate main chr in these
                        # do del and stuff up here, but keep+__ sign somehow
                    
                        # count chromosome loss properly here for special cases
                        # sex chromosome count check how it interacts with xmod ymod
                        if (
                            !all(grepl("\\+", temp_table[, 4]))
                            & any(!grepl("^t\\(|^ins\\(|^inv\\(", temp_table[, 4]))
                        ) {
                            # takes mainchr off of list of chromosomes because it has been
                            # acccounted for (idems only)
                            if (grepl("^ids$", Cyto_sample[length(Cyto_sample)])) {
                                if (length(Mainchr[-(grep("X|Y", Mainchr))]) > 0) {
                                    tempindexchromtracker <- sapply(
                                        gsub("\\$", "", Mainchr),
                                        function(x) {
                                            if (!grepl("X|Y", x)) {
                                                which.max(clone_chrom_tracker == as.numeric(x))
                                            }
                                        }
                                    )
                          
                                    # fix up if its a list
                                    if (is.list(tempindexchromtracker)) {
                                        tempindexchromtracker <- unlist(tempindexchromtracker)

                                    }
                          
                                    if (
                                        !is.null(tempindexchromtracker)
                                        && length(tempindexchromtracker) > 0
                                    ) {
                                        clone_chrom_tracker <-
                                            clone_chrom_tracker[-1 * tempindexchromtracker]

                                    }

                                }

                            }
                      
                            # check how this handles X chromosomes
                            # if more than 2 chromosomes involved and its not a translocation an 
                            # insertion or an inversion, something is "deleted"
                            if (
                                any(grepl("-[:alpha:]", temp_table[, 4]))
                                & any(!grepl("X\\$|Y\\$", Mainchr))
                            ) {
                                # deltot <- deltot + (1 * multi)

                            } else {
                                if (length(Mainchr) != length(grep("X|Y", Mainchr))) {
                                    deltot <- deltot + 
                                        (length(grep("X|Y", Mainchr, invert = T)) - 1) * multi

                                }
                        
                            }
                      
                        }
                    
                        # insertions affect count if in derivative chromosome or the like
                    
                        # take any translocations and insertions and stores for later
                        # if (any(grepl("(t|ins)\\(", ex_table[, 4])))
                        # {
                        #   transloctable <-
                        #   rbind(transloctable, ex_table[grep("(t|ins)\\(", ex_table[, 4]), ])
                        # }
                    
                        if (any(grepl("LongDer", temp_table[, 4]))) {
                            mod <- TRUE
                            temp_table[
                                intersect(
                                    grep("LongDer", temp_table[, 4]),
                                    grep(
                                        paste(Mainchr, collapse = "|", sep = ''),
                                        temp_table[, 1],
                                        invert = T
                                    )
                                ),
                                4
                            ] <- "Gain"

                            if (nrow(ex_table) > 0) {
                                ex_table[grep("LongDer", ex_table[, 4]), 4] <- "Loss"

                            }

                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        # Think about how this is handling issue
                        #   Really check on this
                        if (
                            any(
                                grepl("((t)|(ins))\\(", temp_table[, 4])
                                & !grepl("^((t)|(ins))\\(", temp_table[, 4])
                            )
                        ) {
                            mod <- TRUE
                            temp_table[
                                intersect(
                                    intersect(
                                        grep(
                                            paste("chr", Mainchr, sep = "", collapse = "|"),
                                            temp_table[, 1],
                                            invert = T
                                        ),
                                        grep("(t|ins)\\(", temp_table[, 4])
                                    ),
                                    grep("^((t)|(ins))\\(", temp_table[, 4], invert = T)
                                ),
                                4
                            ] <- "Gain"

                            if (nrow(ex_table) > 0) {
                                ex_table[
                                    intersect(
                                        intersect(
                                            grep(
                                                paste("chr", Mainchr, sep = "", collapse = "|"),
                                                ex_table[, 1]
                                            ),
                                            grep("((t)|(ins))\\(", ex_table[, 4])
                                        ),
                                        grep("^((t)|(ins))\\(", ex_table[, 4], invert = T)
                                    ),
                                    4
                                ] <- "Loss"

                            }

                        }
                    
                        # Take deletion areas out of translocation adjacent things
                        if (
                            any(
                                grepl("del", temp_table[, 4])
                                & !grepl("^del", temp_table[, 4])
                                & !grepl("multi[[:digit:]]*del", temp_table[, 4])
                            )
                            | any(
                                grepl("add", temp_table[, 4])
                                & !grepl("^add", temp_table[, 4])
                                & !grepl("multi[[:digit:]]*add", temp_table[, 4])
                            )
                        ) {
                            if (
                                nrow(temp_table)
                                    > length(
                                        union(
                                            intersect(
                                                grep("del", temp_table[, 4]),
                                                intersect(
                                                    grep("^del", temp_table[, 4], invert = T),
                                                    grep("multi[[:digit:]]*del", temp_table[, 4], invert = T)
                                                )
                                            ), 
                                            intersect(
                                                grep("add", temp_table[, 4]),
                                                intersect(
                                                    grep("^add", temp_table[, 4], invert = T),
                                                    grep("multi[[:digit:]]*add", temp_table[, 4], invert = T)
                                                )
                                            )
                                        )
                                    )
                            ) {
                                temp_table <- mod_merge$mergeDeletions(temp_table, Mainchr)
                                deletions <- T
                                # if it is a vector, convert
                                if (is.vector(temp_table)) {
                                    temp_table <- as.data.frame(as.list(temp_table))
                                    temp_table[, 2:3] <- as.numeric(temp_table[, 2:3])
                                    colnames(temp_table) <- c("Chr", "Start", "End", "Type")

                                } else if (nrow(temp_table) == 4 & ncol(temp_table) == 1) {
                                    # if its a translocated vector, fix now
                                    temp_table <- t(temp_table)
                                    temp_table <- as.data.frame(as.list(temp_table))
                                    temp_table[, 2:3] <- as.numeric(temp_table[, 2:3])
                                    colnames(temp_table) <- c("Chr", "Start", "End", "Type")

                                }
                        
                                additionTable <- mod_utils$detectAdd(
                                    temp_table[, 4],
                                    if (nrow(ex_table) > 0) {
                                        ex_table[, 4]

                                    } else {
                                        NULL

                                    }
                                )
                        
                            }
                      
                        }
                    
                        # if temp table is empty, skip all this
                        if (nrow(temp_table) > 0) {
                      
                        }

                        if ((
                            any(grepl("((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))", temp_table[, 4]))
                            & length(temp_table[, 4])
                                > length(
                                    c(
                                        which(grepl("del", temp_table[, 4])),
                                        which(grepl("add", temp_table[, 4]))
                                    )
                                )
                            & length(temp_table[, 4]) > 1
                            & (
                                length(temp_table[, 4])
                                - length(
                                    c(
                                        which(grepl("del", temp_table[, 4])),
                                        which(grepl("add", temp_table[, 4]))
                                    )
                                )
                            ) > 1
                        ) || (
                            any(grepl("((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))", temp_table[, 4]))
                            && deletions
                        )) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                temp_table[
                                    intersect(
                                        grep(
                                            paste(paste("chr", Mainchr, sep = ""), collapse = "|"),
                                            temp_table[, 1],
                                            invert = T
                                        ),
                                        grep(
                                            "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                                            temp_table[, 4]
                                        )
                                    ),
                                    4
                                ] <- "Gain"

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[
                                    intersect(
                                        grep(
                                            paste(paste("chr", Mainchr, sep = ""), collapse = "|"),
                                            ex_table[, 1]
                                        ),
                                        grep(
                                            "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                                            ex_table[, 4]
                                        )
                                    ),
                                    4
                                ] <- "Loss"

                            }
                      
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }
                      
                            temp_table <- rbind(temp_table, ex_table)
                      
                        }
                    
                        if (any(grepl("i\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                original_temp_table <- temp_table
                                temp_table[grep("i\\(.*", temp_table[, 4]), 4] <- "Gain"
                                additionTable[[1]] <- c(
                                    additionTable[[1]],
                                    additionTable[[1]][grep("i\\(.*", original_temp_table[, 4])]
                                )
                                temp_table <- rbind(
                                    temp_table,
                                    original_temp_table[grep("i\\(.*", original_temp_table[, 4]),]
                                )
                                # add to additionTable in light of new addition

                            }

                            if (nrow(ex_table) > 0) {
                                ex_table[grep("i\\(.*", ex_table[, 4]), 4] <- "Loss"

                            }

                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }

                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            temp_table <- rbind(temp_table, ex_table)

                        }

                        # get whats included, dup, then rest on chromosome is deletion
                    
                        if (any(grepl("ider\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            # Have to be able to handle these cases
                            # Do inclusion and exclusion based on deleted case
                            if (nrow(temp_table) > 0) {
                                original_temp_table <- temp_table
                                temp_table[
                                    intersect(
                                        grep("ider\\(.*", temp_table[, 4]),
                                        grep("del\\(.*|add\\(.*", temp_table[, 4], invert = T)
                                    ),
                                    4
                                ] <- "Gain"

                                additionTable[[1]] <- c(
                                    additionTable[[1]],
                                    additionTable[[1]][
                                        intersect(
                                            grep("ider\\(.*", temp_table[, 4]),
                                            grep("del\\(.*|add\\(.*", temp_table[, 4], invert = T)
                                        )
                                    ]
                                )

                                temp_table <- rbind(
                                    temp_table,
                                    original_temp_table[grep("ider\\(.*", original_temp_table[, 4]),]
                                )
                        
                                # activate deletions ? Think about this
                                # JAN Thinking time
                            }

                            if (nrow(ex_table) > 0) {
                                ex_table[grep("ider\\(.*", ex_table[, 4]), 4] <- "Loss"

                            }
                      
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            temp_table <- rbind(temp_table, ex_table)
                      
                        }
                    
                        if (any(grepl("idic\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            # Need to either break it up or replace,
                            # think about this one too
                            # deleted <- ex_table[grep("idic\\(.*&del\\(", ex_table[, 4]), 4]
                            if (nrow(temp_table) > 0) {
                                original_temp_table <- temp_table
                                temp_table[
                                    intersect(
                                        grep("idic\\(.*", temp_table[, 4]),
                                        grep("dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", temp_table[, 4], invert = T)
                                    ),
                                    4
                                ] <- "Loss"

                                additionTable[[1]] <- c(
                                    additionTable[[1]],
                                    additionTable[[1]][
                                        intersect(
                                            grep("del\\(.*|add\\(.*", original_temp_table[, 4], invert = T),
                                            grep("idic\\(.*", original_temp_table[, 4])
                                        )
                                    ]
                                )

                                temp_table <- rbind(
                                    temp_table,
                                    original_temp_table[
                                        intersect(
                                            grep("del\\(.*|add\\(.*", original_temp_table[, 4], invert = T),
                                            grep("idic\\(.*", original_temp_table[, 4])
                                        )
                                        ,
                                    ]
                                )

                            }

                            if (nrow(ex_table) > 0) {
                                original_ex_table <- ex_table
                                ex_table[grep("idic\\(.*", ex_table[, 4]), 4] <- "Gain"
                                additionTable[[2]] <- c(
                                    additionTable[[2]],
                                    additionTable[[2]][grep("idic\\(.*", original_ex_table[, 4])]
                                )
                                ex_table <- rbind(
                                    ex_table,
                                    original_ex_table[grep("idic\\(.*", original_ex_table[, 4]),]
                                )

                            }

                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }

                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        # figure out why you did this
                        # make sure this is reversed for long form
                        if (
                            any(
                                which(
                                    grepl("dic\\(.*", temp_table[, 4])
                                    & !grepl("idic\\(.*", temp_table[, 4])
                                )
                            )
                        ) {

                            mod <- TRUE
                            if (any(grepl("long", temp_table[, 4]))) {

                                if (nrow(temp_table) > 0) {
                                    temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                                }

                                if (nrow(ex_table) > 0) {
                                    ex_table[grep("dic\\(.*", ex_table[, 4]) , 4] <- "Loss"
                                    ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                                }
                        
                            } else {
                                if (nrow(temp_table) > 0) {
                                    temp_table[
                                        intersect(
                                            grep("dic\\(.*", temp_table[, 4]),
                                            grep(
                                                "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                                                temp_table[, 4],
                                                invert = T
                                            )
                                        ),
                                        4
                                    ] <- "Loss"
                          
                                    temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
                                
                                }

                                if (nrow(ex_table) > 0) {
                                    ex_table[
                                        intersect(
                                            intersect(
                                                grep("dic\\(.*", ex_table[, 4]),
                                                grep(
                                                    "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                                                    ex_table[, 4],
                                                    invert = T
                                                )
                                            ),
                                            grep("del\\(.*|add\\(.*", ex_table[, 4])
                                        ),
                                        4
                                    ] <- gsub(
                                        "del\\(.*|add\\(.*",
                                        "",
                                        ex_table[
                                            intersect(
                                                intersect(
                                                    grep("dic\\(.*", ex_table[, 4]),
                                                    grep(
                                                        "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                                                        ex_table[, 4],
                                                        invert = T
                                                    )
                                                ),
                                                grep("del\\(.*|add\\(.*", ex_table[, 4])
                                            ),
                                            4
                                        ]
                                    )
                                    ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")
                          
                                }

                            }
                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        # make sure this is reversed for long form
                        if (any(grepl("trc\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            if (any(grepl("long", temp_table[, 4]))) {
                                if (nrow(ex_table) > 0) {
                                    ex_table[
                                        intersect(
                                            grep("trc\\(.*", ex_table[, 4]),
                                            grep(
                                                paste("chr", Mainchr[1], "chr", Mainchr[3], sep = '|'),
                                                ex_table[, 1]
                                            )
                                        ),
                                        4
                                    ] <- "Loss"

                                }

                            } else {
                                if (nrow(temp_table) > 0) {
                                    temp_table[
                                        intersect(
                                            grep("trc\\(.*", temp_table[, 4]),
                                            grep(
                                                paste("chr", Mainchr[1], "chr", Mainchr[3], sep = '|'),
                                                temp_table[, 1]
                                            )
                                        ),
                                        4
                                    ] <- "Loss"

                                }

                            }
                      
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[
                                    intersect(
                                        grep("trc\\(.*", ex_table[, 4]),
                                        grep("chr", Mainchr[2], ex_table[, 1])
                                    ),
                                    4
                                ] <- "Loss"
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }
                      
                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        if (any(grepl("rob\\(", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(ex_table) > 0) {
                                ex_table[grep("rob\\(", ex_table[, 4]), 4] <- "Loss"

                            }
                      
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }

                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        if (any(grepl("^r\\(.*|^\\+r\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(ex_table) > 0) {
                                ex_table[grep("r\\(.*", ex_table[, 4]), 4] <- "Loss"

                            }

                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        if (any(grepl("trp\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            # JP: ex_table not in rbind??
                            temp_table <- rbind(temp_table, temp_table[grep("trp\\(.*", temp_table[, 4]),])
                      
                        }
                    
                        if (any(grepl("qdp\\(.*", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }

                            # JP: why is this duplicated in rbind?
                            temp_table <- rbind(
                                temp_table,
                                temp_table[grep("qdp\\(.*", temp_table[, 4]),],
                                temp_table[grep("qdp\\(.*", temp_table[, 4]),]
                            )

                        }
                    
                    
                        # Add code here for +gain +loss (+gain=gain, +loss ="")
                        if (any(grepl("del", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                temp_table[grep("del", temp_table[, 4]), 4] <- "Loss"
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }

                            if (nrow(ex_table) > 0) {
                                if (is.vector(ex_table)) {
                                    ex_table[, 4] <- ""

                                } else {
                                    ex_table[, 4] <- ""

                                }
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }
                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        if (any(grepl("add", temp_table[, 4]))) {
                            mod <- TRUE
                            if (nrow(temp_table) > 0) {
                                temp_table[grep("add", temp_table[, 4]), 4] <- "Loss"
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }

                            if (nrow(ex_table) > 0) {
                                if (is.vector(ex_table)) {
                                    ex_table[, 4] <- ""

                                } else {
                                    ex_table[, 4] <- ""

                                }
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }
                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        # keep X and Y chromosomes here for ex table for later processing, may have to 
                        # modify later for autosomes, only if modifications aren't triggered
                        if (mod == FALSE) {
                            if (is.vector(ex_table)) {
                                ex_table[4] <- ""

                            } else {
                                if (nrow(ex_table) > 0) {
                                    ex_table[, 4] <- ""

                                }

                            }

                            if (nrow(temp_table) > 0) {
                                temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")

                            }
                      
                            if (nrow(ex_table) > 0) {
                                ex_table[, 4] <- paste(additionTable[[2]], ex_table[, 4], sep = "")

                            }
                            temp_table <- rbind(temp_table, ex_table)

                        }
                    
                        # handle stuff for minus as well
                        # if +Gain
                        if (nrow(temp_table) > 0) {
                            temp_table[grep("^-Gain|^-$", temp_table[, 4]), 4] <- "Loss"
                            temp_table[grep("\\+Loss", temp_table[, 4]), 4] <- "+Loss"
                            temp_table[grep("\\+Gain", temp_table[, 4]), 4] <- "Gain"
                            temp_table[
                                intersect(
                                    grep("\\+", temp_table[, 4]),
                                    grep("\\+Loss", temp_table[, 4], invert = T)
                                ),
                                4
                            ] <- "Gain"

                            # check other permutations of this/ dont think insertions belong here
                            temp_table[grep("dup|qdp|tan|trp|\\+$", temp_table[, 4]), 4] <- "Gain"
                      
                            Plus_Loss <- temp_table[grep("\\+Loss", temp_table[, 4]), ]
                      
                            # cut off intersections early
                            # Think aboiut this
                            if (
                                length(Plus_Loss) > 0
                                && nrow(Plus_Loss) > 0
                                && any(grepl("Gain", temp_table[, 4]) & any(grepl("Loss", temp_table[, 4])))
                            ) {
                                temp_table <- mod_merge$mergeTable(temp_table)

                            }
                      
                            temp_table[grep("\\+Loss", temp_table[, 4]), 4] <- ""

                        }
                    }
                }
            }
          
            # combine running file and temp new file together
            if (length(temp_table) != 0) {
                temp_table <- apply(temp_table, 2, as.character)
                sample_table <- rbind(sample_table, temp_table)

            }

        } # for (j in startcol:length(Cyto_sample))
        
    
    
        # Total count, ensure no wild stuff happened
        # have to account for ranges,
        # have to count markers
        val <- strsplit(Cyto_sample[1], "<")[[1]] # number of chromosomes indicated in first value
        val <- paste(strsplit(val, "[[:alpha:]]+")[[1]], sep = "", collapse = "")

        if (grepl("~|-", val)) {
            val <- unlist(strsplit(val, "~|-"))
            if (all(!is.na(as.numeric(val)))) {
                val <- as.numeric(val)
        
                if (as.numeric(val[2]) >= as.numeric(val[1])) {
                    val <- seq(from = val[1], to = val[2], by = 1)

                } else {
                    val <- seq(from = val[2], to = val[1], by = 1)

                }

            } else {
                val <- as.numeric(val)
                val <- val[2]

            }
      
        }
    
        val <- as.numeric(val)
        # print(paste("val 1: ", val))
        # prelim ploidy gueeses for xideal and yideal from val
        # take x and y completely out of this
        val <- val - xcount - ycount - (xadd + yadd - xdel - ydel)
        if (normX + normY != 2) {
            if (normX + normY > 2) {
                val <- val + (normX + normY - 2)

            } else {
                val <- val + (2 - normX + normY)

            }

        }
        # print(paste("val 2: ", val))
        
        
        # original val
        orgval <- strsplit(Cyto_sample[1], "<")[[1]] # number of chromosomes indicated in first value
        orgval <- paste(strsplit(orgval, "[[:alpha:]]+")[[1]], sep = "", collapse = "")

        if (grepl("~|-", orgval)) {
            orgval = unlist(strsplit(orgval, "~|-"))
            if (all(!is.na(as.numeric(orgval)))) {
                if (all(!is.na(as.numeric(orgval)))) {
                    orgval <- as.numeric(orgval)
                    if (as.numeric(orgval[2]) >= as.numeric(orgval[1])) {
                        orgval <- seq(from = orgval[1], to = orgval[2], by = 1)

                    } else {
                        orgval <- seq(from = orgval[2], to = orgval[1], by = 1)

                    }

                }

            } else {
                orgval <- as.numeric(orgval)
                orgval <- orgval[2]

            }
            
        }
    
        if (any(is.na(val))) {
            Dump_table <- rbind(
                Dump_table,
                c(as.vector(Con_data), "Error in unclear chrom number")
            )

        } else {
            idealval = 46 + addtot - deltot - xcount - ycount -
                (xadd + yadd - xdel - ydel)
            #idealval = 46 + addtot - deltot - xcount - ycount - (xadd + yadd - xdel -
            #                                                       ydel)            
            # print(paste(addtot, deltot, xcount, ycount, xadd, yadd, xdel, ydel))
            # print(paste(normX, normY))
            # print(paste("idealval 0: ", idealval))
            if (normX + normY != 2) {
                if (normX + normY > 2) {
                    idealval = idealval + (normX + normY - 2)

                } else {
                    idealval = idealval + (2 - normX + normY)

                }

            }
            # print(paste("idealval 1: ", idealval))
            
            diffval <- val - idealval
            val_remainder <- diffval %% 22
            val_remainder_2 <- diffval %% 23
            val_divider <- diffval / 22 ## for stuff over diploidy
            val_divider_2 <- diffval / 23
            val_temp <- 0
      
            if (!is.null(val_remainder)) {
                if (any(val_remainder == 0)) {
                    val_divider <- val_divider[grep(TRUE, val_remainder == 0)]
                    val_remainder <- val_remainder[grep(TRUE, val_remainder == 0)]
          
                } else if (any(val_remainder_2 == 0)) {
                    val_divider <- val_divider_2[grep(TRUE, val_remainder_2 == 0)]
                    val_remainder <- val_remainder_2[grep(TRUE, val_remainder_2 == 0)]

                } else {
                    val_remainder <- val_remainder[1]
                    val_divider <- val_divider[1]

                }
        
                # print(c(val, idealval, val_remainder, val_divider))
        
        
                if (val_remainder == 0) {
                    if (val_divider > 0) {
                        for (f in 1:(val_divider)) {
                            temp_table <- data.frame(
                                ref_table[, 1],
                                rep(0, nrow(ref_table)),
                                ref_table[, 2],
                                rep("Gain", nrow(ref_table))
                            )
                            temp_table <- temp_table[
                                grep("chrY|chrX|chrM", temp_table[, 1], invert = T),
                            ]
              
                            temp_table[, 4] <- as.character(temp_table[, 4])
                            temp_table <- as.matrix(temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                            sample_table <- rbind(sample_table, temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                        }
            
                        ploidy <- val_divider + 2
            
                        #print(c("val_div>0", Cyto_sample))
                    }
          
                    if (val_divider < 0) {
                        if (val_divider == -1) {
                            temp_table <- data.frame(
                                ref_table[, 1],
                                rep(0, nrow(ref_table)),
                                ref_table[, 2],
                                rep("Loss", nrow(ref_table))
                            )

                            temp_table <- temp_table[
                                grep("chrY|chrX|chrM", temp_table[, 1], invert = T),
                            ]
              
                            temp_table[, 4] <- as.character(temp_table[, 4])
                            temp_table <- as.matrix(temp_table)
                            colnames(temp_table) <- colnames(sample_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                            sample_table <- rbind(sample_table, temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
              
                            ploidy = 1
              
                        }

                        #print(c("val_div<0", Cyto_sample))
                    }
          
                } else if (
                    grepl("^ids$", Cyto_sample[length(Cyto_sample)])
                    && any(
                        diffval / length(clone_chrom_tracker) > 0
                        & diffval %% length(clone_chrom_tracker) == 0
                    )
                ) {
          
                    # if unaccounted chromosom number == diff val, add them,
                    # let this estimate for uncertainty
          
                    diffval <- diffval[
                        intersect(
                            grep(TRUE, diffval / length(clone_chrom_tracker) > 0),
                            grep(TRUE, diffval %% length(clone_chrom_tracker) == 0)
                        )
                    ]
                    ploidy <- (diffval / length(clone_chrom_tracker)) + 2
          
                    for (k in 1:diffval / length(clone_chrom_tracker)) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[
                            grep(
                                paste("chr", clone_chrom_tracker, collapse = '|'),
                                temp_table[, 1]
                            ),
                        ]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                    }
          
                } else if (
                    guess == T
                    & (
                        any(diffval < (-14))
                        | any(diffval < 31 & diffval > 14)
                        | any(diffval < 53 & diffval > 38)
                        | any(
                            diffval > 56
                            & (floor(val_divider) * 22 + 3 < diffval)
                            & (ceiling(val_divider) * 23 - 3 > diffval)
                        )
                    )
                ) {
                    # clone_chrom_tracker <- rep(1:22, 2)
                    if (length(diffval) > 1) {
                        # take the one that fits and first one that is true
                        new_diffval <- diffval[
                            which(
                                diffval < (-14)
                                | (diffval < 31 & diffval > 14)
                                | (diffval < 53 & diffval > 38)
                                | (
                                    diffval > 56
                                    & (floor(val_divider) * 22 + 3 < diffval)
                                    & (ceiling(val_divider) * 23 - 3 > diffval)
                                )
                            )
                        ][1]

                    } else {
                        new_diffval <- diffval

                    }

                    # think about how to implement this beyond teraploidy
                    # doesnt work
                    # print(new_diffval)
                    if (new_diffval < (-14)) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Loss", nrow(ref_table))
                        )
                        temp_table <- temp_table[1:22, ]

                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        ploidy <- 1

                    }


                    if (new_diffval < 31 & new_diffval > 14) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[1:22, ]

                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        ploidy <- 3

                    }

                    if (new_diffval < 53 & new_diffval > 38) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[1:22, ]

                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])

                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[1:22, ]

                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        ploidy <- 4

                    }

                    if (
                        new_diffval > 56
                        & floor(val_divider) * 22 + 3 < new_diffval
                        & ceiling(val_divider) * 22 - 3 > new_diffval
                    ) {
                        end <- round(val_divider)
                        ploidy <- end + 2
                        for (k in 1:end) {
                            temp_table <- data.frame(
                                ref_table[, 1],
                                rep(0, nrow(ref_table)),
                                ref_table[, 2],
                                rep("Gain", nrow(ref_table))
                            )
                            temp_table <- temp_table[1:22, ]

                            temp_table[, 4] <- as.character(temp_table[, 4])
                            temp_table <- as.matrix(temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                            sample_table <- rbind(sample_table, temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                        }

                    }


                    # time to guess according to ranges if guessing is true

                    # put in dump table
                    Dump_table <- rbind(
                        Dump_table,
                        c(
                            as.vector(Con_data),
                            "Warning in some chromosomes unaccounted for"
                        )
                    )

                    # print(c("unaccounted", Cyto_sample))


                } else if (val_divider > 1) {
                    # if doesnt match initially, try to see if bottom few match
                    val_temp <- floor(val_divider)
                    ploidy <- val_temp + 2
                    if (val_temp > 0) {
                        for (f in 1:(val_temp)) {
                            temp_table <- data.frame(
                                ref_table[, 1],
                                rep(0, nrow(ref_table)),
                                ref_table[, 2],
                                rep("Gain", nrow(ref_table))
                            )
                            temp_table <- temp_table[
                                grep("chrY|chrX|chrM", temp_table[, 1], invert = T),
                            ]
              
                            temp_table[, 4] <- as.character(temp_table[, 4])
                            temp_table <- as.matrix(temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                            sample_table <- rbind(sample_table, temp_table)
                            sample_table[, 4] <- as.character(sample_table[, 4])
                        }
            
                    }
          
                    # print(c("val_div>0", Cyto_sample))
                    # add this to uncertain
          
                    # put in dump table
                    Dump_table <- rbind(
                        Dump_table,
                        c(
                            as.vector(Con_data),
                            "Warning in some chromosomes unaccounted for"
                        )
                    )
                    # print(c("unaccounted", Cyto_sample))
          
                }

            }
                
        }

        # Deal with sex chromosomes here
        # remember to deal with 46,xxx
    
        count_before_extras <- 0 # theoretically original
        ploidy_count <- ploidy # count_before_mods * plody
        count_after_mods <- 0 # counts with deletions
        constitutionalcount <- xconstitutional + yconstitutional
    
        count_before_extras <- xcount + ycount + xmod + ymod
    
        idealTotal <- normX + normY
        sexDev_from_norm <- 0
        difference <- 0
        count_after_mods <- count_before_extras

        # if mitelman specifications are true and there is no sex count, assume
        # constitutionality, just count -X and - Y straight
        if (forMtn == T & (xcount + ycount) == 0) {
            constitutional <- T
            if (ydel >= 1) {
                for (f in 1:ydel) {
                    temp_table <- data.frame(
                        ref_table[, 1],
                        rep(0, nrow(ref_table)),
                        ref_table[, 2],
                        rep("Loss", nrow(ref_table))
                    )
                    temp_table <- temp_table[grep("chrY", temp_table[, 1]), ]
            
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table) <- colnames(sample_table)
            
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                }

            }

            if (xdel >= 1) {
                for (f in 1:xdel) {
                    temp_table <- data.frame(
                        ref_table[, 1],
                        rep(0, nrow(ref_table)),
                        ref_table[, 2],
                        rep("Loss", nrow(ref_table))
                    )
                    temp_table <- temp_table[grep("chrX", temp_table[, 1]), ]
            
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table) <- colnames(sample_table)
            
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                }

            } 

        } else if (constitutionalcount > 0 && constitutional == F) {
            # if the x count is consitutional and we dont want to count the constitutional state,
            # just print -xs as is 
            if (ydel >= 1) {
                for (f in 1:ydel) {
                    temp_table <- data.frame(
                        ref_table[, 1],
                        rep(0, nrow(ref_table)),
                        ref_table[, 2],
                        rep("Loss", nrow(ref_table))
                    )
                    temp_table <- temp_table[grep("chrY", temp_table[, 1]), ]
          
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table) <- colnames(sample_table)
          
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                }

            }

            if (xdel >= 1) {
                for (f in 1:xdel) {
                    temp_table <- data.frame(
                        ref_table[, 1],
                        rep(0, nrow(ref_table)),
                        ref_table[, 2],
                        rep("Loss", nrow(ref_table))
                    )
                    temp_table <- temp_table[grep("chrX", temp_table[, 1]), ]
          
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table) <- colnames(sample_table)
          
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                }

            }

        } else {
            # do the complicated calculations
            # if we need to take into account ploidy
            # handle 69,xx,-y (would be one xgain, one y loss)
            if (
                idealTotal * (ploidy - 2) / 2 + idealTotal != count_before_extras
                | (
                    count_before_extras + xdel + ydel > idealTotal * (ploidy - 2) / 2 + idealTotal
                    & xdel + ydel > 0
                )
            ) {
                if (ploidy > 2) {
                    idealTotal <- idealTotal * (ploidy - 2) / 2 + idealTotal
          
                } else if (ploidy == 1) {
                    idealTotal <- 1

                }
        
                if (constitutionalcount > 0) {
                    count_after_mods <- count_before_extras - xdel - ydel

                } else if (idealTotal == count_before_extras + xdel + ydel) {
                    # think about this, (46,xx,-x,-x)
                    count_after_mods = count_before_extras

                } else if (count_before_extras + xdel + ydel > idealTotal) {
                    # if theres overflow, determine if there is an addition and subtraction
                    # needed to be calculated, or if raw counts should count as total
                    if (
                        count_before_extras + xdel + ydel > idealTotal
                        & xdel + ydel > 0
                    ) {
                        # if the raw counts are greater than idealTotal, assume raw counts
                        # are the total and deletions are subtractions from that value
                        if (count_before_extras >= idealTotal) {
                            count_after_mods <- count_before_extras
              
                            # set equal to ploidy, then subtract difference
                            count_after_mods <- count_after_mods - xdel - ydel

                        } else {
                            # calculate difference
                            count_after_mods <- count_before_extras
              
                            # set equal to ploidy, then subtract difference
                            difference <- xdel + ydel + count_after_mods - ploidy
                            count_after_mods <- count_after_mods - difference

                        }

                    } else if (count_before_extras < idealTotal) {
                        count_after_mods = (ploidy - 2) * count_before_extras
            
                        # set equal to ploidy, then subtract difference
                        difference <- xdel + ydel - count_after_mods
                        count_after_mods <- count_after_mods + difference

                    } else {
                        count_after_mods <- count_before_extras + xdel + ydel
            
                    }
          
                }
        
                if (constitutional == F) {
                    # add a gain/loss to counter to negate loss
                    # if haploidy, add a gain to counteract
                    if (constitutionalcount > 2) {
                        ynew <- 0
            
                        if (yconstitutional > 0) {
                            for (f in 1:yconstitutional - 1) {
                                temp_table <- data.frame(
                                    ref_table[, 1],
                                    rep(0, nrow(ref_table)),
                                    ref_table[, 2],
                                    rep("Loss", nrow(ref_table))
                                )
                                temp_table <- temp_table[grep("chrY", temp_table[, 1]),]
                
                                temp_table[, 4] <- as.character(temp_table[, 4])
                                temp_table <- as.matrix(temp_table)
                                colnames(temp_table) <- colnames(sample_table)
                
                                sample_table[, 4] <- as.character(sample_table[, 4])
                                sample_table <- rbind(sample_table, temp_table)
                                sample_table[, 4] <- as.character(sample_table[, 4])
                            }

                        }
            
                        if (xconstitutional > 0) {
                            for (f in 1:(xconstitutional - xcount)) {
                                temp_table <- data.frame(
                                    ref_table[, 1],
                                    rep(0, nrow(ref_table)),
                                    ref_table[, 2],
                                    rep("Loss", nrow(ref_table))
                                )
                                temp_table <- temp_table[grep("chrX", temp_table[, 1]),]
                
                                temp_table[, 4] <- as.character(temp_table[, 4])
                                temp_table <- as.matrix(temp_table)
                                colnames(temp_table) <- colnames(sample_table)
                
                                sample_table[, 4] <- as.character(sample_table[, 4])
                                sample_table <- rbind(sample_table, temp_table)
                                sample_table[, 4] <- as.character(sample_table[, 4])
                            }

                        }

                    }
          
                }

                # if polyploidy, loss ploidy-2
                # for loop for each constitutional value
            }
      
            sexDev_from_norm <- count_after_mods - 2
      
            if (sexDev_from_norm > 0) {
                # Add gain of sex chrom for loop
                # Copy paste old code

                ynew = 0
                for (f in 1:sexDev_from_norm) {
                    # Include if y and x discordance
                    if ((ycount + ymod) > (ynew + normY)) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrY", temp_table[, 1]), ]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        ynew = ynew + 1

                    } else {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrX", temp_table[, 1]), ]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
            
                    }
                }
        
            } else if (sexDev_from_norm < 0) {
                ynew = 0
                for (f in 1:(-1 * sexDev_from_norm)) {
                    if ((ycount + ymod) < (ynew + normY)) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Loss", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrY", temp_table[, 1]), ]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        ynew = ynew + 1
            
                    } else {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Loss", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrX", temp_table[, 1]), ]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
            
                    }
                }
        
            } else {
                # If sexDev == 0, make sure to take into account -y and XX canceling out
                if ((ydel) >= 1) {
                    for (f in 1:ydel) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Loss", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrY", temp_table[, 1]),]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
                    }
          
                    if (xcount + xmod > 1) {
                        temp_table <- data.frame(
                            ref_table[, 1],
                            rep(0, nrow(ref_table)),
                            ref_table[, 2],
                            rep("Gain", nrow(ref_table))
                        )
                        temp_table <- temp_table[grep("chrX", temp_table[, 1]),]
            
                        temp_table[, 4] <- as.character(temp_table[, 4])
                        temp_table <- as.matrix(temp_table)
                        colnames(temp_table) <- colnames(sample_table)
            
                        sample_table[, 4] <- as.character(sample_table[, 4])
                        sample_table <- rbind(sample_table, temp_table)
                        sample_table[, 4] <- as.character(sample_table[, 4])
            
                    }

                }
        
            }

        }
    
        # Something is wrong here
        if (
            any(is.na(sample_table[, 2]))
            | any(is.na(sample_table[, 3]))
        ) {
            sample_table <- sample_table[
                -union(which(is.na(sample_table[, 2])),
                which(is.na(sample_table[, 3]))),
            ]
        }
    
        if (!is.vector(sample_table)) {
            sample_table <- sample_table[rowSums(!is.na(sample_table)) > 0, ]
        }
    
        sample_table <- as.data.frame(sample_table, row.names = FALSE)
        if (is.vector(sample_table)) {
            sample_table[1] <- as.character(sample_table[1])
            sample_table[2] <- as.integer(as.numeric(as.character(sample_table[2])))
            sample_table[3] <- as.integer(as.numeric(as.character(sample_table[3])))
            sample_table[4] <- as.character(sample_table[4])
            sample_table <- t(sample_table)
            sorted_sample_table <- sample_table

        } else if (ncol(sample_table) == 1) {
            sample_table <- t(sample_table)
            sample_table[, 1] <- as.character(sample_table[, 1])
            sample_table[, 2] <- as.integer(as.numeric(as.character(sample_table[, 2])))
            sample_table[, 3] <- as.integer(as.numeric(as.character(sample_table[, 3])))
            sample_table[, 4] <- as.character(sample_table[, 4])
            sorted_sample_table <- sample_table

        } else if (nrow(sample_table) > 1) {
            sample_table[, 1] <- as.character(sample_table[, 1])
            sample_table[, 2] <- as.integer(as.numeric(as.character(sample_table[, 2])))
            sample_table[, 3] <- as.integer(as.numeric(as.character(sample_table[, 3])))
            sample_table[, 4] <- as.character(sample_table[, 4])
            # Eliminate duplicates and gain/loss with same coordinates
            sorted_sample_table <- mod_merge$mergeTable(sample_table)

        }
    
        # Correct format for sample table
        # Something is  wrong here
        if (is.vector(sorted_sample_table)) {
            sorted_sample_table[1] <- as.character(sorted_sample_table[1])
            sorted_sample_table[2] <- as.integer(as.character(sorted_sample_table[2]))
            sorted_sample_table[3] <- as.integer(as.character(sorted_sample_table[3]))
            sorted_sample_table[4] <- as.character(sorted_sample_table[4])
            sorted_sample_table <- t(sorted_sample_table)

        } else if (ncol(sorted_sample_table) == 1) {
            sorted_sample_table <- t(sorted_sample_table)
            sorted_sample_table <- as.data.frame(sorted_sample_table)
            sorted_sample_table[, 1] <- as.character(sorted_sample_table[, 1])
            sorted_sample_table[, 2] <- as.integer(as.character(sorted_sample_table[, 2]))
            sorted_sample_table[, 3] <- as.integer(as.character(sorted_sample_table[, 3]))
            sorted_sample_table[, 4] <- as.character(sorted_sample_table[, 4])

        } else if (nrow(sorted_sample_table) > 1) {
            sorted_sample_table <- as.data.frame(sorted_sample_table)
            sorted_sample_table[, 1] <- as.character(sorted_sample_table[, 1])
            sorted_sample_table[, 2] <- as.integer(as.character(sorted_sample_table[, 2]))
            sorted_sample_table[, 3] <- as.integer(as.character(sorted_sample_table[, 3]))
            sorted_sample_table[, 4] <- as.character(sorted_sample_table[, 4])

        }
    
        return(list(sorted_sample_table, Dump_table, transloctable))
    
    }

}

