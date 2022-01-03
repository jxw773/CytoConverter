
##function for separating normal data
##take into acc same chrom insestion
##add cen into here (pter qter analouge)

mod_utils <- modules::use('modules/utils.R')
mod_cytobands <- modules::use('modules/cytobands.R')

colparse <- function(
        Cyto_ref_table,
        ref_table,
        coln,
        xmod,
        ymod,
        transloctable,
        addtot,
        Cyto,
        guess_q,
        constitutional,
        forMtn
) {
    Cyto_sample <- Cyto
  
    # for or statements, take first statement
    Cyto_sample[coln] <- gsub("or.*$", "", Cyto_sample[coln])
    
    # if we are guessing ? marks
    if (guess_q == T & any(grepl("\\?", Cyto_sample[coln]))) {
        Cyto_sample[coln] <- gsub("\\?", "", Cyto_sample[coln])
    }
  
    # if we are counting constitutional, substitute c
    if (constitutional == T) {
        Cyto_sample[coln] <- gsub("c$", "", Cyto_sample[coln])
    }
  
  
    # Derivative chromosomes with translocations are a loss (on native chromosome)
    #   gain (on new chromosome)
    # figure out what add bool is and consaolidate that
    # string splits by ;, takes into account multiple chromosomes odd values are chromosomes even
    #   are the positions, for derivaties, first one needs to be treated differently


    ## JP: test and temp are the same, why have both?
    test <- strsplit(
        gsub(
            "[\\(\\)]",
            "",
            regmatches(
                Cyto_sample[coln],
                gregexpr("\\(.*?\\)",  Cyto_sample[coln])
            )[[1]]
        ),
        ";"
    )
    temp <- strsplit(
        gsub(
            "[\\(\\)]",
            "",
            regmatches(
                Cyto_sample[coln],
                gregexpr("\\(.*?\\)",  Cyto_sample[coln])
            )[[1]]
        ),
        ";"
    )

    # loop to go through everything and extract table if its in short form
    coord <- data.frame()
    excoord <- data.frame()
    derMods <- Cyto_sample[coln]
    addBool <- ""

    # right now this is doing both der stuff and adition stuff, want it just to do + stuff make
    #   dermods do other stuff
    Mainchr <- vector()
    Allchr <- vector()
    
    # remeber what this is supposed to do
    derMods <- strsplit(derMods, "\\)")[[1]]
    derModsBackup <- strsplit(derMods, "\\)")[[1]]
  
    # skip second +t if first one is triggered
    plusT <- F
  
    # regins
    regtranschrom = NULL
  
    if (length(temp) > 0) {
        # chooses main chr for exclusion
        Mainchr <- temp[[1]]
        Mainchr <- paste(Mainchr, "$", sep = "")
        Mainchr <- gsub("p|q", "", Mainchr)
    }
  
    # If it straight up describes derivative makeup afterwads
    if (
        grepl(
            "der\\([0-9;]+\\)\\(|rec\\([0-9;]+\\)\\(",
            Cyto_sample[coln]
        )
        & !grepl("ider", Cyto_sample[coln])
    ) {
        addBool <- paste(addBool, "LongDer", sep = "")

    } else if (
        grepl("der|rec|dic", Cyto_sample[coln])
        & !grepl("ider|idic", Cyto_sample[coln])

    ) {
        if (grepl("der|rec", Cyto_sample[coln])) {
            addBool <- paste(addBool, derMods[1], sep = "")

        }

        # Something is going wrong where when i changed this
        #   only want more than one
        if (
            grepl(
                "der\\([[:alnum:]]+(;[[:alnum:]])*\\)[[:alpha:]]+|rec\\([[:alnum:]]+(;[[:alnum:]])*\\)[[:alpha:]]+",
                Cyto_sample[coln]
            )
        ) {
            derMods <- strsplit(
                Cyto_sample[coln],
                "(der|rec)\\([[:digit:]XY]+(;[[:digit:]XY])*\\)"
            )[[1]][2]
      
            derMods <- strsplit(derMods, "\\)")[[1]]
            temp <- utils::tail(temp, length(temp) - 1)

        } else if (
            grepl("der|rec", Cyto_sample[coln])
            && !is.na(derMods)
            && length(derMods) > 1
            && grepl("^r\\(", derMods[2])
        ) {
            derMods <- utils::tail(derMods, length(derMods) - 1)
            temp <- utils::tail(temp, length(temp) - 1)
      
        } else if (
            grepl(
                "der\\([[:digit:]]+;[[:digit:]]+\\)t\\(|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\(",
                Cyto_sample[coln]
            )
        ) {
            # handle der(11;13)t( and dic(11;13)t( esq cases here
      
            # remember to search t(11;13 if its referring to previous call of translocation
            # figure out how to do this
      
            # this will only work with two translocations and the translocation must immediantly
            #   follow the der or dic
      
            if (
                grepl(
                    "der\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)\\(|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)\\(",
                    Cyto_sample[coln]
                )
            ) {
                # translocation is fully described
                if (grepl("der\\(", Cyto_sample[coln])) {
                    # replace der with dic
                    derMods[1] <- gsub("der", "dic", derMods[1])
                    # delete der so its treated like a pure dic
                    addBool[1] <- gsub("der", "dic", addBool[1])

                }
        
                # delete translocation from both temp and derMods
                derMods <- derMods[-2]
                temp <- temp[-2]

            } else if (
                grepl(
                    "der\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)",
                    Cyto_sample[coln]
                )
            ) {
                # if the translocation needs lookup
        
                # see if translocation is found
                # else this thing crashes
                if (
                    length(transloctable) > 0
                    && grepl(
                        derMods[2],
                        names(transloctable),
                        perl = F,
                        fixed = T
                    )
                ) {
                    # last chunk of translocation indicating position
                    addOnTrans <- names(transloctable)[
                        grep(
                            derMods[2],
                            names(transloctable),
                            perl = F,
                            fixed = T
                        )
                    ]
                    addOnTrans <- strsplit(addOnTrans, "\\)")[[1]][2]
                    addOnTransTemp <- strsplit(gsub("\\(|\\)", "", addOnTrans), ";")

                    # if the string is longer than just the translocation
                    if (length(derMods) > 2) {
                        derMods <- c(derMods[1:2], addOnTrans, derMods[3:length(derMods)])
                        temp <- c(temp[1:2], addOnTransTemp, temp[3:length(temp)])

                    } else {
                        # just the der and the translocation
                        derMods <- c(derMods[1:2], addOnTrans)
                        temp <- c(temp[1:2], addOnTransTemp)

                    }
          
                    if (grepl("der\\(", Cyto_sample[coln])) {
                        # replace der with dic
                        derMods[1] <- gsub("der", "dic", derMods[1])
                        # delete der so its treated like a pure dic
                        addBool[1] <- gsub("der", "dic", addBool[1])
            
                    }
          
                    # delete translocation from both temp and derMods
                    derMods <- derMods[-2]
                    temp <- temp[-2]

                } else {
                    print("translocation undefined")
                    return(NULL)
                }
            }
        }
    }
  
  
  
    # if temp length is greater than 2 (derivative chromosome, then, dermods 1 is master indicator,
    # not including ders above)
    if (
        length(temp) > 2
        & !(grepl("der|rec", Cyto_sample[coln]))
        || grepl("ider", Cyto_sample[coln])
    ) {
        addBool <- paste(addBool, derMods[1], sep = '')

    }
  
    # if(grepl("\\+.*\\(.*",derMods[1]))
    # {
    # adds +etc on if it starts with something with +something
    # this part is being messed up
    # handle differently if we are not counting constitutional
  
    # check if there is a X3 etc value, if there is pick that up, store, add to addtot,
    # make it process through twice later
    multi = 1
    if (constitutional == F) {
        temp_cyto <- gsub("(c$)|(c\\?$)", "", Cyto_sample[coln])
        if (grepl("\\)X|\\)x", temp_cyto)) {
            multi = unlist(strsplit(temp_cyto, ")X|)x"))[2]
            if (grepl("-|~", multi)) {
                multi <- unlist(strsplit(multi, "~|-"))[1]
        
            }
            # multi=gsub("[^0-9]","",multi)
      
            multi = as.numeric(multi)
            addBool <- paste(addBool, "multi", multi, sep = "")
        }
    
    } else {
        if (grepl("\\)X|\\)x", Cyto_sample[coln])) {
            multi = unlist(strsplit(Cyto_sample[coln], ")X|)x"))[2]
            if (grepl("-|~", multi)) {
                multi <- unlist(strsplit(multi, "~|-"))[1]
        
            }
            # multi=gsub("[^0-9]","",multi)
      
            multi = as.numeric(multi)
            addBool <- paste(addBool, "multi", multi, sep = "")
        }
    }

    # increment X in presence of x modifications here
    if (any(grepl("X", Mainchr))) {
        xmod <- xmod + length(grep("X", Mainchr))
    }

    if (any(grepl("Y", Mainchr))) {
        ymod <- ymod + length(grep("Y", Mainchr))
    }
  
    # make sure you count + properly for ? marks or constitutional
    if (
        (
            any(grepl("\\?|\\~", Cyto_sample[coln]))
            & any(grepl("\\+", Cyto_sample[coln]))
        )
        | (
            constitutional == F
            & multi > 1
            & grepl("(c$)|(c\\?$)", Cyto_sample[coln])
        )
    ) {
        if (
            constitutional == F
            & multi > 2
            & !grepl("\\+", Cyto_sample[coln])
            &  grepl("(c$)|(c\\?$)", Cyto_sample[coln])
        ) {
            addtot <- addtot + (1 * (multi - 2))
      
        } else {
            addtot <- addtot + (1 * multi)

        }
    }
  
    # if guess is true, try to process ? marks
    # think about how this can affect counting  + and \\?

    # not processing things like t(9;22)(p?;q10)
  
    if (length(temp) == 1) {
        if (grepl("(t|ins)\\(", derMods[1])) {
            transchrom <- stringr::str_extract(
                Cyto_sample[coln],
                paste(
                    gsub("\\?", "\\\\?",
                        gsub("\\+", "\\\\+",
                            gsub("\\(", "\\\\(", derMods[1])
                        )
                    ),
                    "\\)\\(.+?\\)",
                    sep = ''
                )
            )

            # if this is not labled
            if (is.na(transchrom)) {
                transchrom <- stringr::str_extract(
                    Cyto_sample[coln],
                    paste(
                        gsub("\\?", "\\\\?",
                            gsub("\\+", "\\\\+",
                                gsub("\\(", "\\\\(", derMods[1])
                            )
                        ),
                        "\\)",
                        sep = ''
                    )
                )
            }
      
            regtranschrom <- gsub("\\+", "",
                gsub("\\)", "\\\\)",
                    gsub("\\(", "\\\\(", transchrom)
                )
            )
        }
    }
  
  
    if (
        (
            length(test) > 1
            | any(grepl("p|q", temp))
            | grepl(
                "(9;22)|(22;9)",
                paste(unlist(temp), collapse = ';', sep = ";")
            )
            | (
                !is.null(regtranschrom)
                && any(
                    grepl(regtranschrom, names(transloctable))
                )
            )
        )
        & !any(grepl("\\?|\\~", temp))
        & !(
            guess_q == F
            & grepl("\\?", Cyto_sample[coln])
        )
    ) {

        # goes by steps of 2, odd indexes indicate chromosomes, even indicate positions
        # length_temp<-(if((length(temp) / 2)==0.5){1}else{length(temp)/2})
        lengthcount = 1
        repeat {
            if (
                lengthcount > (
                    if (!is.integer((length(temp) / 2))) {
                        ceiling(length(temp) / 2)
                    } else {
                        length(temp) / 2
                    }
                )
            ) {
                break
            }

            if (lengthcount > 60) {
                print("while loop not terminating")
                break
            }
      
            Allchr <- c(
                Allchr, as.vector(paste(temp[[(lengthcount * 2 - 1)]], "$", sep = ""))
            )
      
            # handle those t(__;__) and ins (__;__) here with no follow up
            if (
                ((lengthcount * 2) - 1) == length(temp)
                || if (length(temp) > (lengthcount * 2 - 1)) {
                    all(grepl("^[[:digit:]]+$", temp[[lengthcount * 2]]))
                } else {
                    FALSE
                }
            ) {
                if (
                    any(grepl("[pq]", temp[[lengthcount * 2 - 1]]))
                ) {
                    if (
                        any(grepl("p", temp[[lengthcount * 2 - 1]]))
                    ) {
                        arm = "p10"
                    }

                    if (any(grepl("q", temp[[lengthcount * 2 - 1]]))) {
                        arm = "q10"
                    }
                    temp <- c(
                        if ((lengthcount * 2 - 1) != 1) {
                            temp[1:(lengthcount * 2 - 1)]
                        } else {
                            temp[lengthcount * 2 - 1]
                        },
                        arm,
                        if (length(temp) >= 2) {
                            temp[[(lengthcount * 2):length(temp)]]
                        }
                    )
                    temp[[lengthcount * 2 - 1]] <- gsub("p|q", "", temp[[lengthcount * 2 - 1]])
                }
        
                if (
                    any(grepl("^r\\(", derMods[lengthcount * 2 - 1]))
                    & forMtn == F
                ) {
                    rindex <- lengthcount * 2 - 1
          
                    rChr <- paste("chr", temp[[rindex]], "$", sep = "")
                    tempRingPosition <- NULL
                    tempRing <- NULL
          
                    for (c in 1:length(rChr)) {
                        tempRingTable <- matrix(nrow = 0, ncol = 5)
                        chr <- Cyto_ref_table[grep(rChr[c], Cyto_ref_table[, 1]), ]
                        tempRingTable <- rbind(tempRingTable, rbind(chr[1, ], chr[nrow(chr), ]))

                        # JP ??
                        tempRingTable[, 4]
            
                        tempRingPosition <- c(
                            tempRingPosition,
                            paste(tempRingTable[, 4], collapse = "")
                        )
            
                    }
          
                    tempRing <- list(unlist(temp[[rindex]]), tempRingPosition)
                    derModsRing <- c(
                        paste(
                            "r(",
                            paste(unlist(temp), sep = "", collapse = ";"),
                            sep = ""
                        ),
                        paste(
                            "(",
                            paste(tempRingPosition, collapse = ";"), sep = ""
                        )
                    )
                    tempstorage <- NULL
                    tempdermods <- NULL
          
                    if (length(temp) > (lengthcount * 2)) {
                        tempstorage <- temp[(lengthcount * 2):length(temp)]
                        tempdermods <- derMods[(lengthcount * 2):length(temp)]
                    }
          
                    temp[[lengthcount * 2 - 1]] <- tempRing[[1]]
          
                    if (length(temp) == 1) {
                        temp <- list(unlist(temp), tempRing[[2]])
                    } else {
                        temp[[lengthcount * 2]] <- tempRing[[2]]
                    }

                    temp <- temp[1:(lengthcount * 2)]
                    derMods[(lengthcount * 2) - 1] <- derModsRing[1]
          
                    if (length(derMods) == 1) {
                        derMods <- c(derMods, derModsRing[2])
                    } else {
                        derMods[lengthcount * 2] <- derModsRing[2]
                    }
          
                    derMods <- derMods[1:(lengthcount * 2)]
          
                    if (!is.null(tempstorage)) {
                        temp <- c(temp, tempstorage)
                        derMods <- c(derMods, tempdermods)
                    }

                }
        
        
                if (grepl("(t|ins)\\(", derMods[lengthcount * 2 - 1])) {
                    # Consider adding t(num;num) option
          
                    transchrom <- stringr::str_extract(
                        Cyto_sample[coln],
                        gsub("\\?", "\\\\?",
                            paste(
                                gsub("\\+", "\\\\+",
                                     gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                                ),
                                "\\)\\(.+?\\)",
                                sep = ''
                            )
                        )
                    )

                    # If this is not labled
                    if (is.na(transchrom)) {
                        transchrom <- stringr::str_extract(
                            Cyto_sample[coln],
                            gsub("\\?", "\\\\?",
                                paste(
                                    gsub("\\+", "\\\\+",
                                        gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                                    ),
                                    "\\)",
                                    sep = ''
                                )
                            )
                        )
                    }
          
                    regtranschrom <- gsub("\\+", "",
                        gsub("\\)", "\\\\)",
                            gsub("\\(", "\\\\(", transchrom)
                        )
                    )

                    temptrans <- data.frame()
          
                    # Make stuff now
                    if (
                        !any(grepl(regtranschrom, names(transloctable)))
                        & (
                            (lengthcount * 2) > length(temp)
                            || all(!grepl("p|q", temp[[lengthcount * 2]]))
                        )
                        & all(grepl("9|22", transchrom))
                    ) {
                        temptrans <- rbind(
                            temptrans,
                            cbind(
                                paste("der(", "9", ")", sep = ''),
                                paste(
                                    "der(", "9", ";", "22", ")(", "q34.1p24.3", ";", "q13.33q11", ")",
                                    sep = '',
                                    collapse = ";"
                                )
                            )
                        )

                        temptrans <- rbind(
                            temptrans,
                            cbind(
                                paste("der(", "22", ")", sep = ''),
                                paste(
                                    "der(", "9", ";", "22", ")(", "q34.3q34.1", ";", "q11.2p13", ")",
                                    sep = '',
                                    collapse = ";"
                                )
                            )
                        )
            
                        transloctable <- c(transloctable, list(temptrans))
                        names(transloctable)[length(transloctable)] <-
                            gsub("\\+", "", transchrom)
                    }
          
                    if (
                        grepl("^\\++t\\(", derMods[lengthcount * 2 - 1])
                        & !grepl("der", Cyto_sample[coln])
                    ) {
                        plusT = T
            
                        transderplus <- transloctable[grep(regtranschrom, names(transloctable))]
                        dertransextract <- transderplus[[1]][
                            grep(
                                paste(
                                    "der\\(",
                                    temp[[(lengthcount * 2) - 1]],
                                    "\\)",
                                    sep = "",
                                    collapse = "|"
                                ),
                                transderplus[[1]][, 1]
                            ),
                            2
                        ]
            
                        transtemp <- unlist(
                            lapply(
                                as.list(dertransextract),
                                function(x) {
                                    strsplit(
                                        gsub("[\\(\\)]", "",
                                            regmatches(
                                                x,
                                                gregexpr("\\(.*?\\)",  x)
                                            )[[1]]
                                        ),
                                        ";"
                                    )
                                }
                            ),
                            recursive = F
                        )
            
                        transdermod <- unlist(
                            lapply(
                                dertransextract,
                                function(x) {
                                    y <- paste("+", x, sep = "")
                                    strsplit(y, "\\)")
                                }
                            )
                        )

                        # Make temp storage for rest of temp
                        tempstorage <- NULL
                        tempdermods <- NULL
                        if (length(temp) > (lengthcount * 2)) {
                            tempstorage <- temp[[((lengthcount * 2) + 1):length(temp)]]
                            tempdermods <- derMods[((lengthcount * 2) + 1):length(temp)]
                        }
            
                        temp[[lengthcount * 2 - 1]] <- transtemp[[1]]
            
                        if (length(temp) == 1) {
                            temp <- list(unlist(temp), transtemp[[2]])
                        } else {
                            temp[[lengthcount * 2]] <- transtemp[[2]]
                        }

                        temp <- c(temp[1:(lengthcount * 2)], transtemp[3:length(transtemp)])
            
                        derMods[(lengthcount * 2) - 1] <- transdermod[1]
                        if (length(derMods) == 1) {
                            derMods <- c(derMods, transdermod[2])
                        } else {
                            derMods[lengthcount * 2] <- transdermod[2]
                        }
            
                        derMods <- c(
                            derMods[1:(lengthcount * 2)],
                            transdermod[3:length(transdermod)]
                        )
            
                        if (!is.null(tempstorage)) {
                            temp <- c(temp, tempstorage)
                            derMods <- c(derMods, tempdermods)
                        }

                    } else {
                        transchrom <- stringr::str_extract(
                            Cyto_sample[coln],
                            paste(
                                gsub("\\?", "\\\\?",
                                    gsub("\\+", "\\\\+",
                                        gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                                    )
                                ),
                                "\\)",
                                sep = ''
                            )
                        )

                        selectedTransTable <- as.matrix(
                            transloctable[
                                grep(
                                    gsub("\\)", "\\\\)",
                                        gsub("\\(", "\\\\(", transchrom)
                                    ),
                                    names(transloctable)
                                )
                            ][[1]]
                        )

                        transDer <- selectedTransTable[
                            grep(
                                paste(
                                    "der\\(",
                                    paste(
                                        sapply(
                                            Mainchr,
                                            function(x) {
                                                substr(x, 0, nchar(x) - 1)
                                            }
                                        ),
                                        sep = '',
                                        collapse = ';'
                                    ),
                                    "\\)",
                                    sep = ''
                                ),
                                selectedTransTable[, 1]
                            ), 
                        ][2]

                        transtemp <- strsplit(
                            gsub("[\\(\\)]", "",
                                regmatches(transDer, gregexpr("\\(.*?\\)",  transDer))[[1]]
                            ),
                            ";"
                        )

                        # This isnt working
                        temp <- c(
                            if ((lengthcount * 2 - 1) != 1) {
                                temp[1:(lengthcount * 2 - 1)]
                            } else {
                                transtemp[1]
                            },
                            transtemp[2],
                            if (
                                (length(temp) >= lengthcount * 2)
                                && (
                                    !grepl("p|q|->|:", derMods[length(temp)])
                                    || IsOdd(length(temp))
                                )
                            ) {
                                temp[lengthcount * 2:length(temp)]
                            }
                        )
                    }

                    # Instead of below, for now, replace with equivilent derivative chromosome
          
                    # Get stuff in transloctable that matches same chromosom
                    # transtable <-
                    #  cbind(transloctable[grep(paste("chr", Mainchr, sep = '', collapse = '|'),
                    #        transloctable[, 1]), 1:3], derMods[(lengthcount * 2-1)])
                    # coord <- rbind(coord, transtable)
                    # temp <- temp[[-((lengthcount * 2-1))]]
                    # derMods <- derMods[-(lengthcount)]
                    # lengthcount <- lengthcount - 1
                }
            }
      
            if (any(grepl("[pq]", temp[[lengthcount * 2 - 1]]))) {
                # get chromosome, at 10 to other half, get coordinates for that
                temp[[lengthcount * 2 - 1]] <- gsub("p|q", "", temp[[lengthcount * 2 - 1]])
            }

            if (grepl("t\\(", derMods[lengthcount * 2 - 1])) {
                # if translocation, do this, shouldnt be placed, should be placed one loop above
                # if not in the table, add to table
                # names are just t(1;12)
        
                # if its t(9;22) or t(22;9) , add to table if not already found
        
                transchrom <- stringr::str_extract(
                    Cyto_sample[coln],
                    paste(
                        gsub("\\?", "\\\\?",
                            gsub("\\+", "\\\\+",
                                gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                            )
                        ),
                        "\\)\\(.+?\\)",
                        sep = ''
                    )
                )

                # if this is not labled
                if (is.na(transchrom)) {
                    transchrom <- stringr::str_extract(
                        Cyto_sample[coln],
                        paste(
                            gsub("\\?", "\\\\?",
                                gsub("\\+", "\\\\+",
                                    gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                                )
                            ),
                            "\\)",
                            sep = ''
                        )
                    )
                }

                regtranschrom <- gsub("\\+", "",
                    gsub("\\)", "\\\\)",
                        gsub("\\(", "\\\\(", transchrom)
                    )
                )

                if (!any(grepl(regtranschrom, names(transloctable)))) {

                    # entire translocation only if nessesary
                    # stringr::str_extract(mem, "t\\([[:digit:]]+(;[[:digit:]])*?\\)\\(.+?\\)")

                    temptrans <- data.frame()
          
                    # Make stuff now

                    if (
                        !any(grepl(regtranschrom, names(transloctable)))
                        & (
                            (lengthcount * 2) > length(temp)
                            || all(!grepl("p|q", temp[[(lengthcount * 2)]]))
                        )
                        & all(grepl("9|22", transchrom))
                    ) {
                        temptrans <- rbind(
                            temptrans,
                            cbind(
                                paste("der(", "9", ")", sep = ''),
                                paste(
                                    "der(", "9", ";", "22", ")(", "q34.1p24.3", ";", "q13.33q11", ")",
                                    sep = '',
                                    collapse = ";"
                                )
                            )
                        )
            
                        temptrans <- rbind(
                            temptrans,
                            cbind(
                                paste("der(", "22", ")", sep = ''),
                                paste(
                                    "der(", "9", ";", "22", ")(", "q34.3q34.1", ";", "q11.2p13", ")",
                                    sep = '',
                                    collapse = ";"
                                )
                            )
                        )

                    } else {
                        for (o in 1:length(temp[[lengthcount * 2 - 1]])) {
                            tempcurvec <- mod_cytobands$getCytoBands(
                                Cyto_ref_table,
                                Cyto_sample,
                                lengthcount,
                                o,
                                temp,
                                coln,
                                derMods,
                                forMtn
                            )

                            currentvec <- tempcurvec[[1]]
                            earlyReturn <- tempcurvec[[2]]
              
                            if (earlyReturn == T) {
                                return(NA)

                            } else {
                                # vector of o cytobands # includes everything but cytoband on der o indicated,
                                # double check how your doing this (going one band before, is this right?)

                                # above is all to find vectors of o that are excluded
                                if (o == 1) {
                                    temptrans <- rbind(
                                        temptrans,
                                        cbind(
                                            paste("der(", temp[[lengthcount * 2 - 1]][o], ")", sep = ''),
                                            paste(
                                                "der(",
                                                temp[[lengthcount * 2 - 1]][o],
                                                ";",
                                                temp[[lengthcount * 2 - 1]][length(temp[[lengthcount * 2 - 1]])],
                                                if (length(currentvec) > 1) {
                                                    ";"
                                                },
                                                rep(temp[[lengthcount * 2 - 1]][o], length(currentvec) - 1),
                                                ")(",
                                                currentvec[1],
                                                ";",
                                                temp[[lengthcount * 2]][length(temp[[lengthcount * 2 - 1]])],
                                                if (length(currentvec) > 1) {
                                                    ";"
                                                },
                                                currentvec[-1],
                                                ")",
                                                sep = '',
                                                collapse = ";"
                                            ),
                                            ""
                                        )
                                    )

                                } else{
                                    temptrans <- rbind(
                                        temptrans,
                                        cbind(
                                            paste("der(", temp[[(lengthcount * 2 - 1)]][o], ")", sep = ''),
                                            paste(
                                                "der(",
                                                temp[[lengthcount * 2 - 1]][o],
                                                ";",
                                                temp[[lengthcount * 2 - 1]][o - 1],
                                                if (length(currentvec) > 1) {
                                                    ";"
                                                },
                                                rep(temp[[lengthcount * 2 - 1]][o], length(currentvec) - 1),
                                                ")(",
                                                currentvec[1],
                                                ";",
                                                temp[[lengthcount * 2]][o - 1],
                                                if (length(currentvec) > 1) {
                                                    ";"
                                                },
                                                currentvec[-1],
                                                ")",
                                                sep = '',
                                                collapse = ";"
                                            ),
                                            ""
                                        )
                                    )
                  
                                }

                            }

                        }

                    }

                    transloctable <- c(transloctable, list(temptrans))
                    names(transloctable)[length(transloctable)] <- gsub("\\+", "", transchrom)

                }

                # For extending der(X)t(X;y) type
                if (
                    !grepl("^\\++t\\(", derMods[lengthcount * 2 - 1])
                    & grepl("der", Cyto_sample[coln])
                ) {

                    # this is a factor, not a table
                    transDer <- as.matrix(
                        transloctable[[grep(regtranschrom, names(transloctable))]][
                            grep(
                                paste(
                                    "der\\(",
                                    sapply(
                                        Mainchr,
                                        function(x) {
                                            substr(x, 0, nchar(x) - 1)
                                        }
                                    ),
                                    "\\)",
                                    sep = "",
                                    collapse = "|"
                                ),
                                transloctable[[grep(regtranschrom, names(transloctable))]][, 1]
                            ),
                        ]
                    )[2]

                    # If it is NA, get chromosome from previous translocation and use that instead
                    if (is.na(transDer) & lengthcount > 1) {
                        transDer <- as.matrix(
                            transloctable[[grep(regtranschrom, names(transloctable))]][
                                grep(
                                    paste(
                                        "der\\(",
                                        "(",
                                        paste(temp[[(lengthcount - 1) * 2 - 1]], collapse = "|"),
                                        ")",
                                        "\\)",
                                        sep = "",
                                        collapse = "|"
                                    ),
                                    transloctable[[grep(regtranschrom, names(transloctable))]][, 1]
                                ),
                            ]
                        )[2]
                    }

                    # Translate derivative into temp stuff
                    transtemp <- strsplit(
                        gsub("[\\(\\)]", "",
                            regmatches(transDer, gregexpr("\\(.*?\\)",  transDer))[[1]]
                        ),
                        ";"
                    )
                    temp[[lengthcount * 2 - 1]] <- transtemp[[1]]
                    temp[[lengthcount * 2]] <- transtemp[[2]]

                }
        
                # fix +t() here
                if (
                    grepl("^\\++t\\(", derMods[lengthcount * 2 - 1])
                    & !grepl("der", Cyto_sample[coln])
                    & plusT == F
                ) {
                    # print(lengthcount)
                    transderplus <- transloctable[grep(regtranschrom, names(transloctable))]
                    dertransextract <- transderplus[[1]][
                        grep(
                            paste(
                                "der\\(",
                                temp[[(lengthcount * 2) - 1]],
                                "\\)",
                                sep = "",
                                collapse = "|"
                            ),
                            transderplus[[1]][, 1]
                        ),
                        2
                    ]
          
                    transtemp <- unlist(
                        lapply(
                            as.list(dertransextract),
                            function(x) {
                                strsplit(
                                    gsub("[\\(\\)]", "",
                                        regmatches(x, gregexpr("\\(.*?\\)",  x))[[1]]
                                    ),
                                    ";"
                                )
                            }
                        ),
                        recursive = F
                    )
          
                    transdermod <- unlist(
                        lapply(
                            dertransextract,
                            function(x) {
                                y <- paste("+", x, sep = "")
                                strsplit(y, "\\)")
                            }
                        )
                    )
          
                    # Make temp storage for rest of temp
                    tempstorage <- NULL
                    tempdermods <- NULL
                    if (length(temp) > (lengthcount * 2)) {
                        tempstorage <- temp[[((lengthcount * 2) + 1):length(temp)]]
                        tempdermods <- derMods[((lengthcount * 2) + 1):length(temp)]
                    }
          
                    temp[[lengthcount * 2 - 1]] <- transtemp[[1]]
          
                    if (length(temp) == 1) {
                        temp <- list(unlist(temp), transtemp[[2]])

                    } else{
                        temp[[lengthcount * 2]] <- transtemp[[2]]
            
                    }

                    temp <- c(temp[1:(lengthcount * 2)], transtemp[3:length(transtemp)])
          
                    derMods[(lengthcount * 2) - 1] <- transdermod[1]
                    if (length(derMods) == 1) {
                        derMods <- c(derMods, transdermod[2])

                    } else {
                        derMods[lengthcount * 2] <- transdermod[2]

                    }

                    derMods <- c(derMods[1:(lengthcount * 2)], transdermod[3:length(transdermod)])
          
                    if (!is.null(tempstorage)) {
                        temp <- c(temp, tempstorage)
                        derMods <- c(derMods, tempdermods)
            
                    }

                }

            }
      
            if (grepl("ins\\(", derMods[lengthcount * 2 - 1])) {
                currentvec <- vector()
                inschrom <- stringr::str_extract(
                    Cyto_sample[coln],
                    paste(
                        gsub("\\?", "\\\\?",
                            gsub("\\+", "\\\\+",
                                gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                            )
                        ),
                        "\\)\\(.+?\\)",
                        sep = ''
                    )
                )
        
                if (is.na(inschrom)) {
                    inschrom <- stringr::str_extract(Cyto_sample[coln],
                        paste(
                            gsub("\\?", "\\\\?",
                                gsub("\\+", "\\\\+",
                                    gsub("\\(", "\\\\(", derMods[lengthcount * 2 - 1])
                                )
                            ),
                            "\\)",
                            sep = ''
                        )
                    )
                }

                reginschrom <- gsub("\\+", "",
                    gsub("\\)", "\\\\)",
                        gsub("\\(", "\\\\(", inschrom)
                    )
                )

                # If insertion is on a single chromosome , do this , if for some reason the call is
                #   separated but the chromosome is not
                if (length(temp[[lengthcount * 2 - 1]]) < length(temp[[lengthcount * 2]])) {
                    temp[[lengthcount * 2 - 1]] <- rep(
                        temp[[lengthcount * 2 - 1]],
                        length(temp[[lengthcount * 2]])
                    )
                }
        
                # Think about this harder
                # If insertion is one whole sequence withint a chromosome
                # If within is true, dont make another derivative chromosome
                within = F
                if (
                    length(temp[[lengthcount * 2 - 1]]) == 1
                    && any(
                        length(temp[[lengthcount * 2 - 1]]) <
                            stringr::str_count(temp[[lengthcount * 2]], "p|q")
                    )
                ) {
                    temp[[lengthcount * 2 - 1]] <- rep(temp[[lengthcount * 2 - 1]], 2)
          
                    # location to split string (2nd q or p)
                    splitloc <- str_locate_all(temp[[lengthcount * 2]], "p|q")[[1]][, 1][2]
          
                    temp[[lengthcount * 2]] <- c(
                        str_sub(temp[[lengthcount * 2]], 1, splitloc - 1),
                        str_sub(
                            temp[[lengthcount * 2]],
                            splitloc,
                            stringr::str_length(temp[[lengthcount * 2]])
                        )
                    )
                    within = T
                    
                }
        
                if (!any(grepl(reginschrom, names(transloctable)))) {
                    # entire translocation only if nessesary
                    tempins <- data.frame()
          
                    currentvec <- mod_cytobands$getCytoBands(
                        Cyto_ref_table,
                        Cyto_sample,
                        lengthcount,
                        1,
                        temp,
                        coln,
                        derMods,
                        forMtn
                    )[[1]]
          
                    # fix this later , insertion is simply stopping where it left of
                    # instead of listing other half of the chromosome
          
                    tempins <- rbind(
                        tempins,
                        cbind(
                            paste("der(", temp[[(lengthcount * 2 - 1)]][1], ")", sep = ''),
                            paste(
                                "der(",
                                temp[[lengthcount * 2 - 1]][1],
                                ";",
                                temp[[lengthcount * 2 - 1]][2],
                                ";",
                                temp[[lengthcount * 2 - 1]][1],
                                ")(",
                                temp[[lengthcount * 2]][1],
                                ";",
                                temp[[lengthcount * 2]][2],
                                ";",
                                currentvec[1],
                                ")",
                                sep = '',
                                collapse = ";"
                            ),
                            ""
                        )
                    )
          
                    # change this to be a del chromosome
                    currentvec <- mod_cytobands$getCytoBands(
                        Cyto_ref_table,
                        Cyto_sample,
                        lengthcount,
                        2,
                        temp,
                        coln,
                        derMods,
                        forMtn
                    )[[1]]
                
                    tempins <- rbind(
                        tempins,
                        cbind(
                            paste("der(", temp[[lengthcount * 2 - 1]][2], ")", sep = ''),
                            paste(
                                "der(",
                                temp[[(lengthcount * 2 - 1)]][2],
                                ";",
                                temp[[(lengthcount * 2 - 1)]][2],
                                ")(",
                                currentvec[1],
                                ";",
                                currentvec[2],
                                ")",
                                sep = '',
                                collapse = ";"
                            ),
                            ""
                        )
                    )
          
                    transloctable <- c(transloctable, list(tempins))
                    names(transloctable)[length(transloctable)] <- inschrom
          
                    # dont take into account second derivative chromosome if there
                    # is only 1 chromosome where the insersion occures
                    if (within) {
                        transloctable[[length(transloctable)]] <-
                            transloctable[[length(transloctable)]][-2, ]
                    }
          
                }
        
                if (!grepl("^\\+*ins\\(", Cyto_sample[coln]))  {
                    # this is a factor, not a table
                    insDer <- as.matrix(
                        transloctable[[grep(reginschrom, names(transloctable))]][
                            grep(
                                paste(
                                    "der\\(",
                                    substr(Mainchr, 0, nchar(Mainchr) - 1),
                                    "\\)",
                                    sep = ""
                                ),
                                transloctable[[grep(reginschrom, names(transloctable))]][, 1]
                            ),
                        ]
                    )[2]

                    instemp <-  strsplit(
                        gsub(
                            "[\\(\\)]",
                            "",
                            regmatches(insDer, gregexpr("\\(.*?\\)",  insDer))[[1]]
                        ),
                        ";"
                    )

                    temp[[lengthcount * 2 - 1]] <- instemp[[1]]
                    temp[[lengthcount * 2]] <- instemp[[2]]
                }
            }

            # check for labling long der again
            if (
                grepl(
                    "der\\([0-9;]+\\)\\(|rec\\([0-9;]+\\)\\(",
                    derMods[lengthcount * 2 - 1]
                )
                & !grepl("ider", Cyto_sample[coln])
            ) {
                addBool <- paste(addBool, "LongDer", sep = "")
            }

            for (i in 1:length(temp[[lengthcount * 2 - 1]])) {
                if (length(test) > 1 | plusT) {
                    if (!is.na(stringr::str_length(temp[[lengthcount * 2]][i]))) {
                        chr_table <- Cyto_ref_table[
                            grep(
                                paste(
                                    paste("chr", temp[[lengthcount * 2 - 1]][i], sep = ""),
                                    "$",
                                    sep = ""
                                ),
                                Cyto_ref_table
                            ),
                        ]

                        # put handling ? and ~ here
                        # put satilites and stuff here
                        # if long form
                        # deal with centromere stuff grep acen, take 2nd last
                        # der long form handling
            
                        # convert all - to ~ and cut out outer part of ~ if it is between 2 numbers
                        # need to do this better
                        if (grepl("~|(- &!(->))", temp[[lengthcount * 2]][i])) {
                            # maybe add function to make this less conservative
                            # right now it takes earlier one
                            # make this reg expression exclude ->
                            temp[[lengthcount * 2]][i] <- gsub("-", "~", temp[[lengthcount * 2]][i])
                            if (grepl("[pq][[:digit:]]+~", temp[[lengthcount * 2]][i])) {
                                # make sure this handles before and end
                                temp[[lengthcount * 2]][i] <- unlist(
                                    strsplit(temp[[lengthcount * 2]][i], "~")
                                )[1]
                            }
                        }

                        # if long form
                        if (
                            any(
                                grepl("::", temp[[lengthcount * 2]][i])
                                | grepl("~>", temp[[lengthcount * 2]][i])
                                | grepl("->", temp[[lengthcount * 2]][i])
                            )
                        ) {

                            # find p or q, take stuff before take stuff after, before is chromosomes
                                # after is positions, this will return 2 objects, must take into account
                            # parse data according to ::, in front of p and q are chromosomes, 
                                # if qter or pter, do stuff, afterward is position, make table
                                # of things included, then make list of stuff excluded
                            # ask tom about this one
                            # only splits first one
                            longform_table <- strsplit(
                                strsplit(temp[[lengthcount * 2]][i], "::")[[1]],
                                "(~>)|(->)"
                            )

                            # take away any front loaded :
                            longform_table <- lapply(
                                longform_table,
                                function(x) {
                                    gsub(':', '', x)
                                }
                            )
                
                            in_table = data.frame()
                
                            # mark for dic, trc
                            if (grepl("dic|trc", derMods[lengthcount])) {
                                addBool <- paste("long", addBool, sep = '')
                            }
            
                            # get data for each read of something->something
                            for (j in 1:length(longform_table)) {
                                stringdx <- str_locate_all(pattern = "p|q", longform_table[[j]])
                                chr_name_long <- substr(longform_table[[j]][1], 0, stringdx[[1]][1] - 1)
                                positions <- as.vector(
                                    cbind(
                                        substr(
                                            longform_table[[j]][1],
                                            stringdx[[1]][1],
                                            nchar(longform_table[[j]][1])
                                        ),
                                        substr(
                                            longform_table[[j]][2],
                                            stringdx[[2]][1],
                                            nchar(longform_table[[j]][2])
                                        )
                                    )
                                )
                    
                                if (nchar(chr_name_long) != 0) {
                                    chr_table_2 <- Cyto_ref_table[
                                        grep(
                                            paste(
                                                paste("chr", chr_name_long, sep = ""),
                                                "$",
                                                sep = ""
                                            ),
                                            Cyto_ref_table
                                        ),
                                    ]
                                } else {
                                    chr_table_2 <- chr_table
                                }
                    
                                # account for terminal ends
                                if (any(grepl("pter", positions))) {
                                    positions[grep("pter", positions)] <- chr_table_2[1, 4]
                                }
                                if (any(grepl("qter", positions))) {
                                    positions[grep("qter", positions)] <- chr_table_2[nrow(chr_table_2), 4]
                                }
                                # be careful on centromeneter
                                if (any(grepl("cen", positions))) {
                                    # likey will have to be careful about this one
                                    positions[grep("cen", positions)] <- chr_table_2[
                                        grep("acen", chr_table_2[, 5]),
                                    ][1, 4]
                                }
                    
                                # account for p10/q10
                                # double check naming convention for this one
                                # have to change the one for q
                                positions[grep("p10", positions)] <-
                                    chr_table_2[grep("acen", chr_table_2[, 5]),][1, 4]
                                positions[grep("q10", positions)] <-
                                    chr_table_2[grep("acen", chr_table_2[, 5]),][2, 4]
                    
                                # make sure positions is in order to be processed correctly
                                positions <- mod_utils$positionSorter(positions)
                                if (length(unlist(strsplit(positions, "-|~"))) > 1) {
                                    positions <- unlist(
                                        lapply(
                                            strsplit(positions, "-|~"),
                                            function(x) {
                                                x[1]
                                            }
                                        )
                                    )
                                }

                                # put stuff in table
                                positions_table <- matrix(
                                    chr_table_2[
                                        grep(
                                            paste(positions, collapse = "|", sep = "|"),
                                            chr_table_2[, 4]
                                        ),
                                    ],
                                    ncol = 5
                                )
                    
                                if (
                                    is.vector(positions_table)
                                    | nrow(positions_table) == 1
                                    | ncol(positions_table) == 1
                                ) {
                                    if (ncol(positions_table) == 1) {
                                        positions_table <- t(positions_table)
                                    }
                                    in_table <- rbind(
                                        in_table,
                                        cbind(
                                            chr_table_2[1, 1],
                                            positions_table[2],
                                            positions_table[3],
                                            derMods[lengthcount * 2 - 1]
                                        )
                                    )
                                } else {
                                    in_table <- rbind(
                                        in_table,
                                        cbind(
                                            chr_table_2[1, 1],
                                            positions_table[1, 2],
                                            positions_table[nrow(positions_table), 3],
                                            derMods[lengthcount * 2 - 1]
                                        )
                                    )
                                }
                            } # for (j in 1:length(longform_table))

                            # add to coord
                    
                            # combine piecewise with total count
                            coord <- rbind(coord, in_table)
                            # make note where the breaks are
                            # excoord<-rbind(excoord,) ##constant exclusion
            
                        } else {
                            in_table <- data.frame()
                            positions <- strsplit(
                                gsub("q", ",q",
                                    gsub("p", ",p", temp[[lengthcount * 2]][i])
                                ),
                                ","
                            )[[1]][2:length(
                                strsplit(
                                    gsub("q", ",q",
                                        gsub("p", ",p", temp[[lengthcount * 2]][i])
                                    ),
                                    ","
                                )[[1]]
                            )]

                            # have to change the one to q
                            positions[grep("p10", positions)] <-
                                chr_table[grep("acen", chr_table[, 5]),][1, 4]
                            positions[grep("q10", positions)] <-
                                chr_table[grep("acen", chr_table[, 5]),][2, 4]
                            # check order, test this
                            postions <- mod_utils$positionSorter(positions)
                            if (length(unlist(strsplit(positions, "-|~"))) > 1) {
                                positions <- unlist(
                                    lapply(
                                        strsplit(positions, "-|~"),
                                        function(x) {
                                            x[1]
                                        }
                                    )
                                )
                            }
                    
                            positions_table <- matrix(
                                chr_table[grep(paste(positions, collapse = "|"), chr_table[, 4]),],
                                ncol = 5
                            )
                    
                            # if only one q or p , add on end point
                            # fix this so if p or q doesnt rely on order (p always being before q)
                            if (length(positions) == 1) {
                                if (any(grepl("p", positions))) {
                                    if (is.vector(positions_table)) {
                                        in_table <- rbind(
                                            in_table,
                                            cbind(
                                                chr_table[1, 1],
                                                "0",
                                                positions_table[3],
                                                derMods[lengthcount * 2 - 1]
                                            )
                                        )
                                    } else {
                                        in_table <- rbind(
                                            in_table,
                                            cbind(
                                                chr_table[1, 1],
                                                "0",
                                                positions_table[nrow(positions_table), 3],
                                                derMods[lengthcount * 2 - 1]
                                            )
                                        )
                                    }
                                }

                                if (any(grepl("q", positions))) {
                                    if (is.vector(positions_table)) {
                                        in_table <- rbind(
                                            in_table,
                                            cbind(
                                                chr_table[1, 1],
                                                positions_table[2],
                                                ref_table[
                                                    grep(
                                                        paste(chr_table[1, 1], "$", sep = ""), ref_table
                                                    ),
                                                ][2],
                                                derMods[lengthcount * 2 - 1]
                                            )
                                        )
                                    } else {
                                        in_table <- rbind(
                                            in_table,
                                            cbind(
                                                chr_table[1, 1],
                                                positions_table[1, 2],
                                                ref_table[
                                                    grep(
                                                        paste(chr_table[1, 1], "$", sep = ""), ref_table
                                                    ),
                                                ][2],
                                                derMods[lengthcount * 2 - 1]
                                            )
                                        )
                                    }
                                }
                            } else {
                                in_table <- rbind(
                                    in_table,
                                    cbind(
                                        chr_table[1, 1],
                                        positions_table[1, 2],
                                        positions_table[nrow(positions_table), 3],
                                        derMods[lengthcount * 2 - 1]
                                    )
                                )
                            }

                            # combine piecewise with total count
                            coord <- rbind(coord, in_table)
                            # make note where the breaks are
                            # excoord<-rbind(excoord,) ##constant exclusion
                        }
                    } # if
                } # if
            } # for (i in 1:length(temp[[lengthcount * 2 - 1]]))

            lengthcount <- lengthcount + 1
      
        } # repeat
    }

    # for loop every chr
    # things are getting offset here
    # +der ider duplicate
    if (
        grepl("ider", Cyto_sample[coln])
        & grepl("t\\(", Cyto_sample[coln])
    ) {
        coord <- coord[-1, ]
    }
  
    if (nrow(coord) > 0) {
        # case for vectors
        coord[, 1] <- as.character(coord[, 1])
    
        # convert to numeric
        coord[, 2:3] <- apply(
            coord[, 2:3],
            2,
            function(x) {
                as.numeric(as.character(x))
            }
        )
    
        # sort in order numerically
        # coord[,2:3]<-t(apply(coord[,2:3],1,function(x){sort(x)}))
    
        # if there are two translocations, fix coordinate overlap so only overlap counts,
    
        # if(str_count(Cyto_sample[coln],"t\\(")>1)
        # {
        ##find which chromosomes have more than one reading
        ##index of entries with chromosomes
    
        #add to coord
    
        ##combine piecewise with total count
        ##coord <- rbind(coord, in_table)
        ##make note where the breaks are
        ##excoord<-rbind(excoord,) ##constant exclusion
        ##chromindex<-lapply(unique(coord[,1]),function(x){grep(x,coord[,1])})
        ##which index in chrom index have more than 2 readings
        ##relevantindex<-which(lapply(chromindex,length)>1)
        ##for(z in 1:length(relevantindex)) 
        ##{
        ##  tempcoord<-coord[chromindex[[relevantindex[z]]],]
        ##    }
    
        ##}
    
        # make sure coords are in numerical order
        coord <- coord[order(coord[, 1], coord[, 2], coord[, 3]), ]
    
    
    
        # all chromosomes used
        Chr <- unique(c(Mainchr, Allchr))
    
    
        # make sure translocations are consistant
        for (p in 1:length(Chr)) {
            curChr <- Chr[p]
      
            tempChr = coord[
                grep(
                    paste(
                        "chr",
                        curChr,
                        "$",
                        sep = "",
                        collapse = "|"
                    ),
                    coord[, 1]
                ),
            ]
            tempChr = apply(tempChr, 2, as.character)

            if (
                !is.vector(tempChr)
                && grepl("t\\(", tempChr[, 4])
                && nchar(addBool) > 0
            ) {
                # sub out translocations--delete areas of no overlap, mainchr only,
                    # translocatios only
                
                transloc <- tempChr[grep("t\\(", tempChr[, 4]), ]
                if (
                    (
                        is.data.frame(transloc)
                        | is.matrix(transloc)
                    )
                    && nrow(transloc) > 1
                ) {
                    if (
                        DescTools::Overlap(
                            as.numeric(transloc[1, 2:3]),
                            as.numeric(transloc[2, 2:3])
                        )
                    ) {
                        # if nothing , handle
                        if (
                            as.numeric(transloc[1, 3]) == as.numeric(transloc[2, 2])
                        ) {

                            temptransloc <- c(
                                tempChr[1, 1],
                                as.numeric(transloc[1, 2]),
                                as.numeric(transloc[2, 3]),
                                tempChr[1, 4]
                            )
                            tempChr <- tempChr[-1 * grep("t\\(", tempChr[, 4]), ]

                            tempChr <- rbind(tempChr, temptransloc)
                            coord <- rbind(
                                temptransloc,
                                coord[
                                    -1 * intersect(
                                        grep("t\\(", coord[, 4]),
                                        grep(paste("chr", curChr, sep = ""), coord[, 1])
                                    ), 
                                ]
                            )
                        } else {
                            temptransloc <- c(
                                tempChr[1, 1],
                                mod_utils$mergeIntOverlap(
                                    as.numeric(transloc[1, 2:3]),
                                    as.numeric(transloc[2, 2:3])
                                ),
                                tempChr[1, 4]
                            )
                            tempChr <- tempChr[-1 * grep("t\\(", tempChr[, 4]), ]
                            tempChr <- rbind(tempChr, temptransloc)
                            coord <- rbind(
                                temptransloc,
                                coord[
                                    -1 * intersect(
                                        grep("t\\(", coord[, 4]),
                                        grep(paste("chr", curChr, sep = ""), coord[, 1])
                                    ), 
                                ]
                            )
                        }
                    }
                }
            }

            if (is.data.frame(tempChr) && nrow(tempChr) == 1) {
                tempChr <- as.vector(tempChr)
            }
        } # for (p in 1:length(Chr))
    
    
        # make sure coords are in numerical order again
        coord <- coord[order(coord[, 1], coord[, 2], coord[, 3]), ]
    
    
        # make excoord table
    
        for (p in 1:length(Mainchr)) {
            # need something for $ so chr 1 doesnt net chr 10
            curChr <- Mainchr[p]
      
            tempChr = coord[
                grep(
                    paste(
                        "chr",
                        curChr,
                        "$",
                        sep = "",
                        collapse = "|"
                    ),
                    coord[, 1]
                ),
            ]
            tempChr = apply(tempChr, 2, as.character)
      
            if (length(tempChr) == 0) {

                # go back to here
                excoord = rbind(
                    excoord,
                    cbind(
                        paste("chr", gsub("\\$", "", curChr), sep = ""),
                        0,
                        ref_table[grep(paste("chr", curChr, sep = ""), ref_table),][2],
                        derModsBackup[1]
                    )
                )
        
            } else if (is.vector(tempChr)) {
                if (
                    !identical(gsub(" ", "", as.character(tempChr[2])), "0")
                ) {
                    excoord = rbind(
                        excoord,
                        cbind(tempChr[1], 0, tempChr[2], tempChr[4])
                    )
                }
        
                if (
                    !identical(
                        gsub(" ", "", as.character(tempChr[3])),
                        gsub(" ", "",
                            as.character(
                                ref_table[
                                    grep(paste(tempChr[1], "$", sep = ""), ref_table),
                                ][2]
                            )
                        )
                    )
                ) {
                    excoord = rbind(
                        excoord,
                        cbind(
                            tempChr[1],
                            tempChr[3],
                            ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table),][2], tempChr[4]
                        )
                    )
                }
            } else {
                # make sure none are inside intervals, if it is, delete out of tempChr
                for (z in 1:nrow(tempChr)) {
                    inbetween <- F
                    for (x in 1:(nrow(tempChr))) {
                        if (
                            as.numeric(tempChr[x, 2]) < as.numeric(tempChr[z, 2])
                            & as.numeric(tempChr[x, 3] > tempChr[z, 3])
                        ) {
                            inbetween <- T
                        }
                    }

                    if (inbetween == T) {
                        tempChr[z, ] <- c(0, 0, 0, 0)
                    }
                }
        
                tempChr <- tempChr[which(tempChr[, 1] != "0"), ]
        
        
                # if its a vector, do top part again
                if (is.vector(tempChr)) {
                    if (!identical(gsub(" ", "", as.character(tempChr[2])), "0")) {
                        excoord = rbind(
                            excoord,
                            cbind(tempChr[1], 0, tempChr[2], tempChr[4])
                        )
                    }
          
                    if (
                        !identical(
                            gsub(" ", "", as.character(tempChr[3])),
                            gsub(" ", "", as.character(
                                ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table),][2]
                            ))
                        )
                    ) {
                        excoord = rbind(
                            excoord,
                            cbind(
                                tempChr[1],
                                tempChr[3],
                                ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table),][2], tempChr[4]
                            )
                        )
                    }
                } else {
                    for (z in 1:nrow(tempChr)) {
                        if (z == 1) {
                            if (
                                !identical(
                                    gsub(" ", "", as.character(tempChr[1, 2])),
                                    "0"
                                )
                            ) {
                                excoord = rbind(
                                    excoord,
                                    cbind(tempChr[z, 1], 0, tempChr[z, 2], tempChr[z, 4])
                                )
                            }
              
                            if (
                                nrow(tempChr) > 1
                                && !identical(
                                    gsub(" ", "", as.character(tempChr[z, 3])),
                                    gsub(" ", "", as.character(tempChr[z + 1, 2]))
                                )
                            ) {
                                excoord = rbind(
                                    excoord,
                                    cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4])
                                )
                            }
                        } else if (
                            z == nrow(tempChr)
                            & !identical(
                                gsub(" ", "", as.character(tempChr[z, 3])),
                                gsub(" ", "", as.character(
                                    ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep = "")), ref_table),][2]
                                ))
                            )
                        ) {
                            excoord = rbind(
                                excoord,
                                cbind(
                                    tempChr[z, 1],
                                    tempChr[z, 3],
                                    ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep = "")), ref_table),][2], tempChr[z, 4]
                                )
                            )
              
                        } else if (
                            z < nrow(tempChr)
                            & z > 1
                            && !identical(
                                gsub(" ", "", as.character(tempChr[z, 3])),
                                gsub(" ", "", as.character(tempChr[z + 1, 2]))
                            )
                        ) {
                            # if this interval is not inside another
                            excoord = rbind(
                                excoord,
                                cbind(
                                    tempChr[z, 1],
                                    tempChr[z, 3],
                                    tempChr[z + 1, 2],
                                    tempChr[z, 4]
                                )
                            )
                        }
                    }
                }
            }
        }
    }
  
    if (grepl("ider", Cyto_sample[coln])) {
        if (any(!grepl(Mainchr, coord[, 1]))) {
            coord <- rbind(coord, coord[grep(Mainchr, coord[, 1], invert = T), ])
        }
    }
  
    if (length(coord) > 0) {
        coord[, 2:3] <- apply(
            coord[, 2:3],
            2,
            function(x) {
                as.numeric(as.character(x))
            }
        )
    
        # make it so coord only takes what is in the mainchr
        #  if(!any(grepl("ider\\(",coord[,4])))
        #   {
        #     coord<-coord[grep(paste(Mainchr,collapse="|"),coord[,1]),]
        #   }
        coord[, 4] <- paste(addBool, coord[, 4], sep = "")
    }
  
    # probably have to rework this because of inversions and insersions
    if (
        length(excoord) > 0
        && !is.vector(excoord) && ncol(excoord) > 1
    ) {
        excoord[, 2:3] <- apply(
            excoord[, 2:3],
            2,
            function(x) {
                as.numeric(as.character(x))
            }
        )
    
        # as numeric as character
        # delete any rows in which first number is larger than the second-
            # inversions occur so dont do this
        bads <- vector(length = 0)
    
        for (i in 1:nrow(excoord)) {
            if (as.numeric(excoord[i, 2]) > as.numeric(excoord[i, 3])) {
                bads <- c(bads, i)
            }
        }

        if (length(bads) > 0) {
            if (nrow(excoord[-bads, ]) > 0) {
                excoord <- excoord[-bads, ]
            } else {
                excoord <- data.frame()
            }
        }

        if (is.vector(excoord) && length(excoord) > 0) {
            excoord[4] <- paste(addBool, excoord[4], sep = "")
        } else if (ncol(excoord) == 1) {
            excoord <- t(excoord)
            excoord[4] <- paste(addBool, excoord[4], sep = "")
        } else if (length(excoord) > 0) {
            excoord[, 4] <- paste(addBool, excoord[, 4], sep = "")
        }
    } else if (ncol(excoord) == 1) {
        excoord <- t(excoord)
        ecoord[2:3] <- as.numeric(as.character(excoord[2:3]))
        excoord[4] <- paste(addBool, excoord[4], sep = "")
    } else if (is.vector(excoord) && length(excoord) > 0) {
        ecoord[2:3] <- as.numeric(as.character(excoord[2:3]))
        excoord[4] <- paste(addBool, excoord[4], sep = "")
    }
  
    # make sure this does what i want it to
    # excoord<-unique(excoord)
    # coord<-unique(coord)
  
    # have to fix derivative chromosomes with two translocations,
        # if they overlap in a chromosome, its the overlapping increment that maters
  
    # do multi (X2) processing right now
    # for over X2 times, its a gain
    if (
        (
            nrow(coord) > 0
            && any(grepl("multi", coord[, 4]))
        )
        && (
            constitutional == T
            | (
                constitutional == F
                & !grepl("(c$)|(c\\?$)", coord[, 4])
            )
        )
    ) {
        n <- multi - 1
        multimastercoord <- coord
        multimasterexcoord <- excoord
        for (f in 1:n) {
            if (f > 1) {
                multitemp <- coord[grep("multi", coord[, 4]), ]
                multitemp[, 4] <- paste("+", multitemp[, 4], sep = '')
                multimastercoord <- rbind(multimastercoord, multitemp)
                if (
                    any(
                        nrow(excoord) > 0
                        && grepl("multi", excoord[, 4])
                    )
                ) {
                    multitemp <- excoord[grep("multi", excoord[, 4]), ]
                    multitemp[, 4] <- paste("+", multitemp[, 4], sep = '')
                    multimasterexcoord <- rbind(multimasterexcoord, multitemp)
                }
            } else {
                multimastercoord <- rbind(coord, coord[grep("multi", coord[, 4]), ])
                if (
                    any(
                        nrow(excoord) > 0
                        && grepl("multi", excoord[, 4])
                    )
                ) {
                    multimasterexcoord <- rbind(
                        multimasterexcoord,
                        excoord[grep("multi", excoord[, 4]), ]
                    )
                }
            }
        }

        coord <- multimastercoord
        excoord <- multimasterexcoord
    }
  
    listCoord <- list(
        coord,
        excoord,
        xmod,
        ymod,
        Mainchr,
        multi,
        transloctable,
        addtot
    )

    return(listCoord)
}