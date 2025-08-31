
mod_utils <- modules::use('modules/utils.R')
mod_cytobands <- modules::use('modules/cytobands.R')
mod_merge <- modules::use('modules/merge.R')

##function that parses gains and losses
gainloss<-function(temp_table,
                   original_temp_table, 
                   ex_fusion_table,
                   original_ex_fusion_table,
                   Mainchr,
                   multi,
                   j,
                   cyto_ref_table,
                   ref_table,
                   Cyto_sample,
                   Con_data,
                   transloctable,
                   Dump_table,
                   constitutional,
                   guess,
                   guess_q,
                   guess_by_first_val,
                   forMtn,
                   orOption,
                   sexstimate,
                   normX,  
                   normY,   
                   xcount, 
                   ycount ,  # counts number of y in 2nd slot
                   xadd ,    # counts if +X occurs
                   yadd ,    # counts if +Y occurs
                   xmod ,    # counts modifications that arent whole chromosome add/del for X
                   ymod ,    # counts modifications that arent whole chromosome add/del for Y
                   xdel ,    # counts if -X occures as constitutional
                   ydel ,    # counts if -Y occures as consitutional
                   
                   xconstitutional ,  # shift counts for xc indications (kind of a correction factor)
                   yconstitutional ,  # shift counts for yc indications (kind of a correction factor)
                   
                   idealx ,  # estimate of what the x value should be
                   idealy ,  # estimate of what the y value should be
                   
                   addtot ,  # counts total "new chromosomes"
                   deltot ,  # counts total complete chrom deletions
                   modtot ,  # for idems only, counts modification chromosomes
                   
                   n ,     # ploidy count
                   ploidy , # ploidy non additive ##default 2 for diploid
                   
                   startcol ){
# print(multi)
# for future implementation: switch to searching at ex table if temp table is empty,
# default is temp_table
pointerConditional <- temp_table
deletions<-F

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
  
  
  if (length(ex_fusion_table) > 0) {
    colnames(ex_fusion_table) <- c("Chr", "Start", "End", "Type")
    ex_fusion_table[, 4] <- as.character(ex_fusion_table[, 4])
    ex_fusion_table[, 2:3] <- apply(
      ex_fusion_table[, 2:3],
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
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4]
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
  
  # temp_table[grep("+\(.\)",ex_fusion_table[,4]),4]<-"Gain"
  # ex_fusion_table[,4]<-"Gain"
  # rbind(temp_table,ex_fusion_table)
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
  # if (any(grepl("(t|ins)\\(", ex_fusion_table[, 4])))
  # {
  #   transloctable <-
  #   rbind(transloctable, ex_fusion_table[grep("(t|ins)\\(", ex_fusion_table[, 4]), ])
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[grep("LongDer", ex_fusion_table[, 4]), 4] <- "Loss"
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[
        intersect(
          intersect(
            grep(
              paste("chr", Mainchr, sep = "", collapse = "|"),
              ex_fusion_table[, 1]
            ),
            grep("((t)|(ins))\\(", ex_fusion_table[, 4])
          ),
          grep("^((t)|(ins))\\(", ex_fusion_table[, 4], invert = T)
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
        if (nrow(ex_fusion_table) > 0) {
          ex_fusion_table[, 4]
          
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[
        intersect(
          grep(
            paste(paste("chr", Mainchr, sep = ""), collapse = "|"),
            ex_fusion_table[, 1]
          ),
          grep(
            "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
            ex_fusion_table[, 4]
          )
        ),
        4
      ] <- "Loss"
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[grep("i\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
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
      ] <- "#ider::isochrom::der"
      
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[grep("ider\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  if (any(grepl("idic\\(.*", temp_table[, 4]))) {
    mod <- TRUE
    # Need to either break it up or replace,
    # think about this one too
    # deleted <- ex_fusion_table[grep("idic\\(.*&del\\(", ex_fusion_table[, 4]), 4]
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
    
    if (nrow(ex_fusion_table) > 0) {
      original_ex_fusion_table <- ex_fusion_table
      ex_fusion_table[grep("idic\\(.*", ex_fusion_table[, 4]), 4] <- "Gain"
      additionTable[[2]] <- c(
        additionTable[[2]],
        additionTable[[2]][grep("idic\\(.*", original_ex_fusion_table[, 4])]
      )
      ex_fusion_table <- rbind(
        ex_fusion_table,
        original_ex_fusion_table[grep("idic\\(.*", original_ex_fusion_table[, 4]),]
      )
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
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
      
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[grep("dic\\(.*", ex_fusion_table[, 4]) , 4] <- "Loss"
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
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
      
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[
          intersect(
            intersect(
              grep("dic\\(.*", ex_fusion_table[, 4]),
              grep(
                "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                ex_fusion_table[, 4],
                invert = T
              )
            ),
            grep("del\\(.*|add\\(.*", ex_fusion_table[, 4])
          ),
          4
        ] <- gsub(
          "del\\(.*|add\\(.*",
          "",
          ex_fusion_table[
            intersect(
              intersect(
                grep("dic\\(.*", ex_fusion_table[, 4]),
                grep(
                  "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                  ex_fusion_table[, 4],
                  invert = T
                )
              ),
              grep("del\\(.*|add\\(.*", ex_fusion_table[, 4])
            ),
            4
          ]
        )
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      
    }
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  # make sure this is reversed for long form
  if (any(grepl("trc\\(.*", temp_table[, 4]))) {
    mod <- TRUE
    if (any(grepl("long", temp_table[, 4]))) {
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[
          intersect(
            grep("trc\\(.*", ex_fusion_table[, 4]),
            grep(
              paste("chr", Mainchr[1], "chr", Mainchr[3], sep = '|'),
              ex_fusion_table[, 1]
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
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[
        intersect(
          grep("trc\\(.*", ex_fusion_table[, 4]),
          grep("chr", Mainchr[2], ex_fusion_table[, 1])
        ),
        4
      ] <- "Loss"
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  if (any(grepl("rob\\(", temp_table[, 4]))) {
    mod <- TRUE
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[grep("rob\\(", ex_fusion_table[, 4]), 4] <- "Loss"
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  if (any(grepl("^r\\(.*|^\\+r\\(.*", temp_table[, 4]))) {
    mod <- TRUE
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[grep("r\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  if (any(grepl("trp\\(.*", temp_table[, 4]))) {
    mod <- TRUE
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    
    # JP: ex_fusion_table not in rbind??
    temp_table <- rbind(temp_table, temp_table[grep("trp\\(.*", temp_table[, 4]),])
    
  }
  
  if (any(grepl("qdp\\(.*", temp_table[, 4]))) {
    mod <- TRUE
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
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
    
    if (nrow(ex_fusion_table) > 0) {
      if (is.vector(ex_fusion_table)) {
        ex_fusion_table[, 4] <- ""
        
      } else {
        ex_fusion_table[, 4] <- ""
        
      }
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  if (any(grepl("add", temp_table[, 4]))) {
    mod <- TRUE
    if (nrow(temp_table) > 0) {
      temp_table[grep("add", temp_table[, 4]), 4] <- "Loss"
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      if (is.vector(ex_fusion_table)) {
        ex_fusion_table[, 4] <- ""
        
      } else {
        ex_fusion_table[, 4] <- ""
        
      }
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    temp_table <- rbind(temp_table, ex_fusion_table)
    
  }
  
  # keep X and Y chromosomes here for ex table for later processing, may have to 
  # modify later for autosomes, only if modifications aren't triggered
  if (mod == FALSE) {
    if (is.vector(ex_fusion_table)) {
      ex_fusion_table[4] <- ""
      
    } else {
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[, 4] <- ""
        
      }
      
    }
    
    if (nrow(temp_table) > 0) {
      temp_table[, 4] <- paste(additionTable[[1]], temp_table[, 4], sep = "")
      
    }
    
    if (nrow(ex_fusion_table) > 0) {
      ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
    }
    temp_table <- rbind(temp_table, ex_fusion_table)
    
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

return(list(temp_table,
       xadd ,    # counts if +X occurs
       yadd ,    # counts if +Y occurs
       xmod ,    # counts modifications that arent whole chromosome add/del for X
       ymod ,    # counts modifications that arent whole chromosome add/del for Y
       xdel ,    # counts if -X occures as constitutional
       ydel ,    # counts if -Y occures as consitutional
       xconstitutional ,  # shift counts for xc indications (kind of a correction factor)
       yconstitutional ,  # shift counts for yc indications (kind of a correction factor)
       idealx,  # estimate of what the x value should be
       idealy ,  # estimate of what the y value should be
       addtot ,  # counts total "new chromosomes"
       deltot ,  # counts total complete chrom deletions
       modtot ,  # for idems only, counts modification chromosomes
       n ,     # ploidy count
       ploidy , # ploidy non additive ##default 2 for diploid
       startcol ))
}





##this is going to be the function that parses fusions 

fusion<-function(temp_fusion_table,
                 original_temp_table, 
                 ex_fusion_table,
                 original_ex_table,
                 Mainchr,
                 multi,
                 j,
                 cyto_ref_table,
                 ref_table,
                 Cyto_sample,
                 Con_data,
                 transloctable,
                 Dump_table,
                 constitutional,
                 guess,
                 guess_q,
                 guess_by_first_val,
                 forMtn,
                 orOption,
                 sexstimate,
                 normX,  
                 normY,   
                 xcount, 
                 ycount ,  # counts number of y in 2nd slot
                 xadd ,    # counts if +X occurs
                 yadd ,    # counts if +Y occurs
                 xmod ,    # counts modifications that arent whole chromosome add/del for X
                 ymod ,    # counts modifications that arent whole chromosome add/del for Y
                 xdel ,    # counts if -X occures as constitutional
                 ydel ,    # counts if -Y occures as consitutional
                 
                 xconstitutional ,  # shift counts for xc indications (kind of a correction factor)
                 yconstitutional ,  # shift counts for yc indications (kind of a correction factor)
                 
                 idealx ,  # estimate of what the x value should be
                 idealy ,  # estimate of what the y value should be
                 
                 addtot ,  # counts total "new chromosomes"
                 deltot ,  # counts total complete chrom deletions
                 modtot ,  # for idems only, counts modification chromosomes
                 
                 n ,     # ploidy count
                 ploidy , # ploidy non additive ##default 2 for diploid
                 
                 startcol){
  # print(multi)
  # for future implementation: switch to searching at ex table if temp table is empty,
  # default is temp_fusion_table
  pointerConditional <- temp_fusion_table
  
  if (nrow(temp_fusion_table) > 0) {
    # if its a translocated vector, fix now
    if (nrow(temp_fusion_table) == 4 & ncol(temp_fusion_table) == 1) {
      temp_fusion_table <- t(temp_fusion_table)
      colnames(temp_fusion_table) <- c("Chr", "Start", "End", "Type")
      
    }
    
    # special cases (if, else if, else if ,else if, else )
    # ider, i , ins, dic etc
    # else, grep stuff on 4th col to detirmine addition or deletions, dont include
    #   translations add etc
    # deal with + values differently
    
    # if it goes through special cases, flip mod = true, excoord does not combine at end
    mod <- FALSE
    temp_fusion_table[, 4] <- as.character(temp_fusion_table[, 4])
    colnames(temp_fusion_table) <- c("Chr", "Start", "End", "Type")
    
    temp_fusion_table[, 4] <- as.character(temp_fusion_table[, 4])
    colnames(temp_fusion_table) <- c("Chr", "Start", "End", "Type")
    
    if (length(ex_fusion_table) > 0) {
      colnames(ex_fusion_table) <- c("Chr", "Start", "End", "Type")
      ex_fusion_table[, 4] <- as.character(ex_fusion_table[, 4])
      ex_fusion_table[, 2:3] <- apply(
        ex_fusion_table[, 2:3],
        2,
        function(x) {
          as.numeric(as.character(x))
        }
      )
      
    
    
    temp_fusion_table[, 2:3] <- apply(
      temp_fusion_table[, 2:3],
      2,
      function(x) {
        as.numeric(as.character(x))
      }
    )
    
    
    # see if its an addition
    additionTable <- mod_utils$detectAdd(
      temp_fusion_table[, 4],
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[, 4]
      } else {
        NULL
      }
    )
    
    # if its not a del # think about this one for all cases
    
    
    # temp_fusion_table[grep("+\(.\)",ex_fusion_table[,4]),4]<-"Gain"
    # ex_fusion_table[,4]<-"Gain"
    # rbind(temp_fusion_table,ex_fusion_table)
    # incorporate main chr in these
    # do del and stuff up here, but keep+__ sign somehow
    
    # count chromosome loss properly here for special cases
    # sex chromosome count check how it interacts with xmod ymod
    if (
      !all(grepl("\\+", temp_fusion_table[, 4]))
      & any(!grepl("(^t\\()|(^ins\\()|(^inv\\()", temp_fusion_table[, 4]))
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
        any(grepl("-[:alpha:]", temp_fusion_table[, 4]))
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
    # if (any(grepl("(t|ins)\\(", ex_fusion_table[, 4])))
    # {
    #   transloctable <-
    #   rbind(transloctable, ex_fusion_table[grep("(t|ins)\\(", ex_fusion_table[, 4]), ])
    # }
    
    if (any(grepl("LongDer", temp_fusion_table[, 4]))) {
      mod <- TRUE
      temp_fusion_table[
        intersect(
          grep("LongDer", temp_fusion_table[, 4]),
          grep(
            paste(Mainchr, collapse = "|", sep = ''),
            temp_fusion_table[, 1],
            invert = T
          )
        ),
        4
      ] <- "#translocation::guess"
      
      if (nrow(ex_fusion_table) > 0) {
        ##ex_fusion_table[grep("LongDer", ex_fusion_table[, 4]), 4] <- "Loss"
        
      }
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    ##balanced translocations
  
    if(any(grepl("^ins\\(",temp_fusion_table[, 4]))){
      mod<-TRUE
      ##t_table<-transloctable[[grep(Cyto_sample[j],names(transloctable))]]
      if(length(Mainchr)==1)
      {}else{
            temp_fusion_table[
              intersect(
                grep(
                  paste("chr", Mainchr[2], sep = "", collapse = "|"),
                  temp_fusion_table[, 1]
                ),
                grep("^ins\\(", temp_fusion_table[, 4])),
              
              4
            ] <- paste("#insertion_chrom::inserted_piece|chrom_",(1),sep="")
            
            
            temp_fusion_table[
              intersect(
                grep(
                  paste("chr", Mainchr[1], sep = "", collapse = "|"),
                  temp_fusion_table[, 1]
                ),
                grep("^ins\\(", temp_fusion_table[, 4])),
              
              4
            ] <- paste("#insertion_chrom::recipiant_chrom_insertion_band|chrom_",(1),sep="")
            
            if (nrow(ex_fusion_table) > 0) {
              ex_fusion_table[
                intersect(
                  grep(
                    paste("chr", Mainchr[1], sep = "", collapse = "|"),
                    ex_fusion_table[, 1]
                  ),
                  grep("^ins\\(", ex_fusion_table[, 4])),
                
                4
              ] <- paste("#insertion_chrom::recipiant_chrom|chrom_",(1),sep="")
              
              
              temp_fusion_table<-rbind(temp_fusion_table,ex_fusion_table)
              
            }   
            
    } 
      
    }
    
    if(any(grepl("inv\\(",temp_fusion_table[,4]))){
      temp_fusion_table[
        intersect(
          grep(
            paste("chr", Mainchr[1], sep = "", collapse = "|"),
            temp_fusion_table[, 1]
          ),
          grep("^inv\\(", temp_fusion_table[, 4])),
        
        4
      ] <- paste("#inversion|chrom_",(1),sep="")
    }
    
    if(any(grepl("fra\\(",temp_fusion_table[,4]))){
      temp_fusion_table[
        intersect(
          grep(
            paste("chr", Mainchr[1], sep = "", collapse = "|"),
            temp_fusion_table[, 1]
          ),
          grep("^fra\\(", temp_fusion_table[, 4])),
        
        4
      ] <- paste("#fragile|chrom_",(1),sep="")
    }
    
    if(any(grepl("fis\\(",temp_fusion_table[,4]))){
      temp_fusion_table[
        intersect(
          grep(
            paste("chr", Mainchr[1], sep = "", collapse = "|"),
            temp_fusion_table[, 1]
          ),
          grep("^fis\\(", temp_fusion_table[, 4])),
        
        4
      ] <- paste("#fission_centromere|chrom_",(1),sep="")
    }
    
    
    if(any(grepl("^t\\(",temp_fusion_table[, 4]))){
      mod<-TRUE
      ##t_table<-transloctable[[grep(Cyto_sample[j],names(transloctable))]]
      
      for(p in 1:length(Mainchr))
      {
        if(p==length(Mainchr))
        {
          temp_fusion_table[
            intersect(
              grep(
                paste("chr", Mainchr[p], sep = "", collapse = "|"),
                temp_fusion_table[, 1]
              ),
              grep("^t\\(", temp_fusion_table[, 4])),
            
            4
          ] <- paste("#translocation_balanced|chrom_",(p),sep="")
          
          if (nrow(ex_fusion_table) > 0) {
            ex_fusion_table[
              grep(
                paste("chr", Mainchr[1], sep = "", collapse = "|"),
                ex_fusion_table[, 1]
              ),
              
              4
            ] <- paste("#translocation_balanced|chrom_",(p),sep="")
            
          } 
          
          
        }else{
          
          temp_fusion_table[
            intersect(
              grep(
                paste("chr", Mainchr[p], sep = "", collapse = "|"),
                temp_fusion_table[, 1]
              ),
              grep("^t\\(", temp_fusion_table[, 4])),
            
            4
          ] <- paste("#translocation_balanced|chrom_",(p),sep="")
          
          if (nrow(ex_fusion_table) > 0) {
            ex_fusion_table[
              grep(
                paste("chr", Mainchr[p+1], sep = "", collapse = "|"),
                ex_fusion_table[, 1]
              ),
              
              4
            ] <- paste("#translocation_balanced|chrom_",(p),sep="")
            
          } 
        }
      }
      temp_fusion_table<-rbind(temp_fusion_table,ex_fusion_table)
    }
    
    # Think about how this is handling issue
    #   Really check on this
    if (
      any(
        grepl("t\\(", temp_fusion_table[, 4])
        & !grepl("^t\\(", temp_fusion_table[, 4])
      )
    ) {
      mod <- TRUE
      temp_fusion_table[
        intersect(
          intersect(
            grep(
              paste("chr", Mainchr, sep = "", collapse = "|"),
              temp_fusion_table[, 1]
            ),
            grep("t\\(", temp_fusion_table[, 4])
          ),
          grep("^t\\(", temp_fusion_table[, 4], invert = T)
        ),
        4
      ] <- "#derivative_chromosome::translocation"
      
      
      temp_fusion_table[
        intersect(
          intersect(
            grep(
              paste("chr", Mainchr, sep = "", collapse = "|"),
              temp_fusion_table[, 1],invert=T
            ),
            grep("t\\(", temp_fusion_table[, 4])
          ),
          grep("^t\\(", temp_fusion_table[, 4], invert = T)
        ),
        4
      ] <- "#derivative_chromosome::translocation"
      
      if (nrow(ex_fusion_table) > 0) {
        ##ex_fusion_table[
        ##  intersect(
        ##    intersect(
        ##      grep(
        ##        paste("chr", Mainchr, sep = "", collapse = "|"),
        ##        ex_fusion_table[, 1],invert=T
        ##      ),
        ##      grep("t\\(", ex_fusion_table[, 4])
        ##    ),
        ##    grep("^t\\(", ex_fusion_table[, 4], invert = T)
        ##  ),
        ##  4
        ##] <- "#derivative_chrom::translocation"
        
      }
      
    }
    
    if (
      any(
        grepl("ins\\(", temp_fusion_table[, 4])
        & !grepl("^ins\\(", temp_fusion_table[, 4])
      )
    ) {
      
      ##create variable on order of ins 
      
      ins_chrom<-paste((strsplit(temp_fusion_table[intersect(grep("ins\\(", temp_fusion_table[, 4])
          ,grep("^ins\\(", temp_fusion_table[, 4], invert = T)),4],"ins\\(|;")[[1]][-1]),"$",sep="")
      
      mod <- TRUE
      
      if(Mainchr==ins_chrom[1])
      {
        temp_fusion_table[
          intersect(
            intersect(
              grep(
                paste("chr", Mainchr, sep = "", collapse = "|"),
                temp_fusion_table[, 1]
              ),
              grep("ins\\(", temp_fusion_table[, 4])
            ),
            grep("^ins\\(", temp_fusion_table[, 4], invert = T)
          ),
          4
        ] <- "#derivative_chrom::insertion_recipiant_chrom"
        
        
        temp_fusion_table[
          intersect(
            intersect(
              grep(
                paste("chr", Mainchr, sep = "", collapse = "|"),
                temp_fusion_table[, 1],
                invert = T
              ),
              grep("ins\\(", temp_fusion_table[, 4])
            ),
            grep("^ins\\(", temp_fusion_table[, 4], invert = T)
          ),
          4
        ] <- "#derivative_chrom::insertion_inserted_piece"
        
      }else if(Mainchr==ins_chrom[2]){
        temp_fusion_table[
          intersect(
            intersect(
              grep(
                paste("chr", Mainchr, sep = "", collapse = "|"),
                temp_fusion_table[, 1]
              ),
              grep("ins\\(", temp_fusion_table[, 4])
            ),
            grep("^ins\\(", temp_fusion_table[, 4], invert = T)
          ),
          4
        ] <- "#derivative_chrom::insertion_leftover_chrom"
        
       

      }
      
      if (nrow(ex_fusion_table) > 0) {
        ##ex_fusion_table[
        ##  intersect(
        ##    intersect(
        ##      grep(
        ##        paste("chr", Mainchr, sep = "", collapse = "|",invert=T),
        ##        ex_fusion_table[, 1]
        ##      ),
        ##      grep("ins\\(", ex_fusion_table[, 4])
        ##    ),
        ##    grep("^ins", ex_fusion_table[, 4], invert = T)
        ##  ),
        ##  4
        ## ] <- "#derivative_chrom::insertion_recipiant_chrom"
        
      }
      
    }
    
    # Take deletion areas out of translocation adjacent things
    if (
      any(
        grepl("del", temp_fusion_table[, 4])
        & !grepl("^del", temp_fusion_table[, 4])
        & !grepl("multi[[:digit:]]*del", temp_fusion_table[, 4])
      )
      | any(
        grepl("add", temp_fusion_table[, 4])
        & !grepl("^add", temp_fusion_table[, 4])
        & !grepl("multi[[:digit:]]*add", temp_fusion_table[, 4])
      )
    ) {
      if (
        nrow(temp_fusion_table)
        > length(
          union(
            intersect(
              grep("del", temp_fusion_table[, 4]),
              intersect(
                grep("^del", temp_fusion_table[, 4], invert = T),
                grep("multi[[:digit:]]*del", temp_fusion_table[, 4], invert = T)
              )
            ), 
            intersect(
              grep("add", temp_fusion_table[, 4]),
              intersect(
                grep("^add", temp_fusion_table[, 4], invert = T),
                grep("multi[[:digit:]]*add", temp_fusion_table[, 4], invert = T)
              )
            )
          )
        )
      ) {
        temp_fusion_table <- mod_merge$mergeDeletions(temp_fusion_table, Mainchr)
        deletions <- T
        # if it is a vector, convert
        if (is.vector(temp_fusion_table)) {
          temp_fusion_table <- as.data.frame(as.list(temp_fusion_table))
          temp_fusion_table[, 2:3] <- as.numeric(temp_fusion_table[, 2:3])
          colnames(temp_fusion_table) <- c("Chr", "Start", "End", "Type")
          
        } else if (nrow(temp_fusion_table) == 4 & ncol(temp_fusion_table) == 1) {
          # if its a translocated vector, fix now
          temp_fusion_table <- t(temp_fusion_table)
          temp_fusion_table <- as.data.frame(as.list(temp_fusion_table))
          temp_fusion_table[, 2:3] <- as.numeric(temp_fusion_table[, 2:3])
          colnames(temp_fusion_table) <- c("Chr", "Start", "End", "Type")
          
        }
        
        additionTable <- mod_utils$detectAdd(
          temp_fusion_table[, 4],
          if (nrow(ex_fusion_table) > 0) {
            ex_fusion_table[, 4]
            
          } else {
            NULL
            
          }
        )
        
      }
      
    }
    
    

    if ((
      any(grepl("((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))", temp_fusion_table[, 4]))
      & length(temp_fusion_table[, 4])
      > length(
        c(
          which(grepl("del", temp_fusion_table[, 4])),
          which(grepl("add", temp_fusion_table[, 4]))
        )
      )
      & length(temp_fusion_table[, 4]) > 1
      & (
        length(temp_fusion_table[, 4])
        - length(
          c(
            which(grepl("del", temp_fusion_table[, 4])),
            which(grepl("add", temp_fusion_table[, 4]))
          )
        )
      ) > 1
    ) || (
      any(grepl("((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))", temp_fusion_table[, 4]))
      && deletions
    )) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        temp_fusion_table[
          intersect(
            grep(
              paste(paste("chr", Mainchr, sep = ""), collapse = "|"),
              temp_fusion_table[, 1],
              invert = T
            ),
            grep(
              "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
              temp_fusion_table[, 4]
            )
          ),
          4
        ] <- "#derivative_chrom"
        
      }
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[
      ##    intersect(
      ##      grep(
      ##        paste(paste("chr", Mainchr, sep = ""), collapse = "|"),
      ##        ex_fusion_table[, 1]
      ##      ),
      ##      grep(
      ##        "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
      ##        ex_fusion_table[, 4]
      ##      )
      ##    ),
      ##    4
      ##  ] <- "Loss"
      
      ##}
      
      ##if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
      
      ##}
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
      ##}
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    if (any(grepl("i\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        original_temp_fusion_table <- temp_fusion_table
        temp_fusion_table[grep("i\\(.*", temp_fusion_table[, 4]), 4] <- "#isometric_chrom"
        ##additionTable[[1]] <- c(
        ##  additionTable[[1]],
        ##  additionTable[[1]][grep("i\\(.*", original_temp_fusion_table[, 4])]
        ##)
        temp_fusion_table <- rbind(
          temp_fusion_table,
          original_temp_fusion_table[grep("i\\(.*", original_temp_fusion_table[, 4]),]
        )
        # add to additionTable in light of new addition
        
      }
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##ex_fusion_table[grep("i\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
      
      ##}
      
      ##if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
      
      ##}
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
      ##}
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    # get whats included, dup, then rest on chromosome is deletion
    
    if (any(grepl("ider\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      # Have to be able to handle these cases
      # Do inclusion and exclusion based on deleted case
      if (nrow(temp_fusion_table) > 0) {
        original_temp_fusion_table <- temp_fusion_table
        temp_fusion_table[
          intersect(
            grep("ider\\(.*", temp_fusion_table[, 4]),
            grep("del\\(.*|add\\(.*", temp_fusion_table[, 4], invert = T)
          ),
          4
        ] <- "#derivative_chrom::isometric_chrom"
        
        ##additionTable[[1]] <- c(
        ##  additionTable[[1]],
        ##  additionTable[[1]][
        ##    intersect(
        ##      grep("ider\\(.*", temp_fusion_table[, 4]),
        ##      grep("del\\(.*|add\\(.*", temp_fusion_table[, 4], invert = T)
        ##    )
        ##  ]
        ##)
        
        temp_fusion_table <- rbind(
          temp_fusion_table,
          original_temp_fusion_table[grep("ider\\(.*", original_temp_fusion_table[, 4]),]
        )
        
        # activate deletions ? Think about this
        # JAN Thinking time
      }
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[grep("ider\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
      
      ##}
      
      ##if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
      
      ##}
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
      ##}
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    if (any(grepl("idic\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      # Need to either break it up or replace,
      # think about this one too
      deleted <- ex_fusion_table[grep("idic\\(.*&del\\(", ex_fusion_table[, 4]), 4]
      if (nrow(temp_fusion_table) > 0) {
        original_temp_fusion_table <- temp_fusion_table
        temp_fusion_table[
          intersect(
            grep("idic\\(.*", temp_fusion_table[, 4]),
            grep("dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", temp_fusion_table[, 4], invert = T)
          ),
         4
        ] <- "Loss"
    
      ##  additionTable[[1]] <- c(
      ##    additionTable[[1]],
      ##    additionTable[[1]][
      ##      intersect(
      ##        grep("del\\(.*|add\\(.*", original_temp_fusion_table[, 4], invert = T),
      ##        grep("idic\\(.*", original_temp_fusion_table[, 4])
      ##      )
      ##    ]
      ##  )
      
       temp_fusion_table <- rbind(
         temp_fusion_table,
          original_temp_fusion_table[
            intersect(
              grep("del\\(.*|add\\(.*", original_temp_fusion_table[, 4], invert = T),
              grep("idic\\(.*", original_temp_fusion_table[, 4])
            )
            ,
         ]
        )
      
      }
      
      if (nrow(ex_fusion_table) > 0) {
        original_ex_fusion_table <- ex_fusion_table
        ex_fusion_table[grep("idic\\(.*", ex_fusion_table[, 4]), 4] <- "#concentric_chrom::dicentric_chrom::isometric_chrom"
        ##additionTable[[2]] <- c(
        ##  additionTable[[2]],
        ##  additionTable[[2]][grep("idic\\(.*", original_ex_fusion_table[, 4])]
        ##)
        ex_fusion_table <- rbind(
          ex_fusion_table,
          original_ex_fusion_table[grep("idic\\(.*", original_ex_fusion_table[, 4]),]
        )
        
      }
      
      ## if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
      
      ##  }
      
      ##  if (nrow(ex_fusion_table) > 0) {
      ##    ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
      
      ##  }
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    # figure out why you did this
    # make sure this is reversed for long form
    if (
      any(
        which(
          grepl("dic\\(.*", temp_fusion_table[, 4])
          & !grepl("idic\\(.*", temp_fusion_table[, 4])
        )
      )
    ) {
      
      mod <- TRUE
      if (any(grepl("long", temp_fusion_table[, 4]))) {
        
        if (nrow(temp_fusion_table) > 0) {
          temp_fusion_table[, 4]<-"#centric_chrom::dicentric_chrom"
          ##temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
          
        }
        
        ##if (nrow(ex_fusion_table) > 0) {
        ##  ex_fusion_table[grep("dic\\(.*", ex_fusion_table[, 4]) , 4] <- "Loss"
        ##  ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
        ##}
        
      } else {
        ##if (nrow(temp_fusion_table) > 0) {
        ##  temp_fusion_table[
        ##    intersect(
        ##      grep("dic\\(.*", temp_fusion_table[, 4]),
        ##      grep(
        ##        "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
        ##        temp_fusion_table[, 4],
        ##        invert = T
        ##      )
        ##    ),
        ##    4
        ##  ] <- "Loss"
        
        ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
        ##}
        
        ##include these extras if they are part of a derivative chromosome
        
        if (nrow(ex_fusion_table) > 0) {
          ex_fusion_table[
            intersect(
              intersect(
                grep("dic\\(.*", ex_fusion_table[, 4]),
                grep(
                  "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                  ex_fusion_table[, 4],
                  invert = T
                )
              ),
              grep("del\\(.*|add\\(.*", ex_fusion_table[, 4])
            ),
            4
          ] <- gsub(
            "del\\(.*|add\\(.*",
            "",
            ex_fusion_table[
              intersect(
                intersect(
                  grep("dic\\(.*", ex_fusion_table[, 4]),
                  grep(
                    "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                    ex_fusion_table[, 4],
                    invert = T
                  )
                ),
                grep("del\\(.*|add\\(.*", ex_fusion_table[, 4])
              ),
              4
            ]
          )
          
          
          
          
          ex_fusion_table[
            intersect(
              intersect(
                grep("dic\\(.*", ex_fusion_table[, 4]),
                grep(
                  "idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*",
                  ex_fusion_table[, 4],
                  invert = T
                )
              ),
              grep("del\\(.*|add\\(.*", ex_fusion_table[, 4],invert=T)
            ),
            4
          ] <- "#centric_chrom::dicentric_chrom"
          
          ##ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
          
        }
        
      }
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    # make sure this is reversed for long form
    if (any(grepl("trc\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (any(grepl("long", temp_fusion_table[, 4]))) {
        
         temp_fusion_table[
            intersect(
               grep("trc\\(.*", temp_fusion_table[, 4]),
               grep(
               paste("chr", Mainchr[1], "chr", Mainchr[3],"chr",Mainchr[2], sep = '|'),
                temp_fusion_table[, 1]
              )
            ),
             4
           ] <- "#centric_chrom_tri"
         
        ##if (nrow(ex_fusion_table) > 0) {
        ##  ex_fusion_table[
        ##    intersect(
       ##       grep("trc\\(.*", ex_fusion_table[, 4]),
       ##       grep(
        ##        paste("chr", Mainchr[1], "chr", Mainchr[3], sep = '|'),
        ##        ex_fusion_table[, 1]
        ##      )
        ##    ),
       ##     4
       ##   ] <- "Loss"
          
        } else {
        if (nrow(temp_fusion_table) > 0) {
          temp_fusion_table[
            intersect(
              grep("trc\\(.*", temp_fusion_table[, 4]),
              grep("chr", Mainchr[2], temp_fusion_table[, 1])
            ),
            4
            ] <- "#centric_chrom_tri"
          
          
        }
        
      }
      
      ##if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      ##}
      
      if (nrow(ex_fusion_table) > 0) {
        
        
        ex_fusion_table[
          intersect(
            grep("trc\\(.*", ex_fusion_table[, 4]),
            grep(
              paste("chr", Mainchr[1], "chr", Mainchr[3], sep = '|'),
              ex_fusion_table[, 1]
            )
          ),
          4
          ] <- "#centric_chrom_tri"
        
      
      
    
    
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    if (any(grepl("rob\\(", temp_fusion_table[, 4]))) {
      mod <- TRUE
      temp_fusion_table[grep("rob\\(", temp_fusion_table[, 4]), 4] <- "#derivative_chrom::translocation_robertsonian"
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[grep("rob\\(", ex_fusion_table[, 4]), 4] <- "Loss"
        
      ##}
      
      ##if (nrow(temp_fusion_table) > 0) {
      ##  temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      ##}
      
      ##if (nrow(ex_fusion_table) > 0) {
      ##  ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      ##}
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    if (any(grepl("^r\\(.*|^\\+r\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[grep("r\\(.*", ex_fusion_table[, 4]), 4] <- "Loss"
        
      }
      
      if (nrow(temp_fusion_table) > 0) {
        #temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        temp_fusion_table[grep("r\\(.*", temp_fusion_table[, 4]), 4] <- "#ring_chrom"
        
      }
      
      
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    
    if (any(grepl("trp\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        #temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      }
      
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      
      # JP: ex_fusion_table not in rbind??
      temp_fusion_table <- rbind(temp_fusion_table, temp_fusion_table[grep("trp\\(.*", temp_fusion_table[, 4]),])
      
    }
    
    if (any(grepl("qdp\\(.*", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      }
      
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      
      # JP: why is this duplicated in rbind?
      # quadruple duplication, so its a gain 4 times 
       temp_fusion_table <- rbind(
       temp_fusion_table,
       temp_fusion_table[grep("qdp\\(.*", temp_fusion_table[, 4]),],
       temp_fusion_table[grep("qdp\\(.*", temp_fusion_table[, 4]),]
      )
      
    }
    
    
    # Add code here for +gain +loss (+gain=gain, +loss ="")
    if (any(grepl("del", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        temp_fusion_table[grep("del", temp_fusion_table[, 4]), 4] <- "Loss"
        ##temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      }
      
      if (nrow(ex_fusion_table) > 0) {
        if (is.vector(ex_fusion_table)) {
          ex_fusion_table[, 4] <- ""
          
        } else {
          ex_fusion_table[, 4] <- ""
          
        }
        ##ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    if (any(grepl("add", temp_fusion_table[, 4]))) {
      mod <- TRUE
      if (nrow(temp_fusion_table) > 0) {
        temp_fusion_table[grep("add", temp_fusion_table[, 4]), 4] <- "Loss"
        ##temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      }
      
      if (nrow(ex_fusion_table) > 0) {
        if (is.vector(ex_fusion_table)) {
          ex_fusion_table[, 4] <- ""
          
        } else {
          ex_fusion_table[, 4] <- ""
          
        }
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    # keep X and Y chromosomes here for ex table for later processing, may have to 
    # modify later for autosomes, only if modifications aren't triggered
    if (mod == FALSE) {
      if (is.vector(ex_fusion_table)) {
        ex_fusion_table[4] <- ""
        
      } else {
        if (nrow(ex_fusion_table) > 0) {
          ex_fusion_table[, 4] <- ""
          
        }
        
      }
      
      if (nrow(temp_fusion_table) > 0) {
        temp_fusion_table[, 4] <- paste(additionTable[[1]], temp_fusion_table[, 4], sep = "")
        
      }
      
      if (nrow(ex_fusion_table) > 0) {
        ex_fusion_table[, 4] <- paste(additionTable[[2]], ex_fusion_table[, 4], sep = "")
        
      }
      temp_fusion_table <- rbind(temp_fusion_table, ex_fusion_table)
      
    }
    
    # handle stuff for minus as well
    # if +Gain
    if (nrow(temp_fusion_table) > 0) {
      temp_fusion_table[grep("^-Gain|^-$", temp_fusion_table[, 4]), 4] <- "Loss"
      temp_fusion_table[grep("\\+Loss", temp_fusion_table[, 4]), 4] <- "+Loss"
      temp_fusion_table[grep("\\+Gain", temp_fusion_table[, 4]), 4] <- "Gain"
      temp_fusion_table[
        intersect(
          grep("\\+", temp_fusion_table[, 4]),
          grep("\\+Loss", temp_fusion_table[, 4], invert = T)
        ),
        4
      ] <- "Gain"
      
      # check other permutations of this/ dont think insertions belong here
      ##temp_fusion_table[grep("dup|qdp|tan|trp|\\+$", temp_fusion_table[, 4]), 4] <- "Gain"
      
      Plus_Loss <- temp_fusion_table[grep("\\+Loss", temp_fusion_table[, 4]), ]
      
      # cut off intersections early
      # Think aboiut this
      if (
        length(Plus_Loss) > 0
        && nrow(Plus_Loss) > 0
        && any(grepl("Gain", temp_fusion_table[, 4]) & any(grepl("Loss", temp_fusion_table[, 4])))
      ) {
        temp_fusion_table <- mod_merge$mergeTable(temp_fusion_table)
        
      }
      
      temp_fusion_table[grep("\\+Loss", temp_fusion_table[, 4]), 4] <- ""
      
    }
  }
  }  
  return(temp_fusion_table)
}