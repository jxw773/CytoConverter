## function for getting cytobands for translocations and insertions
## (taking compliment of stuff not included in discription)

getCytoBands <-
  function(Cyto_ref_table,
           Cyto_sample,
           lengthcount,
           o,
           temp,
           coln,
           derMods) {
    ##must take into account acen and qter pter and only one listing (will have to relte to two), reuse later code for this
    chr_table <-
      Cyto_ref_table[grep(paste(paste("chr", temp[[(lengthcount * 2 - 1)]][o], sep =
                                        ""), "$", sep = ""), Cyto_ref_table),]
    
    ##for isoderivative chromosomes, end point is potentially different, make boolean now
    
    isiso = FALSE
    
    ##quit if karyotype returns false
    earlyReturn = F
    
    if (any(grepl("ider", Cyto_sample[coln]) &
            grepl("t\\(", Cyto_sample[coln])))
    {
      isiso = TRUE
      arm <- gsub("[[:digit:]]", "", temp[[(grep("ider", derMods)) + 1]])
      if (arm == "q")
      {
        unused = "p"
      }
      if (arm == "p")
      {
        unused = "q"
      }
    }
    
    
    currentvec <- vector()
    
    
    if (any(grepl("::", temp[[lengthcount * 2]][i]) |
            grepl("~>", temp[[lengthcount * 2]][i]) |
            grepl("->", temp[[lengthcount * 2]][i]))) {
      #parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##find p or q, take stuff before take stuff after, before is chromosomes after is positions, this will return 2 objects, must take into account
      ##parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##only splits first one
      longform_table <-
        strsplit(strsplit(temp[[lengthcount * 2]][i], "::")[[1]], "(~>)|(->)")
      ##take away any front loaded :
      longform_table <-
        lapply(longform_table, function(x) {
          gsub(':', '', x)
        })
      
      
      
      in_table = data.frame()
      
      ##mark for dic, trc
      if (grepl("dic|trc", derMods[lengthcount]))
      {
        addBool <- paste("long", addBool, sep = '')
      }
      
      ##get data for each read of something->something
      for (j in 1:length(longform_table))
      {
        stringdx <- str_locate_all(pattern = "p|q", longform_table[[j]])
        chr_name_long <-
          substr(longform_table[[j]][1], 0, stringdx[[1]][1] - 1)
        positions <-
          as.vector(cbind(
            substr(longform_table[[j]][1],
                   stringdx[[1]][1],
                   nchar(longform_table[[j]][1])),
            substr(longform_table[[j]][2],
                   stringdx[[2]][1],
                   nchar(longform_table[[j]][2]))
          ))
        
        if (nchar(chr_name_long) != 0)
        {
          chr_table_2 <-
            Cyto_ref_table[grep(paste(paste("chr", chr_name_long, sep = ""), "$", sep = ""), Cyto_ref_table),]
        } else {
          chr_table_2 <- chr_table
        }
        
        ##account for terminal ends
        if (any(grepl("pter", positions)))
        {
          positions[grep("pter", positions)] <- chr_table_2[1, 4]
        }
        if (any(grepl("qter", positions)))
        {
          positions[grep("qter", positions)] <-
            chr_table_2[nrow(chr_table_2), 4]
        }
        ##be careful on centromeneter
        if (any(grepl("cen", positions)))
        {
          ##likey will have to be careful about this one
          positions[grep("cen", positions)] <-
            chr_table_2[grep("acen", chr_table_2[, 5]),][1, 4]
        }
        
        ##account for p10/q10
        ##double check naming convention for this one
        ##have to change the one for q
        positions[grep("p10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]),][1, 4]
        positions[grep("q10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]),][2, 4]
        
        ##make sure positions is in order to be processed correctly
        positions <- positionsorter(positions)
        
        ##put stuff in table
        positions_table <-
          matrix(chr_table_2[grep(paste(positions, collapse = "|", sep = "|"),
                                  chr_table_2[, 4]),], ncol = 5)
        
        
        if (is.vector(positions_table))
        {
          positions_table <- t(positions_table)
        }
        currentvec <-
          c(currentvec,
            as.vector(positions_table[, 4])[1],
            as.vector(positions_table[, 4])[length(as.vector(positions_table[, 4]))])
        
      }
      
      
      
      
    } else {
      positions <-
        strsplit(gsub("q", ",q", gsub("p", ",p", temp[[lengthcount * 2]][o])), ",")[[1]][2:length(strsplit(gsub("q", ",q", gsub(
          "p", ",p", temp[[lengthcount * 2]][o]
        )), ",")[[1]])]
      ##have to change the one to q
      positions[grep("p10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]),][1, 4]
      positions[grep("q10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]),][2, 4]
      positions <- positionsorter(positions)
      
      
      
      
      
      ###########################################################################################################
      #####################for mitelman data only###############################################################
      #####check for more than one band per chromosome for translocations######################################
      ##############################################################################################################
      
      if (forMtn == T &
          grepl("t\\(", derMods[lengthcount]) & length(positions) > 1)
      {
        earlyReturn = T
        
      } else{
        ################################################
        ##############################################
        #########################################
        ##if only one q or p , add on end point
        ##something is wrong here
        if (length(positions) == 1)
        {
          if (any(grepl("p", positions)))
          {
            ##currentvec <-
            ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))]
            ##           , 4], chr_table[nrow(chr_table), 4], sep = ''))
            currentvec <-
              c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))] +
                                              1, 4], chr_table[nrow(chr_table), 4], sep = ''))
            
          }
          if (any(grepl("q", positions)))
          {
            currentvec <-
              ##c(currentvec, paste(chr_table[grep(paste(positions,sep="",collapse="|"), chr_table[, 4])[length(grep(paste(positions,sep="",collapse="|"), chr_table[, 4]))] -
              ##                              1, 4], chr_table[1, 4], sep = ''))
              c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[1]
                                            - 1, 4], chr_table[1, 4], sep = ''))
            ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])-1
            ##                              , 4][1], chr_table[1, 4], sep = ''))
            ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])
            ##                          , 4][1], chr_table[1, 4], sep = ''))
          }
        } else{
          ## restrict if positions are at the ends of the chromosome
          pos <-
            grep(paste(positions, sep = "", collapse = "|"), chr_table[, 4])
          if (pos[length(pos)] + 1 > nrow(chr_table) &&
              pos[1] == 1)
          {
            currentvec <- c(currentvec,
                            paste(chr_table[nrow(chr_table), 4], sep = ''),
                            paste(chr_table[1, 4], sep = ''))
            
          } else if (pos[length(pos)] + 1 > nrow(chr_table)) {
            ##if cyto is at the top
            currentvec <-   c(
              currentvec,
              paste(chr_table[nrow(chr_table), 4], sep = ''),
              paste(chr_table[pos[1] - 1, 4], chr_table[1, 4], sep = '')
            )
            
            
          } else if (pos[1] == 1) {
            ##if cyto is on the bottom
            
            currentvec <-    c(
              currentvec,
              paste(chr_table[pos[length(pos)] +
                                1, 4], chr_table[nrow(chr_table), 4], sep = ''),
              paste(chr_table[1, 4], sep = '')
            )
            
            
          } else{
            ##if not
            currentvec <-
              c(
                currentvec,
                paste(chr_table[pos[length(pos)] +
                                  1, 4], chr_table[nrow(chr_table), 4], sep = ''),
                paste(chr_table[pos[1] -
                                  1, 4], chr_table[1, 4], sep = '')
              )
            ##c(
            ##currentvec,
            ##paste(chr_table[grep(paste(positions,collapse=""), chr_table[, 4])[length(grep(paste(positions,collapse=""), chr_table[, 4]))]
            ##                , 4], chr_table[nrow(chr_table), 4], sep = ''),
            ##paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))]
            ##                , 4], chr_table[1, 4], sep = '')
            ##)
          }
        }
      }
      
      
    }
    
    ##handling isoderivatives
    if (isiso)
    {
      currentvec <-
        sapply(currentvec, function(x) {
          gsub(paste(unused, "[[:digit:]]+(\\.[[:digit:]])*", sep = ''),
               temp[[(grep("ider", derMods)) + 1]],
               x)
        })
    }
    
    return(list(currentvec, earlyReturn))
  }
