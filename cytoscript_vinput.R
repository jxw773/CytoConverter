
if(require("stringr")){
  print("stringr is loaded correctly")
} else {
  print("trying to install stringr")
  install.packages("stringr")
  if(require("stringr")){
    print("stringr installed and loaded")
  } else {
    stop("could not install stringr")
  }
}
  
if(require("stringi")){
  print("stringi is loaded correctly")
} else {
    print("trying to install stringr")
    install.packages("stringi")
if(require("stringi")){
      print("stringi installed and loaded")
    } else {
      stop("could not install stringi")
    }  
}

##function wrapper, activate at the end
CytoConverter<-function(in_data)
{
  
  
  ##table with desired output
  Final_table <- matrix(ncol = 5, nrow = 0)
  
  colnames(Final_table) <- c("Sample ID", "Chr", "Start", "End", "Type")
  
  ##convert any single string into a table
  if (is.vector(in_data))
  {
    in_data <- t(matrix(c("sample", in_data)))
  }
  
  ##dump table of stuff containing unprocessed reads
  ##Write this later
  ##get all fish recorded
  Dump_table <- matrix(ncol = 3, nrow = 0)
  ##double check that this does not delete later data potentially
  fish_table <- in_data[grep("ish.*$", in_data[, 2]), ]
  if (is.vector(fish_table))
  {
    Dump_table <- rbind(Dump_table, c(fish_table, "fish reading"))
  } else{
    Dump_table <- rbind(Dump_table, cbind(fish_table, "fish reading"))
  }
  
  ##now take out ish readings
  in_data[, 2] <- gsub("ish.*$", "", as.character(in_data[, 2]))
  
  ##set everything to lowercase, then set x and y to uppercase
  in_data[, 2] <- tolower(in_data[, 2])
  in_data[, 2] <- chartr("x", "X", in_data[, 2])
  in_data[, 2] <- chartr("y", "Y", in_data[, 2])
  
  Con_data = matrix(nrow = 0, ncol = 2)
  ##spliting cell lines
  for (i in 1:nrow(in_data))
  {
    if (!is.na(in_data[i, 2]) & nchar(as.character(in_data[i, 2])) > 0)
    {
      c2 <- strsplit(in_data[i, 2], split = "/")[[1]]
      c1 <- paste(in_data[i, 1], 1:length(c2), sep = "_")
      Con_data <- rbind(Con_data, cbind(c1, c2))
    } else{
      Con_data <-
        rbind(Con_data, cbind(paste(in_data[i, 1], "_1"), in_data[i, 2]))
    }
    
    
    ##make idems here
    if (any(grepl("idem|sl|sdl", Con_data[, 2])))
    {
      idem_index <- grep("idem|sl|sdl", Con_data[, 2])
      temp_data <-
        Con_data[(idem_index[1] - 1):idem_index[length(idem_index)], ]
      
      for (j in 2:nrow(temp_data))
      {
        if (grepl("idem|sl", temp_data[j, 2]))
        {
          prev <- unlist(strsplit(temp_data[1, 2], "\\["))[1]
        } else{
          ##implement this so it can handel two sl1 in sucession and sdl1 sdl2
          prev <- unlist(strsplit(temp_data[j - 1, 2], "\\["))[1]
          
        }
        prev <- unlist(strsplit(prev, ","))
        prev <- prev[2:length(prev)]
        sexchromprev <- prev[grep("^[XY]+", prev)[1]]
        autochromprev <- prev[grep("^[XY]+", prev, invert = T)]
        
        
        clonecount = 1
        if (grepl("idemx|idemX|slx|slX|sdl*x|sdl*X", temp_data[j, 2]))
        {
          ##check if there is a X3 etc value , if there is pick that up, store, add to addtot, make it process through twice later
          
          clonecount = unlist(strsplit(temp_data[j, 2], "idemx|idemX|slx|slX"))[2]
          if (grepl("-|~", clonecount))
          {
            clonecount <- unlist(strsplit(clonecount, "~|-"))[1]
            
          }
          clonecount = as.numeric(unlist(strsplit(
            as.character(clonecount), "\\[|,"
          ))[1])
        }
        cur <- unlist(strsplit(temp_data[j, 2], ","))
        if (!is.na(sexchromprev))
        {
          if (grepl("^[XY]+", cur[2]))
          {
            cur[2] <-
              paste(paste(rep(sexchromprev, clonecount), collapse = ''), cur[2], sep =
                      '')
          } else if (grepl("idem|sl|sdl", cur[2]))
          {
            cur[2] <- paste(rep(sexchromprev, clonecount), collapse = '')
          } else{
            cur <-
              c(cur[1], paste(rep(sexchromprev, clonecount), collapse = ''), cur[2:length(cur)])
            
          }
        }
        if (!is.null(autochromprev))
        {
          if (length(cur) > 2)
          {
            cur <- c(cur[1:2], rep(autochromprev, clonecount), cur[3:length(cur)])
          } else
          {
            cur <- c(cur[1:2], rep(autochromprev, clonecount))
          }
        }
        cur <- cur[grep("idem|sl|sdl", cur, invert = T)]
        temp_data[j, 2] <- paste(c(cur, "ids"), sep = '', collapse = ',')
      }
      Con_data[(idem_index[1] - 1):idem_index[length(idem_index)], ] <-
        temp_data
    }
    
    
    
  }
  
  rownames(Con_data)<-1:nrow(Con_data)
  
  ##specific loc
  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("cytoBand.txt", header = FALSE)
    ), as.character)
  
  ##ref_table <-sapply(as.data.frame(read.delim("GRCh38.d1.vd1.fa.fai",nrows = 25, header = FALSE)), as.character)
  ref_table <-as.data.frame(Cyto_ref_table[sapply(unique(Cyto_ref_table[,1]),function(x){grep(x,Cyto_ref_table[,1])[length(grep(paste(x,"$",sep=""),Cyto_ref_table[,1]))]}),][,1:2])
  ref_table<-apply(ref_table,2,as.character)
  
  
  Con_data[, 2] <- gsub(" ", "", Con_data[, 2])
  ##any rejoins will be ;; now
  ##actually lets depricate this
  ##Con_data[, 2] <- gsub(":", ";", Con_data[, 2])
  Con_data[, 2] <- gsub("crYp", "", Con_data[, 2])
  
  ## taking all reads with ? and ~ as well inot dump table
  unsure_table <- Con_data[grep("\\?|\\~|inc", Con_data[, 2]), ]
  Dump_table <- if(is.vector(unsure_table)){
    rbind(Dump_table, c(unsure_table, "?/~/ or incomplete karyotype detected"))
  }else{
    rbind(Dump_table, cbind(unsure_table, "?/~/ or incomplete karyotype detected"))
    
  }
  
  ##test con data group
  ##Con_data<-rbind(Con_data[521,],Con_data[525,],Con_data[632,],Con_data[701,],Con_data[781,],Con_data[788,],Con_data[76,],Con_data[138,],Con_data[483,],Con_data[5,],Con_data[6,],Con_data[10,])
  ##Con_data[621,],
  ##position match function, given chromosome and position, matches to UCSC table
  ##for a row
  
  ##take compliment of genomic coordinates
  
  ##convert to genomic coordinates (used in parser several times)
  
  ##function for getting cytobands for translocations and insertions (taking compliment of stuff not included in discription)
  getCytoBands <- function(lengthcount, o, temp,coln,derMods) {
    ##must take into account acen and qter pter and only one listing (will have to relte to two), reuse later code for this
    chr_table <-
      Cyto_ref_table[grep(paste(paste("chr", temp[[(lengthcount * 2-1)]][o], sep =
                                        ""), "$", sep = ""), Cyto_ref_table), ]
    
    ##for isoderivative chromosomes, end point is potentially different, make boolean now
    
    isiso=FALSE
    
    if(any(grepl("ider",Cyto_sample[coln])&grepl("t\\(",Cyto_sample[coln])))
    {
      isiso=TRUE
      arm<-gsub("[[:digit:]]","",temp[[(grep("ider",derMods))+1]])
      if(arm=="q")
      {
        unused="p"
      }
      if(arm=="p")
      {
        unused="q"
      }
    }
    
    
    currentvec <- vector()
    
    
    if (any(grepl("::", temp[[lengthcount * 2]][i]) |
            grepl("~>", temp[[lengthcount * 2]][i]))) {
      #parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##find p or q, take stuff before take stuff after, before is chromosomes after is positions, this will return 2 objects, must take into account
      ##parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##only splits first one
      longform_table <-
        strsplit(strsplit(temp[[lengthcount * 2]][i], "::")[[1]], "~>")
      ##take away any front loaded : 
      longform_table <-lapply(longform_table,function(x){gsub(':','',x)})
      
      in_table = data.frame()
      
      ##mark for dic, trc
      if(grepl("dic|trc",derMods[lengthcount]))
      {
        addBool<-paste("long",addBool,sep='')
      }
      
      ##get data for each read of something->something
      for (j in 1:length(longform_table))
      {
        stringdx <- str_locate_all(pattern = "p|q", longform_table[[j]])
        chr_name_long <-
          substr(longform_table[[j]][1], 0, stringdx[[1]][1] - 1)
        positions <-
          as.vector(cbind(
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
          ))
        
        if (nchar(chr_name_long) != 0)
        {
          chr_table_2 <-
            Cyto_ref_table[grep(paste(paste(
              "chr", chr_name_long, sep = ""
            ), "$", sep = ""), Cyto_ref_table), ]
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
          positions[grep("qter", positions)] <- chr_table_2[nrow(chr_table_2), 4]
        }
        ##be careful on centromeneter
        if (any(grepl("cen", positions)))
        {
          ##likey will have to be careful about this one
          positions[grep("cen", positions)] <-
            chr_table_2[grep("acen", chr_table_2[, 5]), ][1, 4]
        }
        
        ##account for p10/q10
        ##double check naming convention for this one
        ##have to change the one for q
        positions[grep("p10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]), ][1, 4]
        positions[grep("q10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]), ][2, 4]
        
        ##make sure positions is in order to be processed correctly
        positions<-positionsorter(positions)
        
        ##put stuff in table
        positions_table <-
          matrix(chr_table_2[grep(paste(positions, collapse = "|", sep = "|"),
                                  chr_table_2[, 4]), ], ncol = 5)
        
        
        
        if (is.vector(positions_table))
        {
          positions_table <- t(positions_table)
        }
        currentvec<-c(currentvec, as.vector(positions_table[,4])[1],as.vector(positions_table[,4])[length(as.vector(positions_table[,4]))])
        
      }
      
      
      
      
    } else {
      positions <-
        strsplit(gsub("q", ",q", gsub("p", ",p", temp[[lengthcount * 2]][o])), ",")[[1]][2:length(strsplit(gsub("q", ",q", gsub(
          "p", ",p", temp[[lengthcount * 2]][o]
        )), ",")[[1]])]
      ##have to change the one to q
      positions[grep("p10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]), ][1, 4]
      positions[grep("q10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]), ][2, 4]
      positions<-positionsorter(positions)
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
            c(currentvec, paste(chr_table[grep(paste(positions,sep="",collapse="|"), chr_table[, 4])[length(grep(paste(positions,sep="",collapse="|"), chr_table[, 4]))] +
                                            1, 4], chr_table[nrow(chr_table), 4], sep = ''))
          
        }
        if (any(grepl("q", positions)))
        {
          currentvec <-
            ##c(currentvec, paste(chr_table[grep(paste(positions,sep="",collapse="|"), chr_table[, 4])[length(grep(paste(positions,sep="",collapse="|"), chr_table[, 4]))] -
            ##                              1, 4], chr_table[1, 4], sep = ''))
            c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[1] 
                                          -1, 4], chr_table[1, 4], sep = ''))
          ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])-1 
          ##                              , 4][1], chr_table[1, 4], sep = ''))
          ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4]) 
          ##                          , 4][1], chr_table[1, 4], sep = ''))
        }
      } else{
        currentvec <-
          c(
            currentvec,
            paste(chr_table[grep(paste(positions,sep="",collapse="|"), chr_table[, 4])[length(grep(paste(positions,collapse="|"), chr_table[, 4]))] +
                              1, 4], chr_table[nrow(chr_table), 4], sep = ''),
            paste(chr_table[grep(paste(positions,sep="",collapse="|"), chr_table[, 4])[1] -
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
    
    ##handling isoderivatives 
    if(isiso)
    {
      currentvec<-sapply(currentvec,function(x){gsub(paste(unused,"[[:digit:]]+(\\.[[:digit:]])*",sep=''),temp[[(grep("ider",derMods))+1]],x)})
    }
    
    return(currentvec)
  }
  
  
  
  ##make sure the order of the bands is in the order we like (p to q highers ps to low, lower qs to high )
  positionsorter<-function(positions){
    if(length(positions)==2)
    {
      if(grepl("q",positions[1])&grepl("p",positions[2]))
      {
        tempstorage<-positions[2]
        positions[2]<-positions[1]
        positions[1]<-tempstorage
      }
      if(grepl("q",positions[1])&grepl("q",positions[2]))
      {
        bands<-sapply(positions,function(x){as.numeric(gsub("q","",x))})
        if(bands[1]>bands[2])
        {
          tempstorage<-positions[2]
          positions[2]<-positions[1]
          positions[1]<-tempstorage
        }
      }
      if(grepl("p",positions[1])&grepl("p",positions[2]))
      {
        bands<-sapply(positions,function(x){as.numeric(gsub("p","",x))})
        if(bands[1]<bands[2])
        {
          tempstorage<-positions[2]
          positions[2]<-positions[1]
          positions[1]<-tempstorage
        }
      }
    }
    return(positions)
  }
  
  
  
  ##function for separating normal data
  ##take into acc same chrom insestion
  ##add cen into here (pter qter analouge)
  
  parser <- function(coln, xmod, ymod, transloctable,addtot)
  {
    ##derivative chromosomes with translocations are a loss (on native chromosome) gain (on new chromosome)
    ##figure out what add bool is and consaolidate that
    ##string splits by ;, takes into account multiple chromosomes  odd values are chromosomes even are the positions, for derivaties, first one needs to be treated differently
    test <-
      strsplit(gsub("[\\(\\)]", "", regmatches(
        Cyto_sample[coln], gregexpr("\\(.*?\\)",  Cyto_sample[coln])
      )[[1]]), ";")
    temp <-
      strsplit(gsub("[\\(\\)]", "", regmatches(
        Cyto_sample[coln], gregexpr("\\(.*?\\)",  Cyto_sample[coln])
      )[[1]]), ";")
    ##loop to go through everything and extract table if its in short form
    coord <- data.frame()
    excoord <- data.frame()
    derMods <- Cyto_sample[coln]
    addBool <- ""
    ##right now this is doing both der stuff and adition stuff, want it just to do + stuff make dermods do other stuff
    Mainchr <- vector()
    Allchr <- vector()
    ##MiscChr<-matrix()
    ##remeber what this is supposed to do
    derMods <- strsplit(derMods, "\\)")[[1]]
    
    if (length(temp) > 0)
    {
      ##chooses main chr for exclusion
      Mainchr <- temp[[1]]
      Mainchr <- paste(Mainchr, "$", sep = "")
      Mainchr<-gsub("p|q","",Mainchr)
      
      
    }
    
    ##this things not working properly
    
    if (length(derMods) > 0)
    {
      ##do this later on, + = +1 addtot, others = -
      ## addtot<-addtot+multi
      ##if(grepl("\\+.*\\(.*",derMods[1]))
      ##{
      ##  addtot<-addtot+1
      ##}
      ##consider gsub
    }
    ##addBool<- paste("+",addBool,sep="")
    ##}
    ##make sure this matches up with the next thing
    ##if (grepl("der\\(|rec\\(", Cyto_sample[coln])&&!grepl("ider\\(.*", Cyto_sample[coln]))
    ##{
    
    ##remeber what this is supposed to do
    ##derMods<-paste(")(",derMods,sep="")
    ##}
    
    ## if it straight up describes derivative makeup afterwads
    if (grepl("der\\([0-9;]+\\)\\(|rec\\([0-9;]+\\)\\(", Cyto_sample[coln]) &
        !grepl("ider", Cyto_sample[coln]))
    {
      addBool <- paste(addBool, "LongDer", sep = "")
    } else if (grepl("der|rec", Cyto_sample[coln]) &
               !grepl("ider", Cyto_sample[coln]))
    {
      addBool <- paste(addBool, derMods[1], sep = "")
      ##something is going wrong where when i changed this
      ##only want more than one
      if (grepl(
        "der\\([[:alnum:]]+(;[[:alnum:]])*\\)[[:alpha:]]+|rec\\([[:alnum:]]+(;[[:alnum:]])*\\)[[:alpha:]]+",
        Cyto_sample[coln]
      ))
      {
        derMods <-
          strsplit(Cyto_sample[coln], "(der|rec)\\([[:digit:]XY]+(;[[:digit:]XY])*\\)")[[1]][2]
        derMods <- strsplit(derMods, "\\)")[[1]]
        temp <- tail(temp, length(temp) - 1)
      }
    }
    ##if temp length is greater than 2 (derivative chromosome, then, dermods 1 is master indicator, not including ders above)
    if(length(temp)>2&(!(grepl("der|rec", Cyto_sample[coln])))||grepl("ider",Cyto_sample[coln]))
    {
      addBool<-paste(addBool,derMods[1],sep='')
    }
    
    ##if(grepl("\\+.*\\(.*",derMods[1]))
    ##{
    ##adds +etc on if it starts with something with +something
    ##this part is being messed up
    multi = 1 ##check if there is a X3 etc value , if there is pick that up, store, add to addtot, make it process through twice later
    if (grepl("\\)X|\\)×", Cyto_sample[coln]))
    {
      multi = unlist(strsplit(Cyto_sample[coln], ")X|)×"))[2]
      if (grepl("-|~", multi))
      {
        multi <- unlist(strsplit(multi, "~|-"))[1]
        
      }
      multi = as.numeric(multi)
      addBool <- paste(addBool, "multi", multi, sep = "")
    }
    
    ##########################################################################################################
    ######################################X chromosomes Y chromosomes #######################################
    #############################################################################################################
    ##increment X in presence of x modifications here
    if (any(grepl("X", Mainchr)))
    {
      xmod <- xmod + length(grep("X", Mainchr))
    }
    if (any(grepl("Y", Mainchr)))
    {
      ymod <- ymod + length(grep("Y", Mainchr))
    }
    
    ##make sure you count + properly for ? marks
    if(any(grepl("\\?|\\~",Cyto_sample[coln])&grepl("\\+",Cyto_sample[coln])))
    {
      addtot<-addtot+(1*multi)
    }
    
    if ((length(test) > 1 | any(grepl("p|q",temp)))  & !any(grepl("\\?|\\~", temp)))
    {
      ##goes by steps of 2, odd indexes indicate chromosomes, even indicate positions
      ##length_temp<-(if((length(temp) / 2)==0.5){1}else{length(temp)/2})
      for (lengthcount in 1:(if(!is.integer((length(temp) / 2))){ceiling(length(temp)/2)}else{length(temp)/2}))
      {
        Allchr <-
          c(Allchr, as.vector(paste(temp[[(lengthcount * 2-1)]], "$", sep = "")))
        
        ##handle those t(__;__) and ins (__;__) here with no follow up
        if (((lengthcount * 2-1)) == length(temp) || if(length(temp) > lengthcount*2-1){
          all(grepl("^[[:digit:]]+$", temp[[lengthcount * 2]]))}else{FALSE})
        {
          #######
          ########
          ########
          ##think about inproving it instead of else then, if if 
          ######
          #####
          #####
          #########################
          if (any(grepl("[pq]", temp[[(lengthcount * 2-1)]])))
          {
            if(grepl("p",temp[[(lengthcount*2-1)]]))
            {
              arm="p10"
            }
            if(grepl("q",temp[[(lengthcount*2-1)]]))
            {
              arm="q10"
            }
            temp<-c(if((lengthcount*2-1)!=1){temp[1:(lengthcount*2-1)]}else{temp[(lengthcount*2-1)]},arm,if(length(temp)>=2){temp[[lengthcount*2:length(temp)]]})
            temp[[(lengthcount*2-1)]]<-gsub("p|q","",temp[[(lengthcount*2-1)]])
            
          }
          if (grepl("(t|ins)\\(", derMods[(lengthcount * 2-1)]))
          {
            transchrom <-
              str_extract(Cyto_sample[coln], paste(gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]),"\\)",sep=''))
            selectedTransTable<-as.matrix(transloctable[grep(gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom)),names(transloctable))][[1]])
            transDer<-selectedTransTable[grep(paste("der\\(",paste(sapply(Mainchr,function(x){substr(x,0,nchar(x)-1)}),sep='',collapse=';'),"\\)",sep=''),selectedTransTable[,1]),][2]
            transtemp <-
              strsplit(gsub("[\\(\\)]", "", regmatches(
                transDer, gregexpr("\\(.*?\\)",  transDer)
              )[[1]]), ";")
            #####
            ##this isnt working 
            ###
            temp<-c(if((lengthcount*2-1)!=1){temp[1:((lengthcount*2-1))]}else{transtemp[1]},transtemp[2],if(length(temp)>=lengthcount*2){temp[lengthcount*2:length(temp)]})
            
            ##instead of below, for now, replace with equivilent derivative chromosome 
            
            ##get stuff in transloctable that matches same chromosom
            ##transtable <-
            ##  cbind(transloctable[grep(paste("chr", Mainchr, sep = '', collapse = '|'),
            ##                           transloctable[, 1]), 1:3], derMods[(lengthcount * 2-1)])
            ##coord <- rbind(coord, transtable)
            ##temp <- temp[[-((lengthcount * 2-1))]]
            ##derMods <- derMods[-(lengthcount)]
            ##lengthcount <- lengthcount - 1
          }
        }
        
        if (any(grepl("[pq]", temp[[(lengthcount * 2-1)]]))) {
          ##get chromosome , at 10 to other half, get coordinates for that
          temp[[(lengthcount*2-1)]]<-gsub("p|q","",temp[[(lengthcount*2-1)]])
        }
        if (grepl("t\\(", derMods[(lengthcount * 2-1)])) {
          ##if translocation , do this , shouldnt be placed , should be placed one loop above
          ## if not in the table, add to table
          ##names are just t(1;12)
          transchrom <-
            str_extract(Cyto_sample[coln], paste(gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]),"\\)\\(.+?\\)",sep=''))
          ##if this is not labled
          if(is.na(transchrom))
          {
            transchrom <-
              str_extract(Cyto_sample[coln], paste(gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]),"\\)",sep=''))
          }
          regtranschrom<-gsub("\\+","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom)))
          if (!any(grepl(regtranschrom, names(transloctable))))
          {
            ##entire translocation only if nessesary
            ##str_extract(mem, "t\\([[:digit:]]+(;[[:digit:]])*?\\)\\(.+?\\)")
            temptrans <- data.frame()
            
            
            
            ##make stuff now
            for (o in 1:length(temp[[(lengthcount * 2-1)]]))
            {
              currentvec <- getCytoBands(lengthcount, o, temp,coln,derMods)
              ##vector of o cytobands ##includes everything but cytoband on der o indicated, double check how your doing this (going one band before, is this right?)
              
              
              ##above is all to find vectors of o that are excluded
              if (o == length(temp[[(lengthcount * 2-1)]]))
              {
                temptrans <-
                  rbind(temptrans, cbind(
                    paste("der(", temp[[(lengthcount * 2-1)]][o], ")", sep = ''),
                    paste(
                      "der(",
                      temp[[(lengthcount * 2-1)]][o],
                      ";",
                      temp[[(lengthcount * 2-1)]][1],
                      if (length(currentvec) > 1) {
                        ";"
                      },
                      rep(temp[[(lengthcount * 2-1)]][o], length(currentvec) - 1),
                      ")(",
                      currentvec[1],
                      ";",
                      temp[[lengthcount * 2]][1],
                      if (length(currentvec) > 1) {
                        ";"
                      },
                      currentvec[-1],
                      ")",
                      sep = '',
                      collapse = ";"
                    ),
                    ""
                  ))
                
              } else{
                temptrans <-
                  rbind(temptrans, cbind(
                    paste("der(", temp[[(lengthcount * 2-1)]][o], ")", sep = ''),
                    paste(
                      "der(",
                      temp[[(lengthcount * 2-1)]][o],
                      ";",
                      temp[[(lengthcount * 2-1)]][o + 1],
                      if (length(currentvec) > 1) {
                        ";"
                      },
                      rep(temp[[(lengthcount * 2-1)]][o], length(currentvec) - 1),
                      ")(",
                      currentvec[1],
                      ";",
                      temp[[lengthcount * 2]][o + 1],
                      if (length(currentvec) > 1) {
                        ";"
                      },
                      currentvec[-1],
                      ")",
                      sep = '',
                      collapse = ";"
                    ),
                    ""
                  ))
                
              }
            }
            
            transloctable <- c(transloctable, list(temptrans))
            names(transloctable)[length(transloctable)] <-
              transchrom
          }
          
          if(!grepl("^\\++t\\(",derMods[(lengthcount*2-1)])&grepl("der",Cyto_sample[coln]))
          {
            ##this is a factor, not a table 
            transDer<-as.matrix(transloctable[[grep(regtranschrom,names(transloctable))]][grep(paste("der\\(",sapply(Mainchr,function(x){substr(x,0,nchar(x)-1)}),"\\)",sep="",collapse="|"),transloctable[[grep(regtranschrom,names(transloctable))]][,1]),])[2]
            ##if it is NA, get chromosome from previous translocation and use that instead
            if(is.na(transDer) & lengthcount > 1)
            {
              transDer<-as.matrix(transloctable[[grep(regtranschrom,names(transloctable))]][grep(paste("der\\(","(",paste(temp[[(lengthcount-1)*2-1]],collapse="|"),")","\\)",sep="",collapse="|"),transloctable[[grep(regtranschrom,names(transloctable))]][,1]),])[2]
            }
            ##translate derivative into temp stuff
            transtemp <-
              strsplit(gsub("[\\(\\)]", "", regmatches(
                transDer, gregexpr("\\(.*?\\)",  transDer)
              )[[1]]), ";")
            temp[[(lengthcount*2-1)]]<-transtemp[[1]]
            temp[[lengthcount*2]]<-transtemp[[2]]
            ##derMods[(lengthcount*2-1)]<-transDer
          }
          
          
          
          
        }
        if (grepl("ins\\(", derMods[(lengthcount * 2-1)])) {
          currentvec<-vector()
          inschrom <-str_extract(Cyto_sample[coln], paste(gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]),"\\)\\(.+?\\)",sep=''))
          
          if(is.na(inschrom))
          {
            inschrom <-
              str_extract(Cyto_sample[coln], paste(gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]),"\\)",sep=''))
          }
          reginschrom<-gsub("\\+","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",inschrom)))
          ## if insertion is on a single chromosome , do this , if for some reason the call is separated but the chromosome is not
          if(length(temp[[(lengthcount * 2-1)]])<length(temp[[lengthcount * 2]]))
          {
            temp[[(lengthcount * 2-1)]]<- rep(temp[[(lengthcount * 2-1)]],length(temp[[lengthcount * 2]]))
          }
          
          ##think about this harder
          ##if insertion is one whole sequence withint a chromosome
          ##if within is true, dont make another derivative chromosome
          within=F
          if(length(temp[[lengthcount*2-1]])==1 && length(temp[[(lengthcount * 2-1)]])<str_count(temp[[lengthcount*2]],"p|q"))
          {
            temp[[(lengthcount * 2-1)]]<- rep(temp[[(lengthcount * 2-1)]],2)
            
            ##location to split string (2nd q or p)
            splitloc<-str_locate_all(temp[[lengthcount*2]],"p|q")[[1]][,1][2]
            
            temp[[lengthcount*2]]<-c(str_sub(temp[[lengthcount*2]],1,splitloc-1),str_sub(temp[[lengthcount*2]],splitloc,str_length(temp[[lengthcount*2]])))
            within=T
          }
          
          if (!any(grepl(inschrom, names(transloctable))))
          {
            ##entire translocation only if nessesary
            tempins <- data.frame()
            
            
            currentvec <- getCytoBands(lengthcount, 1,temp,coln,derMods)
            
            ##fix this later , insertion is simply stopping where it left off instead of listing other half of the chromosome
            
            tempins <-
              rbind(tempins, cbind(
                paste("der(", temp[[(lengthcount * 2-1)]][1], ")", sep = ''),
                paste(
                  "der(",
                  temp[[(lengthcount * 2-1)]][1],
                  ";",
                  temp[[(lengthcount * 2-1)]][2],
                  ";",
                  temp[[(lengthcount * 2-1)]][1],
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
              ))
            
            currentvec <- getCytoBands(lengthcount, 2,temp,coln,derMods)
            
            tempins <-
              rbind(tempins, cbind(
                paste("der(", temp[[(lengthcount * 2-1)]][2], ")", sep = ''),
                paste(
                  "der(",
                  temp[[(lengthcount * 2-1)]][2],
                  ";",
                  temp[[(lengthcount * 2-1)]][2],
                  ")(",
                  currentvec[1],
                  ";",
                  currentvec[2],
                  ")",
                  sep = '',
                  collapse = ";"
                ),
                ""
              ))
            
            transloctable <- c(transloctable, list(tempins))
            names(transloctable)[length(transloctable)] <- inschrom
            
            ##dont take into account second derivative chromosome if there is only 1 chromosome where the insersion occures
            if(within)
            {
              transloctable[[length(transloctable)]]<-transloctable[[length(transloctable)]][-2,]
            }
            
          }
          
          
          
          if(!grepl("^\\+*ins\\(",Cyto_sample[coln]))
          {
            ##this is a factor, not a table 
            insDer<-as.matrix(transloctable[[grep(reginschrom,names(transloctable))]][grep(paste("der\\(",substr(Mainchr,0,nchar(Mainchr)-1),"\\)",sep=""),transloctable[[grep(reginschrom,names(transloctable))]][,1]),])[2]
            instemp <-
              strsplit(gsub("[\\(\\)]", "", regmatches(
                insDer, gregexpr("\\(.*?\\)",  insDer)
              )[[1]]), ";")
            temp[[(lengthcount*2-1)]]<-instemp[[1]]
            
            temp[[lengthcount*2]]<-instemp[[2]]
            
          }
          
        }
        ##check for labling long der again
        if (grepl("der\\([0-9;]+\\)\\(|rec\\([0-9;]+\\)\\(", derMods[(lengthcount*2-1)]) &
            !grepl("ider", Cyto_sample[coln]))
        {
          addBool <- paste(addBool, "LongDer", sep = "")
        }
        
        
        for (i in 1:length(temp[[(lengthcount * 2-1)]]))
        {
          
          if(!is.na(str_length(temp[[(lengthcount * 2)]][i])))
          {
            chr_table <-
              Cyto_ref_table[grep(paste(paste("chr", temp[[(lengthcount * 2-1)]][i], sep =
                                                ""), "$", sep = ""), Cyto_ref_table), ]
            ##put handling ? and ~ here
            ##put satilites and stuff here
            ##if long form
            ##deal with centromere stuff grep acen, take 2nd last
            ##der long form handling
            
            ##convert all - to ~ and cut out outer part of ~ if it is between 2 numbers
            ##need to do this better
            if (grepl("~|-", temp[[lengthcount * 2]][i]) )
            {
              ##maybe add function to make this less conservative
              ##right now it takes earlier one
              ##make this reg expression exclude ->
              temp[[lengthcount * 2]][i] <-
                gsub("-", "~", temp[[lengthcount * 2]][i])
              if (grepl("[pq][[:digit:]]+~", temp[[lengthcount * 2]][i]))
              {
                ##make sure this handles before and end
                temp[[lengthcount * 2]][i] <-
                  unlist(strsplit(temp[[lengthcount * 2]][i], "~"))[1]
              }
              
            }
            
            
            ##if long form
            
            if (any(grepl("::", temp[[lengthcount * 2]][i]) |
                    grepl("~>", temp[[lengthcount * 2]][i]))) {
              #parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
              ##ask tom about this one
              ##find p or q, take stuff before take stuff after, before is chromosomes after is positions, this will return 2 objects, must take into account
              ##parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
              ##ask tom about this one
              ##only splits first one
              longform_table <-
                strsplit(strsplit(temp[[lengthcount * 2]][i], "::")[[1]], "~>")
              ##take away any front loaded : 
              longform_table <-lapply(longform_table,function(x){gsub(':','',x)})
              
              in_table = data.frame()
              
              ##mark for dic, trc
              if(grepl("dic|trc",derMods[lengthcount]))
              {
                addBool<-paste("long",addBool,sep='')
              }
              
              ##get data for each read of something->something
              for (j in 1:length(longform_table))
              {
                stringdx <- str_locate_all(pattern = "p|q", longform_table[[j]])
                chr_name_long <-
                  substr(longform_table[[j]][1], 0, stringdx[[1]][1] - 1)
                positions <-
                  as.vector(cbind(
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
                  ))
                
                if (nchar(chr_name_long) != 0)
                {
                  chr_table_2 <-
                    Cyto_ref_table[grep(paste(paste(
                      "chr", chr_name_long, sep = ""
                    ), "$", sep = ""), Cyto_ref_table), ]
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
                  positions[grep("qter", positions)] <- chr_table_2[nrow(chr_table_2), 4]
                }
                ##be careful on centromeneter
                if (any(grepl("cen", positions)))
                {
                  ##likey will have to be careful about this one
                  positions[grep("cen", positions)] <-
                    chr_table_2[grep("acen", chr_table_2[, 5]), ][1, 4]
                }
                
                ##account for p10/q10
                ##double check naming convention for this one
                ##have to change the one for q
                positions[grep("p10", positions)] <-
                  chr_table_2[grep("acen", chr_table_2[, 5]), ][1, 4]
                positions[grep("q10", positions)] <-
                  chr_table_2[grep("acen", chr_table_2[, 5]), ][2, 4]
                
                ##make sure positions is in order to be processed correctly
                positions<-positionsorter(positions)
                
                ##put stuff in table
                positions_table <-
                  matrix(chr_table_2[grep(paste(positions, collapse = "|", sep = "|"),
                                          chr_table_2[, 4]), ], ncol = 5)
                
                if (is.vector(positions_table))
                {
                  positions_table <- t(positions_table)
                }
                
                in_table <-
                  rbind(
                    in_table,
                    cbind(
                      chr_table_2[1, 1],
                      positions_table[1, 2],
                      positions_table[nrow(positions_table), 3],
                      derMods[(lengthcount * 2-1)]
                    )
                  )
              }
              
              ##make table of all chromosomes used in this, then make table
              coord <- rbind(coord, in_table)
              ##make note where the breaks are
              ##excoord<-rbind(excoord,) ##constant exclusion
              
              
            } else {
              positions <-
                strsplit(gsub("q", ",q", gsub("p", ",p", temp[[lengthcount * 2]][i])), ",")[[1]][2:length(strsplit(gsub(
                  "q", ",q", gsub("p", ",p", temp[[lengthcount * 2]][i])
                ), ",")[[1]])]
              ##have to change the one to q
              positions[grep("p10", positions)] <-
                chr_table[grep("acen", chr_table[, 5]), ][1, 4]
              positions[grep("q10", positions)] <-
                chr_table[grep("acen", chr_table[, 5]), ][2, 4]
              ##check order, test this
              postions<-positionsorter(positions)
              
              positions_table <-
                matrix(chr_table[grep(paste(positions, collapse = "|"), chr_table[, 4]), ], ncol =
                         5)
              
              ################################################
              ##############################################
              #########################################
              ##if only one q or p , add on end point
              
              ##fix this so if p or q doesnt rely on order (p always being before q)
              if (length(positions) == 1)
              {
                if (any(grepl("p", positions)))
                {
                  if (is.vector(positions_table))
                  {
                    coord <-
                      rbind(coord,
                            cbind(chr_table[1, 1], "0", positions_table[3], derMods[lengthcount * 2 -
                                                                                      1]))
                    
                  } else{
                    coord <-
                      rbind(coord,
                            cbind(chr_table[1, 1], "0", positions_table[nrow(positions_table), 3], derMods[lengthcount *
                                                                                                             2 - 1]))
                  }
                }
                if (any(grepl("q", positions)))
                {
                  if (is.vector(positions_table))
                  {
                    coord <-
                      rbind(
                        coord,
                        cbind(
                          chr_table[1, 1],
                          positions_table[2],
                          ref_table[grep(paste(chr_table[1, 1], "$", sep = ""), ref_table), ][2],
                          derMods[(lengthcount * 2-1)]
                        )
                      )
                    
                  } else{
                    coord <-
                      rbind(
                        coord,
                        cbind(
                          chr_table[1, 1],
                          positions_table[1, 2],
                          ref_table[grep(paste(chr_table[1, 1], "$", sep = ""), ref_table), ][2],
                          derMods[(lengthcount * 2-1)]
                        )
                      )
                  }
                }
              } else{
                coord <-
                  rbind(
                    coord,
                    cbind(
                      chr_table[1, 1],
                      positions_table[1, 2],
                      positions_table[nrow(positions_table), 3],
                      derMods[(lengthcount * 2-1)]
                    )
                  )
              }
            }
            
            
          }
        }
      }
    }
    
    
    ##for loop every chr
    ##things are getting offset here
    ##+der ider duplicate
    if(grepl("ider",Cyto_sample[coln]) & grepl("t\\(",Cyto_sample[coln]))
    {
      coord<-coord[-1,]
    }
    
    if (nrow(coord) > 0)
    {
      ##case for vectors
      coord[, 1] <- as.character(coord[, 1])
      
      ##convert to numeric
      coord[, 2:3] <-
        apply(coord[, 2:3], 2, function(x) {
          as.numeric(as.character(x))
        })
      
      ##sort in order numerically
      coord[,2:3]<-t(apply(coord[,2:3],1,function(x){sort(x)}))
      
      ##if there are two translocations, fix coordinate overlap so only overlap counts,
      if(str_count(Cyto_sample[coln],"t\\(")>1)
      {
        ##find which chromosomes have more than one reading
        ##index of entries with chromosomes
        chromindex<-lapply(unique(coord[,1]),function(x){grep(x,coord[,1])})
        ##which index in chrom index have more than 2 readings
        relevantindex<-which(lapply(chromindex,length)>1)
        for(z in 1:length(relevantindex))
        {
          tempcoord<-coord[chromindex[[relevantindex[z]]],]
          
          
        }
        
      }
      ##fix this for different chromosomes
      ##excoord<-Cyto_ref_table[which(as.character(Cyto_ref_table[,1])=="chr16" & as.numeric(as.character(Cyto_ref_table[,2])) > 16800000 & as.numeric(as.character(Cyto_ref_table[,2]))<0),1:3]
      ##somehere here its grepping the wrong thing, seems speficical to vecotrs
      ##something wrong is happening here - 37 long  der()t()del() <- think somethieng here is over experessiong
      ##Chr<-unique(c(Mainchr,Allchr))
      for (p in 1:length(Mainchr))
      {
        ##need something for $ so chr 1 doesnt net chr 10
        curChr <- Mainchr[p]
        
        tempChr = coord[grep(paste(
          "chr",
          curChr,
          "$" ,
          sep = "",
          collapse = "|"
        ),
        coord[, 1]), ]
        tempChr = apply(tempChr, 2, as.character)
        
        
        if (is.vector(tempChr))
        {
          if (!identical(as.character(tempChr[2]), "0"))
          {
            excoord = rbind(excoord, cbind(tempChr[1], 0, tempChr[2], tempChr[4]))
          }
          
          if (!identical(as.character(tempChr[3]),
                         as.character(ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table), ][2])))
            excoord = rbind(excoord,
                            cbind(tempChr[1], tempChr[3], ref_table[grep(paste(tempChr[1], "$", sep =
                                                                                 ""), ref_table), ][2], tempChr[4]))
          
        }
        else{
          for (z in 1:nrow(tempChr))
          {
            if (z == 1 & tempChr[1, 2] != 0)
            {
              excoord = rbind(excoord,
                              cbind(tempChr[z, 1], 0, tempChr[z, 2], tempChr[z, 4]))
              if (nrow(tempChr) > 1)
              {
                excoord = rbind(excoord,
                                cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4]))
              }
              
            } else if (z == 1)
            {
              if (nrow(tempChr) > 1)
              {
                excoord = rbind(excoord,
                                cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4]))
              }
              
            } else if (z == nrow(tempChr) &
                       !identical(as.character(tempChr[z, 3]),
                                  as.character(ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep =
                                                                                 "")), ref_table), ][2]))) {
              excoord = rbind(excoord,
                              cbind(tempChr[z, 1], tempChr[z, 3], ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep =
                                                                                                      "")), ref_table), ][2], tempChr[z, 4]))
              
            } else if (z < nrow(tempChr) & z > 1) {
              excoord = rbind(excoord,
                              cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4]))
            }
            
            
          }
        }
      }
      
    }
    
    if(grepl("ider",Cyto_sample[coln]))
    {
      if(any(!grepl(Mainchr,coord[,1])))
      {
        coord<-rbind(coord,coord[grep(Mainchr,coord[,1],invert=T),])
      }
    }
    
    if (length(coord) > 0)
    {
      coord[, 2:3] <-
        apply(coord[, 2:3], 2, function(x) {
          as.numeric(as.character(x))
        })
      
      ##make it so coord only takes what is in the mainchr
      ## if(!any(grepl("ider\\(",coord[,4])))
      ##  {
      ##    coord<-coord[grep(paste(Mainchr,collapse="|"),coord[,1]),]
      ##  }
      coord[, 4] <- paste(addBool, coord[, 4], sep = "")
    }
    
    ##probably have to rework this because of inversions and insersions
    if (length(excoord) > 0)
    {
      ##excoord<-excoord[grep(paste(Mainchr,collapse="|"),excoord[,1]),]
      excoord[, 2:3] <-
        apply(excoord[, 2:3], 2, function(x) {
          as.numeric(as.character(x))
        })
      
      #as numeric as character
      ##delete any rows in which first number is larger than the second- inversions occur so dont do this
      bads<-vector(length=0)
      
      for(i in 1:nrow(excoord))
      {
        if(as.numeric(excoord[i,2])>as.numeric(excoord[i,3]))
        {
          bads<-c(bads,i)
        }
        
      }
      if(length(bads)>0)
      {
        excoord<-excoord[-bads,]
      }
      excoord[, 4] <- paste(addBool, excoord[, 4], sep = "")
      
    }
    
    ##make sure this does what i want it to
    ##excoord<-unique(excoord)
    ##coord<-unique(coord)
    
    ##do multi (X2) processing right now
    ##for over X2 times, its a gain 
    if (nrow(coord) > 0 && any(grepl("multi", coord[, 4])))
    {
      
      n <- multi - 1
      multimastercoord<-coord
      multimasterexcoord<-excoord
      for (f in 1:n)
      {
        if(f>1)
        {
          multitemp<-coord[grep("multi", coord[, 4]), ]
          multitemp[,4]<-paste("+",multitemp[,4],sep='')
          multimastercoord<- rbind(multimastercoord, multitemp)
          if (any(nrow(excoord) > 0 && grepl("multi", excoord[, 4])))
          {
            multitemp<-excoord[grep("multi", excoord[, 4]), ]
            multitemp[,4]<-paste("+",multitemp[,4],sep='')
            multimasterexcoord <- rbind(multimasterexcoord, multitemp)
          }
          
        }else{
          multimastercoord<- rbind(coord, coord[grep("multi", coord[, 4]), ])
          if (any(nrow(excoord) > 0 && grepl("multi", excoord[, 4])))
          {
            multimasterexcoord <- rbind(multimasterexcoord, excoord[grep("multi", excoord[, 4]), ])
            
          }
        }
      }
      coord<-multimastercoord
      excoord<-multimasterexcoord
    }
    
    ##have to fix derivative chromosomes  with two translocations , if they overlap in a chromosome , its the overlapping increment that maters
    
    
    
    
    
    
    listCoord <- list(coord, excoord, xmod, ymod, Mainchr, multi,transloctable,addtot)
    return(listCoord)
    
  }
  
  ##detects whether there is a + in special chromosomes to handle differently
  detectAdd <- function(temp_table_processed,
                        ex_table_processed)
  {
    ##not working properly because no longer table
    
    temp_table_processed[grep("\\+", temp_table_processed)] <- "+"
    
    ex_table_processed[grep("\\+", ex_table_processed)] <- "+"
    
    temp_table_processed[grep("^-", temp_table_processed)] <- "-"
    
    ex_table_processed[grep("^-", ex_table_processed)] <- "-"
    
    temp_table_processed[grep("\\+|^-", temp_table_processed, invert = T)] <-""
    
    ex_table_processed[grep("\\+|^-", ex_table_processed, invert = T)] <-""
    
    additionList <- list(temp_table_processed, ex_table_processed)
    return(additionList)
  }
  
  
  ##try catch between chromosomes
  colparse <- function(j, xmod, ymod, transloctable,addtot)
  {
    out <- try(parser(j, xmod, ymod, transloctable,addtot))
    
    if (inherits(out, "try-error")) {
      return(paste(geterrmessage(), "in", j, "field"))
    }
    return(out)
  }
  
  ##section to cancel out loss gains
  ### First code that takes two intervals (start, loss), where gain is first, that 
  ### returns a 3-column data frame with
  ### the intersection missing. 
  ### Can have 0, 1, or two rowsAlso has a column for type.
  
  mergeGL<-function(v1,v2){
    labs<-c("Gain", "Loss")
    v1<-as.vector(as.integer(as.numeric(v1[1:2])))
    v2<-as.vector(as.integer(as.numeric(v2[1:2])))
    out<-as.data.frame(matrix(nrow=0,ncol=3))
    colnames(out)<-c("Start","End","Type")
    out<-cbind(as.integer(as.numeric(out[,1])),as.integer(as.numeric(out[,2])),out[,3])
    if(v1[1]>v2[1]){
      labs<-labs[2:1]
      tmp<-v1
      v1<-v2
      v2<-tmp}
    
    if(v1[2]>=v2[2]){
      out<-data.frame(Start=c(v1[1],v2[2]+1), End=c(v2[1]-1, v1[2]), 
                      Type=rep(labs[1],2))} else{if(v1[2]>=v2[1]){
                        out<-data.frame(Start=c(v1[1],v1[2]+1), End=c(v2[1]-1, v2[2]), 
                                        Type=labs)} else{out<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                                                                         Type=labs)}}
    
    bads<-which(out[,1]>=out[,2])
    if(length(bads)>0){out<-out[-bads,]}
    origMat<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                        Type=labs)
    if(nrow(out)==nrow(origMat) && all(out==origMat))
    {
      out<-list(out,TRUE)
    }
    else
    {
      out<-list(out,FALSE)
    }
    return(out)}
  
  ### Now one that takes two two-column matrices of gains and losses
  ### (gains first) and returns a 3-column data frame of merged
  
  mergeGLmat<-function(G, L){
    i<-1
    ##G[,1:2]<-apply(G[,1:2],2,as.numeric)
    ##L[,1:2]<-apply(L[,1:2],2,as.numeric)
    origL<-L
    while(nrow(L)>0 & i<=nrow(G)){
      newL<-matrix(ncol=2, nrow=0)
      j<-1
      modified<-TRUE
      ##print(list(i,modified))
      while(j<=nrow(L) & modified){
        nxt<-mergeGL(G[i,], L[j,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        w<-which(nxt[,3]=="Loss")
        if(length(w)>0){newL<-rbind(newL, nxt[w,],if(nrow(L[-j])<j){L[-j,][j:nrow(L[-j,]),]})}else if(nrow(nxt)>0){newL<-rbind(newL,if(nrow(L[-j])<j){L[-j,][j:nrow(L[-j,]),]})}
        ##print(list(j,modified,newL))
        j<-j+1
      }
      L<-newL
      i<-i+1
    }
    
    
    i<-1
    while(nrow(G)>0 & i<=nrow(origL)){
      newG<-matrix(ncol=2, nrow=0)
      j<-1
      modified<-TRUE
      ##print(list(i,G))
      
      while(j<=nrow(G)& modified){
        nxt<-mergeGL(G[j,], origL[i,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        w<-which(nxt[,3]=="Gain")
        if(length(w)>0){newG<-rbind(newG, nxt[w,],if(nrow(G[-j])<j){G[-j,][j:nrow(G[-j,]),]})}else if(nrow(nxt)>0){newG<-rbind(newG,if(nrow(G[-j])<j){G[-j,][j:nrow(G[-j,]),]})}
        ##print(list(j,modified,newG))
        j<-j+1
      }
      G<-newG
      i<-i+1
    }
    
    out<-data.frame(Start=c(G[,1], L[,1]), End=c(G[,2],L[,2]), 
                    Type=c(rep("Gain", nrow(G)), rep("Loss", nrow(L))))
    return(out)}
  
  ### Now the whole thing:
  bigGLMerge<-function(M){
    
    ##ask what does this do 
    out<-M[-(1:nrow(M)),]
    
    chrs<-unique(M[,1])
    
    for(i in 1:length(chrs)){
      w<-which(M[,1]==chrs[i])
      Msub<-M[w,]
      
      wg<-which(Msub[,4]=="Gain")
      wl<-which(Msub[,4]=="Loss")
      
      if(length(wg)==0 | length(wl)==0){
        out<-rbind(out, Msub)} else{
          nxt<-mergeGLmat(Msub[wg,-1],Msub[wl,-1])
          if(nrow(nxt)>0)
          {
            nxt<-cbind(chrs[i],nxt)
            colnames(nxt)[1]<-"Chr"
          }
          out<-rbind(out, nxt)
        }
    }
    
    return(out)}
  
  ### This function takes a table M (a data frame) with columns:
  ### chromosome, start, stop, and type (Gain or Loss)
  ## and does proper merging
  
  ### First we need a function that takes two intervals (start, loss) that 
  ### returns NA if they don't overlap and merges them if they do.
  ### v1 and v2 are each two-vectors giving the interval  
  
  mergeInt<-function(v1,v2){
    if(is.na(prod(c(v1,v2)))){return(NA)}
    if(v1[1]>v2[1]){
      tmp<-v1
      v1<-v2
      v2<-tmp}
    
    if(v1[2]+1<v2[1]){return(NA)}
    
    return(c(v1[1], max(c(v1[2], v2[2]))))}
  
  ##mergeintL
  mergeIntL<-function(v1,v2){
    if(is.na(prod(c(v1,v2)))){return(NA)}
    if(v1[2]==v2[1]){
      return(NA)}
    
    if(v1[2]==v2[1]){return(NA)}
    
    return(c(v1[1], max(c(v1[2], v2[2]))))}
  
  
  ### Now a function that takes a 3-column table of intervals 
  ### (chr, start, stop) and
  ### returns a 3-column vector of intervals and c(NA,NA, NA)s where
  ### the top column is the merge of it and everything below
  ### that overlaps with it (these get replaced by c(NA,NA,NA)s)
  
  mergeBelow<-function(tbl){
    for(i in 2:nrow(tbl)){
      nxt<-mergeInt(t(tbl[1,2:3])[,1], t(tbl[i,2:3])[,1])
      if(!is.na(nxt[1])&(tbl[1,1]==tbl[i,1])){tbl[1,2:3]<-nxt
      tbl[i,]<-c(NA,NA,NA)}}
    return(tbl)} 
  
  ##for losses, merge consecutive things togehter
  mergeBelowL<-function(tbl){
    for(i in 2:nrow(tbl)){
      nxt<-mergeIntL(t(tbl[1,2:3])[,1], t(tbl[i,2:3])[,1])
      if(!is.na(nxt[1])&(tbl[1,3]==tbl[i,2])){tbl[1,2:3]<-nxt
      tbl[i,]<-c(NA,NA,NA)}}
    return(tbl)} 
  
  ##for merginign loss gain pairs (returns index )
  ##consolidatesimple<-function(gTable,lTable,i){ 
  ##think i should do this via recursion some how
  ##return(which(lTable[,1] == gTable[i,1] & lTable[,2] == gTable[i,2] & lTable[,3] == gTable[i,3]))
  
  ##}
  
  
  ### Now the whole deal
  
  mergeTable<-function(M){
    wg<-which(as.character(M[,4])=="Gain")
    wl<-which(as.character(M[,4])=="Loss")
    M<-M[c(wg,wl),]
    if(length(wg)>0 & length(wl)>0)
    {
      M<- bigGLMerge(M)
    }
    
    
    wg<-which(as.character(M[,4])=="Gain")
    wl<-which(as.character(M[,4])=="Loss")
    Mg<-M[wg,1:3]
    Ml<-M[wl,1:3]
    
    if(length(wg)>1){
      ##Mg<-M[wg,1:3]
      for(j in 1:(nrow(Mg)-1)){
        Mg[j:nrow(Mg),]<-mergeBelow(Mg[j:nrow(Mg),])}
      M[wg,1:3]<-Mg}
    
    if(length(wl)>1){
      Ml<-M[wl,1:3]
      for(j in 1:(nrow(Ml)-1)){
        Ml[j:nrow(Ml),]<-mergeBelowL(Ml[j:nrow(Ml),])}
      M[wl,1:3]<-Ml}
    
    bads<-which(!M[,4]%in%c("Loss", "Gain"))
    
    w<-setdiff(which(!is.na(M[,2])), bads)
    
    
    out<-M[w,]
    if(length(w)==1){out<-t(out)}
    
    return(out)}
  
  
  
  
  rowparse<-function(i){
    sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sample_table) <- c("Chr", "Start", "End", "Type")
    sorted_sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sorted_sample_table) <- c("Chr", "Start", "End", "Type")
    
    temp_table <-
      matrix(byrow = TRUE,
             nrow = 1,
             ncol = 4) ##temporary table for storage for mutations
    colnames(temp_table) <- c("Chr", "Start", "End", "Type")
    #  sample_table<-matrix(byrow=TRUE,nrow=0,ncol=4) ##running table for sample
    #  colnames(sample_table)<-list(paste(Con_data[i,1], "Chr"),"Start","End","Type")
    
    
    ##colnames(transloctable) <- c("Chr", "Start", "End", "Type")
    
    #  rownames(tempTable)<-rep(NULL,nrow(tempTable))
    normX = 0 ##normal number of X chormosomes
    normY = 0 ##normal number of Y chromosomes
    xcount = 0 ##counts number of x in 2nd slot
    ycount = 0##counts number of y in 2nd slot
    xadd = 0 ##counts if +X occures
    yadd = 0 ##counts if +Y occures
    xmod = 0 ##counts modifications that arent whole chromosome add/del for X
    ymod = 0##counds modifications that arent whole chromosome add/del for Y
    xdel = 0 ##counts if -X occures
    ydel = 0 ##counts if -Y occures
    addtot = 0 ##counts total "new chromosomes"
    deltot = 0 ##counts total complete chrom deletions
    modtot = 0 ##for idems only, counts modification chromosomes
    n = 1 ##ploidy count
    startcol = 2
    if (length(Cyto_sample) > 0 &
        grepl("[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
              Con_data[i, 2]))
    {
      ##initiate clonal evolution counter if ids is present
      if (grepl("^ids$", Cyto_sample[length(Cyto_sample)]))
      {
        ##dont count sex chromosomes here
        ##this has 46 chromosomes, as they get deleted or modified (without \\+), knocks em out of this tracker, stuff that remains is either gained or lost (likely gained)
        clone_chrom_tracker <- rep(1:22, 2)
      }
      
      if (grepl("[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
                Con_data[i, 2]))
      {
        ##set normal count for XY chromosomes
        ##think about what 46,X,+Y would mean
        if (any(grepl("Y", Cyto_sample)) &
            (sum(grepl("Y", Cyto_sample), na.rm = TRUE) > sum(grepl("\\+Y", Cyto_sample), na.rm =
                                                              TRUE)))
        {
          normX = 1
          normY = 1
        }
        else
        {
          normX = 2
        }
        
        ##cound number of XY in 2nd slot ##make sure 2nd slot is sex chromosomes, change this code later
        if (!grepl("^[XY]+", Cyto_sample[2]))
        {
          #make sure reg ecpression means what you want it to
        }  else {
          xcount = str_count(Cyto_sample[2], "X")
          ycount = str_count(Cyto_sample[2], "Y")
          ##x,y calculation
        }
      }
      
      ##check for ploidy levels and digit is over 2
      if (grepl("<[[:digit:]]n>", Cyto_sample[1]))
      {
        ##extract range of stuff before n
        n = as.numeric(strsplit(Cyto_sample[1], "<|n>")[[1]][2])
        n = n - 2
        if (n > 0)
        {
          for (p in 1:n)
          {
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Gain", nrow(ref_table)))
            ##colnames(temp_table)<-c( "Chr","Start","End","Type")
            temp_table <-
              temp_table[grep("chrY|chrX|chrM", temp_table[, 1], invert = T), ]
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
            
            ##not accoutning for sex chromosomes (will acount for later)
            addtot <- addtot + 22
            
          }
        }
        ## n= ##change ploidy count for x and y calculations
      }
      
      
      for (j in startcol:length(Cyto_sample))
      {
        temp_table <-
          matrix(byrow = TRUE,
                 nrow = 1,
                 ncol = 4) ##temporary table for storage for mutations
        #addition and deletions of entire chromosomes can be handeled easily and separate from the rest
        if (grepl(
          "mar|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*r\\(*[[:digit:]]*\\)*$|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*neo[[:digit:]]*$",
          Cyto_sample[j]
        ))
        {
          if (grepl("\\+", Cyto_sample[j]))
          {
            tem = 1
            ##figure out how to do this
            if (grepl("-|~", Cyto_sample[j]))
            {
              if (grepl("mar", Cyto_sample[j]))
              {
                Cyto_sample[j] <-
                  paste(unlist(strsplit(Cyto_sample[j], "~|-"))[1], "mar", sep = "")
              }
              if (grepl(
                "^\\+*([[:digit:]]((~|-)[[:digit:]])*)*r\\(*[[:digit:]]*\\)*$",
                Cyto_sample[j]
              ))
                Cyto_sample[j] <-
                  paste(unlist(strsplit(Cyto_sample[j], "~|-"))[1], "r", sep = "")
              
              if (grepl("neo", Cyto_sample[j]))
                Cyto_sample[j] <-
                  paste(unlist(strsplit(Cyto_sample[j], "~|-"))[1], "neo", sep = "")
              
            }
            if (grepl("\\+[[:digit:]]", Cyto_sample[j]))
            {
              tem <-
                as.numeric((strsplit(
                  strsplit(Cyto_sample[j], "mar|r\\(*[[:digit:]]*\\)*$|neo")[[1]][1],
                  "\\+"
                )[[1]][2]))
            }
            addtot <- addtot + tem
          }
        }else if (grepl("^\\+[[:digit:]]+c*$", Cyto_sample[j]) |
                  grepl("\\+X", Cyto_sample[j]) | grepl("\\+Y", Cyto_sample[j]))
        {
          cytoName <- gsub("c", "", substring(Cyto_sample[j], first = 2))
          chr_name <-
            ref_table[grep(paste("chr", as.character(cytoName), "$", sep = ""), ref_table), ]
          temp_table[1, 1] = chr_name[1]
          temp_table[1, 2] = "0"
          temp_table[1, 3] = chr_name[2]
          temp_table[1, 4] = "Gain"
          #if whole additions occur
          if (grepl("X", chr_name[1]))
          {
            xadd <- xadd + 1
          }
          
          if (grepl("Y", chr_name[1]))
          {
            yadd <- yadd + 1
          }
          ##if(!grepl("X|Y",chr_name[1]))
          ##{
          addtot <- addtot + 1
          ##}
        } else if (grepl("^-[[:digit:]]+c*$", Cyto_sample[j]) |
                   grepl("-X", Cyto_sample[j]) | grepl("-Y", Cyto_sample[j]))
        {
          cytoName <- gsub("c", "", substring(Cyto_sample[j], first = 2))
          chr_name <-
            ref_table[grep(paste("chr", as.character(cytoName), "$", sep = ""), ref_table), ]
          temp_table[1, 1] = chr_name[1]
          temp_table[1, 2] = "0"
          temp_table[1, 3] = chr_name[2]
          temp_table[1, 4] = "Loss"
          
          ##if deletions occur in X or Y, up count for the respective mutation
          if (grepl("X", chr_name[1]))
          {
            xdel <- xdel + 1
          }
          
          if (grepl("Y", chr_name[1]))
          {
            ydel <- ydel + 1
          }
          
          ##if(!grepl("X|Y",chr_name[1]))
          ##{
          deltot <- deltot + 1
          ##}        
        } else
        {
          ##for all other cases call this first
          ##check for x modifications in parser
          inc_table <- colparse(j, xmod, ymod, transloctable,addtot)
          if (is.character(inc_table))
          {
            Dump_table <- rbind(Dump_table, c(Con_data[i, ], inc_table))
          }
          else if (is.null(inc_table)) {
            
          }
          else{
            temp_table <- inc_table[[1]]
            ex_table <- inc_table[[2]]
            xmod <- inc_table[[3]]
            ymod <- inc_table[[4]]
            Mainchr <- inc_table[[5]]
            multi <- inc_table[[6]]
            transloctable<-inc_table[[7]]
            addtot<-inc_table[[8]]
            if (nrow(temp_table) > 0)
            {
              ##special cases (if, else if, else if ,else if, else )
              ##ider, i , ins, dic etc
              ##else, grep stuff on 4th col to detirmine addition or deletions, dont include translations add etc
              ##deal with + values differently
              
              ###if it goes through special cases, flip mod = true, excoord does not combine at end
              mod = FALSE
              temp_table[, 4] <- as.character(temp_table[, 4])
              colnames(temp_table) <- c("Chr", "Start", "End", "Type")
              
              
              if (length(ex_table) > 0)
              {
                colnames(ex_table) <- c("Chr", "Start", "End", "Type")
                ex_table[, 4] <- as.character(ex_table[, 4])
                ex_table[, 2:3] <-
                  apply(ex_table[, 2:3], 2, function(x) {
                    as.numeric(as.character(x))
                  })
              }
              temp_table[, 2:3] <-
                apply(temp_table[, 2:3], 2, function(x) {
                  as.numeric(as.character(x))
                })
              
              ## see if its an addition
              
              additionTable <- detectAdd(temp_table[, 4], if(nrow(ex_table)>0){ex_table[, 4]}else{NULL})
              ##if its not a del ##think about this one for all cases
              
              ##need to improve counter , maybe put it way upstream in each j so if it finds +, adds, then take it out of everything else
              ##probbbly doesnt take into account x3 (multi -2)
              if (grepl("\\+", Cyto_sample[j]) ) ##&&!any(grepl("X|Y", Mainchr))
              {
                addtot <- (1 * multi) + addtot
              }
              else if(multi>2){
                addtot<-addtot+multi-2
              }
              ##temp_table[grep("+\(.\)",ex_table[,4]),4]<-"Gain"
              ##ex_table[,4]<-"Gain"
              ##rbind(temp_table,ex_table)
              ##incorporate main chr in these
              ## do del and stuff up here , but keep+__ sign somehow
              
              ##count chromosome loss properly here for special cases
              ##sex chromosome count check how it interacts with xmod ymod
              if (!all(grepl("\\+", temp_table[, 4])) &
                  any(!grepl("^t\\(|^ins\\(|^inv\\(", temp_table[, 4])))
              {
                ##takes mainchr off of list of chromosomes because it has been acccounted for (idems only)
                if (grepl("^ids$", Cyto_sample[length(Cyto_sample)]))
                {
                  if(length(Mainchr[-(grep("X|Y",Mainchr))])>0)
                  {
                    clone_chrom_tracker <-
                      clone_chrom_tracker[-1 * sapply(gsub("\\$", "", Mainchr), function(x) {if(!grepl("X|Y",x))
                        which.max(clone_chrom_tracker == as.numeric(x))
                      })]
                  }
                }
                
                ##check how this handles X chromosomes
                ## if more than 2 chromosomes involved and its not a translocation an insertion or an inversion, something is "deleted"
                if (any(grepl("-[:alpha:]", temp_table[, 4]) &
                        !grepl("X\\$|Y\\$", Mainchr)))
                {
                  ##deltot <- deltot + (1 * multi)
                }
                else
                {
                  if(length(Mainchr)!=length(grep("X|Y", Mainchr)))
                  {
                    deltot <- deltot + (length(grep("X|Y", Mainchr, invert = T)) - 1) * multi
                  }
                  
                }
                
              }
              
              ##insertions affect count if in derivative chromosome or the like
              
              ##take any translocations and insertions and stores for later
              ##if (any(grepl("(t|ins)\\(", ex_table[, 4])))
              ##{
              ##transloctable <-
              ##rbind(transloctable, ex_table[grep("(t|ins)\\(", ex_table[, 4]), ])
              ##}
              
              if (any(grepl("LongDer", temp_table[, 4])))
              {
                mod = TRUE
                temp_table[intersect(grep("LongDer", temp_table[, 4]),grep(paste(Mainchr,collapse="|",sep=''),temp_table[,1],invert=T)), 4]<-"Gain"
                if(nrow(ex_table)>0)
                {
                  ex_table[grep("LongDer", ex_table[, 4]), 4] <- "Loss"
                }
                ##temp_table[, 4] <-
                ##paste(additionTable[[1]], temp_table[, 4], sep = "")
                ##ex_table[, 4] <-
                ##paste(additionTable[[2]], ex_table[, 4], sep = "")
                temp_table <- rbind(temp_table, ex_table)
              }
              
              ##if(any(grepl("\\+der\\(.*",temp_table[,4])))
              ##{
              ##mod=TRUE
              ##temp_table[grep("\\+der\\(.*",temp_table[,4]),4]<-"Gain"
              ##think about how this is handling issue
              ###really check on this
              ##ex_table[grep("\\+der\\(.*",ex_table[,4]),4]<-"Gain"
              ##temp_table[,4]<-paste(additionTable[[1]],temp_table[,4],sep="")
              ##ex_table[,4]<-paste(additionTable[[2]],ex_table[,4],sep="")
              ##temp_table<-rbind(temp_table,ex_table)
              ##}
              if(any(grepl("(t|ins)\\(",temp_table[,4])&!grepl("^(t|ins)\\(",temp_table[,4])))
              {
                mod = TRUE
                temp_table[intersect(intersect(
                  grep(
                    paste("chr", Mainchr, sep = "",collapse="|"),
                    temp_table[, 1],
                    invert = T
                  ),
                  grep(
                    "(t|ins)\\(",
                    temp_table[, 4]
                  )), grep("^(t|ins)\\(",temp_table[,4],invert=T)
                ), 4] <- "Gain"
                if(nrow(ex_table)>0)
                {
                  ex_table[intersect(intersect(
                    grep(paste("chr", Mainchr, sep = "",collapse="|"), ex_table[, 1]),
                    grep(
                      "t\\(",
                      ex_table[, 4]
                    )), grep("^t\\(",ex_table[,4],invert=T)
                  ), 4] <- "Loss"
                }
                ##temp_table[, 4] <-
                ##  paste(additionTable[[1]], temp_table[, 4], sep = "")
                ##ex_table[, 4] <-
                ##  paste(additionTable[[2]], ex_table[, 4], sep = "")
                ##temp_table <- rbind(temp_table, ex_table)
              }
              
              
              if (any(grepl("del", temp_table[, 4])&!grepl("^del", temp_table[, 4])&!grepl("multi[[:digit:]]*del",temp_table[,4])))
              {
                ##deals with if there is one deletion, will deal w multiple cases later
                ##this is a factor problem
                deleted <- temp_table[intersect(grep("del", temp_table[-1, 4]),grep("^del", temp_table[-1, 4],invert=T)), ]
                if(length(deleted)>0 && nrow(deleted)>0)
                {
                  
                  deleted[2:3] <-
                    apply(deleted[ 2:3], 2, function(x) {
                      as.numeric(as.character(x))
                    })
                  
                  
                  ##chr_table<-Cyto_ref_table[grep(paste(paste("chr",temp_table[1,1],sep=""),"$",sep= ""),Cyto_ref_table),]
                  ##if its identical to the end, truncate end
                  ## if its idential to the beggining, truncate beginning
                  ##else, split into lines
                  ##maybe not a vector when changed
                  ##fix this part
                  if (is.vector(deleted))
                  {
                    if (identical(as.numeric(deleted[2]),
                                  as.numeric(temp_table[1, 2])))
                    {
                      temp_table[1, 2] <- deleted[3]
                    } else if (identical(as.numeric(deleted[3]),
                                         as.numeric(temp_table[1, 3])))
                    {
                      temp_table[1, 3] <- deleted[2]
                    } else
                    {
                      ##rearranging new data
                      ##double check if this can be an else
                      temp_table <- matrix(nrow = 0, ncol = 4)
                      temp_table <- rbind(temp_table, deleted)
                      temp_table[1, 3] <- deleted[2]
                      temp_table[2, 2] <- deleted[3]
                      
                    }
                    
                    
                  } else
                  {
                    ##thinkg about what your doing here
                    for(k in 1:nrow(deleted))
                    {
                      ##chromosome of deletion
                      tempChr<-paste(deleted[k,1],"$",sep='')
                      tempChrTableIndex<-intersect(grep(tempChr,temp_table[,1]),grep("del",temp_table[,4],invert=T))  
                      if(length(tempChrTableIndex)>0)
                      {
                        
                        for(l in 1:length(tempChrTableIndex))
                        {
                          if (identical(as.numeric(deleted[k, 2]),
                                        as.numeric(temp_table[tempChrTableIndex[1], 2])))
                          {
                            temp_table[1, 2] <- deleted[k,3]
                          } else if (identical(as.numeric(deleted[k,3]),
                                               as.numeric(temp_table[1, 3])))
                          {
                            temp_table[1, 3] <- deleted[k, 2]
                          } else
                          {
                            ##rearranging new data
                            ##double check if this can be an else
                            temp_table <- matrix(nrow = 0, ncol = 4)
                            temp_table <- rbind(temp_table, deleted)
                            temp_table[1, 3] <- deleted[k, 2]
                            temp_table[2, 2] <- deleted[k, 3]
                            
                          }
                        }
                      }
                    }
                    
                    
                  }
                }
              }
              
              
              
              if (any(grepl(
                "^\\++(der|rec)\\(.*|^(der|rec)\\(.*",
                temp_table[, 4]
              )))
              {
                mod = TRUE
                #temp_table[intersect(
                #  grep(
                #    paste("chr", Mainchr, sep = ""),
                #    temp_table[, 1],
                #    invert = T
                #  ),
                #  grep(
                #    "^\\++(der|rec)\\(.*|^(der|rec)\\(.*",
                #    temp_table[, 4]
                #  )
                #), 4] <- "Gain"
                #ex_table[intersect(
                #  grep(paste("chr", Mainchr, sep = ""), ex_table[, 1]),
                #  grep(
                #    "^\\++(der|rec)\\(.*|^(der|rec)\\(.*",
                #    ex_table[, 4]
                #  )
                #), 4] <- "Loss"
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
                
              }
              
              if (any(grepl("i\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                temp_table[grep("i\\(.*", temp_table[, 4]), 4] <- "Gain"
                if(nrow(ex_table)>0)
                {
                  ex_table[grep("i\\(.*", ex_table[, 4]), 4] <- "Loss"
                }
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
              }
              ##get whats included , dup, then rest on chromosome is deletion
              
              if (any(grepl("ider\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                ##have to be able to handle these cases
                ##do inclusion and exclusion based on deleted case
                
                temp_table[grep("ider\\(.*", temp_table[, 4]), 4] <-
                  "Loss"
                if(nrow(ex_table)>0)
                {
                  ex_table[grep("ider\\(.*", ex_table[, 4]), 4] <- "Gain"
                }
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
                
              }
              ##fix this stuff dic
              
              if (any(grepl("idic\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                ##need to either break it up or replace,
                ##think about this one too
                ##deleted <-
                ##ex_table[grep("idic\\(.*&del\\(", ex_table[, 4]), 4]
                
                temp_table[grep("idic\\(.*", temp_table[, 4]), 4] <- "Loss"
                if(nrow(ex_table)>0)
                {
                  ex_table[grep("idic\\(.*", ex_table[, 4]), 4] <- "Gain"
                }
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
              }
              
              ##figour out why you did this 
              ##make sure this is reversed for long form
              if (any(grepl("dic\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                if(any(grepl("long",temp_table[,4])))
                {
                  temp_table[, 4] <-
                    paste(additionTable[[1]], temp_table[, 4], sep = "")
                  if(nrow(ex_table)>0)
                  {
                    ex_table[grep("dic\\(.*", ex_table[, 4]) ,4]<-"Loss"
                    ex_table[, 4] <-
                      paste(additionTable[[2]], ex_table[, 4], sep = "")
                  }
                  
                }else
                {
                  temp_table[grep("dic\\(.*", temp_table[, 4]), 4] <- "Loss"
                  
                  temp_table[, 4] <-
                    paste(additionTable[[1]], temp_table[, 4], sep = "")
                  if(nrow(ex_table)>0)
                  {
                    ex_table[, 4] <-
                      paste(additionTable[[2]], ex_table[, 4], sep = "")
                  }
                }
                temp_table <- rbind(temp_table, ex_table)
              }
              
              ##make sure this is reversed for long form 
              if (any(grepl("trc\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                if(any(grepl("long",temp_table[,4])))
                {
                  ex_table[intersect(grep("trc\\(.*", ex_table[, 4]),grep(paste(Mainchr[1],Mainchr[3],sep='|'),ex_table[,4])),4]<-"Loss"
                  
                }else
                {
                  temp_table[intersect(grep("trc\\(.*", temp_table[, 4]),grep(paste(Mainchr[1],Mainchr[3],sep='|'),temp_table[,4])),4]<-"Loss"
                }
                
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[intersect(grep("trc\\(.*", ex_table[, 4]),grep(Mainchr[2],ex_table[,4])), 4] <- "Loss"
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                
                temp_table <- rbind(temp_table, ex_table)
              }
              
              if (any(grepl("rob\\(", temp_table[, 4])))
              {
                mod = TRUE
                if(nrow(ex_table)>0)
                {
                  ex_table[grep("rob\\(", ex_table[, 4]), 4] <- "Loss"
                }
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0){
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
                
              }
              
              if (any(grepl("^r\\(.*|^\\+r\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                if(nrow(ex_table)>0)
                {
                  
                  ex_table[grep("r\\(.*", ex_table[, 4]), 4] <- "Loss"
                }
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
                
              }
              
              if (any(grepl("trp\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table<-rbind(temp_table,temp_table[grep("trp\\(.*", temp_table[, 4]), ]) 
                
              }
              
              if (any(grepl("qdp\\(.*", temp_table[, 4])))
              {
                mod = TRUE
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0){
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table<-rbind(temp_table,temp_table[grep("qdp\\(.*", temp_table[, 4]), ],temp_table[grep("qdp\\(.*", temp_table[, 4]), ]) 
              }
              
              
              
              ##add code here for +gain +loss (+gain=gain, +loss ="")
              if (any(grepl("del", temp_table[, 4])))
              {
                mod = TRUE
                temp_table[grep("del", temp_table[, 4]), 4] <-
                  "Loss"
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table) > 0)
                {
                  if (is.vector(ex_table))
                  {
                    ex_table[, 4] <- ""
                  } else
                  {
                    ex_table[, 4] <- ""
                  }
                  
                  
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table<-rbind(temp_table,ex_table)
              } 
              
              
              ##keep X and Y chromosomes here for ex table for later processing, may have to modify later for autosomes, only if modifications aren't triggered
              if (mod == FALSE)
              {
                
                if (is.vector(ex_table))
                {
                  ex_table[ 4] <- ""
                } else
                {
                  if(nrow(ex_table)>0)
                  {
                    ex_table[, 4] <- ""
                  }
                }
                
                temp_table[, 4] <-
                  paste(additionTable[[1]], temp_table[, 4], sep = "")
                if(nrow(ex_table)>0)
                {
                  ex_table[, 4] <-
                    paste(additionTable[[2]], ex_table[, 4], sep = "")
                }
                temp_table <- rbind(temp_table, ex_table)
              }
              
              
              
              ##handel stuff for minus as well
              temp_table[grep("^-Gain|^-$", temp_table[, 4]), 4] <-
                "Loss"
              temp_table[grep("\\+Gain", temp_table[, 4]), 4] <- "Gain"
              temp_table[grep("\\+Loss", temp_table[, 4]), 4] <- ""
              temp_table[grep("\\+", temp_table[, 4]), 4] <- "Gain"
              ##check other permutations of this/ dont think insertions belong here
              temp_table[grep("dup|qdp|tan|trp|\\+$", temp_table[, 4]), 4] <-
                "Gain"
            }
          }
        }
        
        ##combine running file and temp new file together
        if (length(temp_table) != 0)
        {
          temp_table <- apply(temp_table, 2, as.character)
          sample_table <- rbind(sample_table, temp_table)
        }
      }
      
      
      ##think about how final table works, make another structure for it
      
      ##stores temporary sex chromosomes (x or y depending on situation)
      temp_sex_chrom <- matrix(nrow = 0, ncol = 4)
      ##stores everything else to remerge
      temp_auto_chrom <- matrix(nrow = 0, ncol = 4)
      ##grand total x y choromosome counts
      ##check for constitutional anomolies in X and Y chromosomes
      ##check for XXY combos
      ##something here is not counting correctly
      
      if (xcount + xadd + xdel + xmod > normX)
      {
        if (xadd + xmod - xdel + xcount > normX &
            length(sample_table) > 0 && grepl("chrX", sample_table))
        {
          temp_sex_chrom <- sample_table[grep("chrX", sample_table[, 1]), ]
          temp_auto_chrom <-
            sample_table[grep("chrX", sample_table[, 1], invert = T), ]
          
          ##modify  so X gains->gains, loss->nothing, otherwise ->Gain
          if (is.vector(temp_sex_chrom))
          {
            ##think hard about where in here would need to be counting stuff
            if (!grepl("Gain", temp_sex_chrom[4]) &
                !grepl("Loss", temp_sex_chrom[4]))
            {
              temp_sex_chrom[intersect(
                grep("Gain", temp_sex_chrom[4], invert = T),
                grep("Loss", temp_sex_chrom[4], invert = T)
              )] <- "Gain"
            }
            if (grepl("Loss", temp_sex_chrom[4]))
            {
              temp_sex_chrom[4] <- ""
            }
            
          } else{
            temp_sex_chrom[intersect(
              grep("Gain", temp_sex_chrom[, 4], invert = T),
              grep("Loss", temp_sex_chrom[, 4], invert = T)
            ), 4] <- "Gain"
            temp_sex_chrom[grep("Loss", temp_sex_chrom[, 4]), 4] <- ""
            
          }
          sample_table <- rbind(temp_auto_chrom, temp_sex_chrom)
          addtot <- xmod + addtot
        } else if (xcount > normX)
        {
          for (z in 1:(xcount - normX)) {
            chr_name <- ref_table[grep("chrX", ref_table), ]
            sample_table <-
              rbind(sample_table,
                    cbind(chr_name[1], "0", chr_name[2], "Gain"))
            addtot <- addtot + 1
            xadd <- xadd + 1
          }
        }
        
        
      } else if (xcount + xadd + xdel + xmod < normX)
      {
        ##add -X to table
        for (z in 1:(normX - (xcount + xadd + xdel + xmod)))
        {
          chr_name <- ref_table[grep("chrX", ref_table), ]
          sample_table <-
            rbind(sample_table, cbind(chr_name[1], "0", chr_name[2], "Loss"))
          deltot <- deltot + 1
          xdel <- xdel + 1
        }
        
      }
      
      ##for y choromosome
      
      if (ycount + yadd + ydel + ymod > normY)
      {
        if (ycount + yadd - ydel + ymod > normY &
            length(sample_table) > 0 && grepl("chrY", sample_table))
        {
          temp_sex_chrom <- sample_table[grep("chrY", sample_table[, 1]), ]
          temp_auto_chrom <-
            sample_table[grep("chrY", sample_table[, 1], invert = T), ]
          
          ##modify  so Y gains->gains, loss->nothing, otherwise ->Gain
          ##but if only one of 2 are gains, one needs to stay the same whi;e the other stays constant
          if (is.vector(temp_sex_chrom))
          {
            if (!grepl("Gain", temp_sex_chrom[4]) &
                !grepl("Loss", temp_sex_chrom[4]))
            {
              temp_sex_chrom[intersect(
                grep("Gain", temp_sex_chrom[4], invert = T),
                grep("Loss", temp_sex_chrom[4], invert = T)
              )] <- "Gain"
            }
            if (grepl("Loss", temp_sex_chrom[4]))
            {
              temp_sex_chrom[4] <- ""
            }
            
          } else{
            temp_sex_chrom[intersect(
              grep("Gain", temp_sex_chrom[, 4], invert = T),
              grep("Loss", temp_sex_chrom[, 4], invert = T)
            ), 4] <- "Gain"
            temp_sex_chrom[grep("Loss", temp_sex_chrom[, 4]), 4] <- ""
            
          }
          ###one of these is a duplication, the other is not.. how to fix this
          addtot <- ymod + addtot
          sample_table <- rbind(temp_auto_chrom, temp_sex_chrom)
          
        }
        ##add Y chromosome if still not accounted for
        if (ycount > normY)
        {
          for (z in 1:(ycount - normY)) {
            chr_name <- ref_table[grep("chrY", ref_table), ]
            sample_table <-
              rbind(sample_table,
                    cbind(chr_name[1], "0", chr_name[2], "Gain"))
            addtot <- addtot + 1
            yadd <- yadd + 1
          }
        }
      } else if (ycount + yadd + ydel + ymod < normY)
      {
        ##add -Y to table
        for (z in 1:(normY - (ycount + yadd + ydel + ymod)))
        {
          chr_name <- ref_table[grep("chrY", ref_table), ]
          sample_table <-
            rbind(sample_table, cbind(chr_name[1], "0", chr_name[2], "Loss"))
          deltot <- deltot + 1
          ydel <- ydel + 1
        }
      }
    }
    
    ##total count, ensure no wild stuff happened
    ##have to account for ranges,
    ##have to count markers
    val = strsplit(Cyto_sample[1], "<")[[1]] ##numer of chromosomes indicated in first value
    val = paste(strsplit(val, "[[:alpha:]]+")[[1]],
                sep = "",
                collapse = "")
    if (grepl("~|-", val))
    {
      val = unlist(strsplit(val, "~|-"))
      if (all(!is.na(as.numeric(val))))
      {
        val = seq(from = val[1],
                  to = val[2],
                  by = 1)
      }
      else{
        val = val[2]
      }
      
    }
    
    val = as.numeric(val)
    
    ##take x and y completelu out of this 
    val = val -xcount -ycount-xadd-yadd+xdel+ydel
    if (any(is.na(val)))
    {
      Dump_table <-
        rbind(Dump_table, c(as.vector(Con_data[i, ]), "unclear chrom number"))
    } else
    {
      idealval = 46 + addtot - deltot- xcount - ycount -(xadd+yadd-xdel-ydel)
      diffval = val - idealval
      val_remainder = diffval %% 22
      val_remainder_2=diffval %% 23
      val_divider = diffval / 22 ## for stuff over diploidy
      val_divider_2=diffval/23
      if (!is.null(val_remainder))
      {
        if (any(val_remainder == 0))
        {
          val_remainder <- val_remainder[grep(TRUE, val_remainder == 0)]
          val_divider <- val_divider[grep(TRUE, val_remainder == 0)]
        }else if(any(val_remainder_2==0)){
          val_remainder<- val_remainder_2[grep(TRUE, val_remainder_2 == 0)]
          val_divider <- val_divider_2[grep(TRUE, val_remainder_2 == 0)]
        }else{
          val_remainder = val_remainder[1]
          val_divider = val_divider[1]
        }
        
        print(c(val, idealval, val_remainder, val_divider))
        
        
        if (val_remainder == 0 )
        {
          
          
          if (val_divider > 0)
          {
            for (f in 1:(val_divider))
            {
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table <-
                temp_table[grep("chrY|chrX|chrM", temp_table[, 1], invert = T), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              
            }
            ##if x and y need to be included in the table for polyploidys
            if(xcount+ycount+xmod+ymod<(val_divider+2))
            {
              for(f in 1:(val_divider+2-(xcount+ycount+xmod+ymod)))
              {
                if(ycount+ymod>0)
                {
                  temp_table <-
                    data.frame(ref_table[, 1],
                               rep(0, nrow(ref_table)),
                               ref_table[, 2],
                               rep("Gain", nrow(ref_table)))
                  temp_table <-
                    temp_table <-
                    temp_table[grep("chrY|chrX", temp_table[, 1]), ]
                  
                  temp_table[, 4] <- as.character(temp_table[, 4])
                  temp_table <- as.matrix(temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  sample_table <- rbind(sample_table, temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  
                }else{
                  
                  
                  temp_table <-
                    rbind(temp_table[grep("chrX", temp_table[, 1]), ],temp_table[grep("chrX", temp_table[, 1]), ])
                  
                  temp_table[, 4] <- as.character(temp_table[, 4])
                  temp_table <- as.matrix(temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  sample_table <- rbind(sample_table, temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  
                }
              }
            }
            print(c("val_div>0", Cyto_sample))
          }
          
          if (val_divider < 0 )
          {
            if (val_divider == -1)
            {
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Loss", nrow(ref_table)))
              temp_table <-
                temp_table[grep("chrY|chrX|chrM", temp_table[, 1], invert = T), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              
              
              
              
            }
            
            print(c("val_div<0", Cyto_sample))
          }
          
        } else if (grepl("^ids$", Cyto_sample[length(Cyto_sample)]))
        {
          ##if unaccounted chromosom number == diff val, add them
          if (any(diffval / length(clone_chrom_tracker) > 0 &
                  diffval %% length(clone_chrom_tracker) == 0))
          {
            
            diffval<-diffval[intersect(grep(TRUE,diffval / length(clone_chrom_tracker) > 0),grep(TRUE,diffval %% length(clone_chrom_tracker) == 0))]
            
            
            for (k in 1:diffval / length(clone_chrom_tracker))
            {
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table <-
                temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                                temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            }
            
            ##for when x y count need sto be added back
            if(xcount+ycount+xmod+ymod<(val_divider+2))
            {
              for(f in 1:(val_divider+2-(xcount+ycount+xmod+ymod)))
              {
                if(ycount+ymod>0)
                {
                  temp_table <-
                    data.frame(ref_table[, 1],
                               rep(0, nrow(ref_table)),
                               ref_table[, 2],
                               rep("Gain", nrow(ref_table)))
                  temp_table <-
                    temp_table <-
                    temp_table[grep("chrY|chrX", temp_table[, 1]), ]
                  
                  temp_table[, 4] <- as.character(temp_table[, 4])
                  temp_table <- as.matrix(temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  sample_table <- rbind(sample_table, temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  
                }else{
                  
                  
                  temp_table <-
                    rbind(temp_table[grep("chrX", temp_table[, 1]), ],temp_table[grep("chrX", temp_table[, 1]), ])
                  
                  temp_table[, 4] <- as.character(temp_table[, 4])
                  temp_table <- as.matrix(temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  sample_table <- rbind(sample_table, temp_table)
                  sample_table[, 4] <- as.character(sample_table[, 4])
                  
                }
              }
              ##if its negative all these were lost
              if (any(diffval / length(clone_chrom_tracker) == -1))
              {
                temp_table <-
                  data.frame(ref_table[, 1],
                             rep(0, nrow(ref_table)),
                             ref_table[, 2],
                             rep("Loss", nrow(ref_table)))
                temp_table <-
                  temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                                  temp_table[, 1]), ]
                
                temp_table[, 4] <- as.character(temp_table[, 4])
                temp_table <- as.matrix(temp_table)
                sample_table[, 4] <- as.character(sample_table[, 4])
                sample_table <- rbind(sample_table, temp_table)
                sample_table[, 4] <- as.character(sample_table[, 4])
              }
              
            }
          }
        }else{
          ##put in dump table
          Dump_table <-
            rbind(Dump_table,
                  c(
                    as.vector(Con_data[i, ]),
                    "some chromosomes unaccounted for"
                  ))
          print(c("unaccounted", Cyto_sample))
        }
      }
    }
    
    ##something is wrong here 
    sample_table <- sample_table[rowSums(!is.na(sample_table)) > 0, ]
    sample_table<-as.data.frame(sample_table,row.names = FALSE)
    if (is.vector(sample_table))
    {
      sample_table[1]<-as.character(sample_table[1])
      sample_table[2]<-as.integer(as.numeric(as.character(sample_table[2])))
      sample_table[3]<-as.integer(as.numeric(as.character(sample_table[3])))
      sample_table[4]<-as.character(sample_table[4])
      sample_table <- t(sample_table)
      sorted_sample_table<-sample_table
    }else if(ncol(sample_table)==1){
      sample_table <- t(sample_table)
      sample_table[,1]<-as.character(sample_table[,1])
      sample_table[,2]<-as.integer(as.numeric(as.character(sample_table[,2])))
      sample_table[,3]<-as.integer(as.numeric(as.character(sample_table[,3])))
      sample_table[,4]<-as.character(sample_table[,4])
      sorted_sample_table<-sample_table
    }else if(nrow(sample_table)>1){
      sample_table[,1]<-as.character(sample_table[,1])
      sample_table[,2]<-as.integer(as.numeric(as.character(sample_table[,2])))
      sample_table[,3]<-as.integer(as.numeric(as.character(sample_table[,3])))
      sample_table[,4]<-as.character(sample_table[,4])
      ## eliminate duplicates and gain/loss with same coordinates
      sorted_sample_table<-mergeTable(sample_table)
    }
    
    ##sorted_sample_table<-sample_table
    
    ##correct format for sample table
    ##something is very wrong here
    if (is.vector(sorted_sample_table))
    {
      sorted_sample_table[1]<-as.character(sorted_sample_table[1])
      sorted_sample_table[2]<-as.numeric(as.character(sorted_sample_table[2]))
      sorted_sample_table[3]<-as.numeric(as.character(sorted_sample_table[3]))
      sorted_sample_table[4]<-as.character(sorted_sample_table[4])
      sorted_sample_table <- t(sorted_sample_table)
    }else if(ncol(sorted_sample_table)==1){
      sorted_sample_table <- t(sorted_sample_table)
      sorted_sample_table<-as.data.frame(sorted_sample_table)
      sorted_sample_table[,1]<-as.character(sorted_sample_table[,1])
      sorted_sample_table[,2]<-as.numeric(as.character(sorted_sample_table[,2]))
      sorted_sample_table[,3]<-as.numeric(as.character(sorted_sample_table[,3]))
      sorted_sample_table[,4]<-as.character(sorted_sample_table[,4])
    }else if(nrow(sorted_sample_table)>1){
      sorted_sample_table<-as.data.frame(sorted_sample_table)
      sorted_sample_table[,1]<-as.character(sorted_sample_table[,1])
      sorted_sample_table[,2]<-as.numeric(as.character(sorted_sample_table[,2]))
      sorted_sample_table[,3]<-as.numeric(as.character(sorted_sample_table[,3]))
      sorted_sample_table[,4]<-as.character(sorted_sample_table[,4])
    }
    
    return(list(sorted_sample_table,Dump_table,transloctable))
  }
  ##try catch for each sample
  samparse <- function(i)
  {
    out <- try(rowparse(i))
    
    if (inherits(out, "try-error")) {
      return(paste(geterrmessage(), "in", i, "sample"))
    }
    return(out)
  }
  
  ##Main program  qdx qd
  ##rewrite for new data per row
  transloctable <- data.frame()
  
  for (i in 1:nrow(Con_data))
  {
    ##Cyto_sample=gsub("ish.*$","",as.character(Con_data[i,2]))
    Cyto_sample <-
      unlist(strsplit(as.character(Con_data[i, 2]), split = ",|\\["))
    ##if X or Y in first col, move it
    if(grepl("X|Y",Cyto_sample[1]))
    {
      Cyto_sample[1]<-gsub(substr(Cyto_sample[1],regexec("X|Y",Cyto_sample[1])[[1]][1],nchar(Cyto_sample[1])),paste(",", substr(Cyto_sample[1],regexec("X|Y",Cyto_sample[1])[[1]][1],nchar(Cyto_sample[1])),sep=""),Cyto_sample[1])
      Cyto_sample <-
        unlist(strsplit(as.character(Cyto_sample), split = ",|\\["))
    }
    ##take away extra accidental commas 
    if(length(which(str_length(Cyto_sample)==0))>0)
    {
      Cyto_sample<-Cyto_sample[-1*which(str_length(Cyto_sample)==0)]
    }
    ##if idem or sl, cancel out any - details 
    if(any(grepl("ids",Cyto_sample)))
    {
      ##index of anything thats - in a clonal evolution step
      mut_gone_index<-grep("-[[:alpha:]]+\\(&!?",Cyto_sample)
      if(length(mut_gone_index)>0)
      {
        mut_list<-sapply(Cyto_sample[mut_gone_index],function(x){gsub("\\)","\\\\)",gsub("\\(","\\\\(",substr(x,2,nchar(x))))})
        index_cancel<-sapply(mut_list,function(x){grep(x,Cyto_sample)[1]})
        Cyto_sample<-Cyto_sample[-1*c(mut_gone_index,index_cancel)]
      }
    }
    
    ##check to make sure the input is roughly a karyotype before processing
    if(grepl("[[:digit:]]+((~|-)[[:digit:]]+)*(<[[:digit:]]+n>)*,",paste(Cyto_sample,collapse=',',sep='')))
    {
      tottable<-samparse(i)
      if(is.character(tottable))
      {
        Dump_table <- rbind(Dump_table, c(Con_data[i, ], tottable))
        transloctable<-data.frame()
      }else{
        
        sorted_sample_table<-tottable[[1]]
        Dump_table<-tottable[[2]]
        ##master table with names
        temp2 <-
          cbind(rep(as.character(Con_data[i, 1]), nrow(sorted_sample_table)), sorted_sample_table)
        
        ##sort temp2 to remove gains and losses that are cooresponding
        
        
        colnames(temp2) <- colnames(Final_table)
        Final_table <- rbind(Final_table, temp2)
        ##correct data formates from factors into chr int int chr
        if(nrow(Final_table)==1)
        {
          Final_table<-as.data.frame(Final_table)
          Final_table[,1]<-as.character(Final_table[,1])   
          Final_table[,2]<-as.character(Final_table[,2])   
          Final_table[,3]<-as.integer(as.numeric(as.character(Final_table[,3])))
          Final_table[,4]<-as.integer(as.numeric(as.character(Final_table[,4])))
          Final_table[,5]<-as.character(Final_table[,5])
          
        }
        
        ##make list to store translocations and insertions
        ##if clone line, then, carry over transloctable to next one
        if(i<nrow(Con_data))
        {
          curname<-strsplit(Con_data[i,1],"_")[[1]]
          nextname<-strsplit(Con_data[i+1,1],"_")[[1]]
          if(curname[-length(curname)]==nextname[-length(nextname)] & as.numeric(curname[length(curname)])+1 ==as.numeric(nextname[length(nextname)]))
          {
            ##keep previous transloctable if sample names the same and cell line is in correct order
            transloctable<-tottable[[3]]
          }else
          {
            transloctable <- data.frame()
          }
        }
      }
    }else
    {
      ##Dump_table<-rbind(Dump_table,c(Con_data[i, ], "karyotype is incorrect"))
    }
  }
  
  
  
  ##take away other stuff
  ##if(is.vector(Final_table))
  ##{
  ##  Final_table<-Final_table[grep("Loss|Gain",Final_table[5]),]
  
  ##}
  ##else
  ##{
  ##  Final_table<-Final_table[grep("Loss|Gain",Final_table[,5]),]
  ##}
  ##Final_table<-na.omit(Final_table)
  ##take away NAS grep(NA,Final_table)
  ##i([[:alnum:]+])q
  
  ##remove duplicates and redundancys (same region loss/gain within group do this first, then duplicates)
  
  ##counting function for % mutated
  out <- Con_data
  
  col3 <- rep("unknown", nrow(out))
  
  IDs <- as.character(out[, 1])
  
  IDs <- lapply(strsplit(IDs, split = "_"), function(x) {
    x[-length(x)]
  })
  
  IDs <- unlist(lapply(IDs, function(x) {
    paste(x, collapse = "_")
  }))
  IDs <- unique(IDs)
  
  for (i in 1:length(IDs)) {
    ##get all ids with same base name
    ##gsub away all regular expressions 
    w <- grep(gsub("\\+","\\\\+",gsub("\\*","\\\\*",gsub("\\!","\\\\!",paste("^", IDs[i], "_", sep = "")))), out[, 1])
    ##if (length(w) > 1) {
    nxt <- unlist(lapply(strsplit(as.character(out[w, 2]), split = "\\["),
                         function(x) {
                           x[2]
                         }))
    
    ##nxt <- gsub("\\](,ids)*", "", nxt)
    
    ##deal with cp more gracefully
    is_cp<-grepl("cp",nxt)
    nxt <- as.numeric(gsub("\\](,ids)*|cp", "", nxt))
    ##if its cp just skip it--get more sophisticated processing later, that one unknonwn other stuff can be calculated
    ##maybe modifify later to have total, something sometihng etc 
    ##can make more complicated, so that cp stuff is unknown but added up to second sample, it creates a percentage for that
    if (any(is.na(nxt)))
    {
      col3[w] <- "unknown"
    }
    else{
      ##nxt <- nxt / sum(nxt)
      if(is_cp)
      {
        col3[w] <- paste("1-",nxt," of ",sum(nxt),sep='')
      }else
      {
        col3[w] <- paste(nxt,"of",sum(nxt))
      }
    }
    ##}
  }
  ##something is gonig wrong when its not just a normal cell and theres no more (says unknown instead of 1)
  ##turn ids back into over all list
  IDs <- as.character(Con_data[, 1])
  
  ##compare to final table
  if (is.vector(Final_table))
  {
    colvect <- vector(length = 1)
    
  } else{
    colvect <- vector(length = length(Final_table[, 1]))
  }
  
  for (i in 1:length(IDs))
  {
    if (is.vector(Final_table))
    {
      colvect[grep(IDs[i], Final_table[1])] <- col3[i]
      
    } else{
      colvect[grep(IDs[i], Final_table[, 1])] <- col3[i]
    }
  }
  
  ##
  
  Final_Final_table <- cbind(Final_table, colvect)
  colnames(Final_Final_table) <-
    c("Sample ID", "Chr", "Start", "End", "Type", "Cells Present")
  ##for ease of testing
  if (is.vector(Final_Final_table))
  {
    Final_Final_table <-
      Final_Final_table[grep("Loss|Gain", Final_Final_table[5]), ]
    
  } else
  {
    Final_Final_table <-
      Final_Final_table[grep("Loss|Gain", Final_Final_table[, 5]), ]
  }
  Final_Final_table <- na.omit(Final_Final_table)
  
  if (is.vector(Final_Final_table))
  {
    Final_Final_table <- t(Final_Final_table)
    colnames(Final_Final_table) <-
      c("Sample ID", "Chr", "Start", "End", "Type", "Percent Present")
  }
  
  result<-list(Final_Final_table,Dump_table)
  names(result)<-list("Result","Error_log")
return(result)
}