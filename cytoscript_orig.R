##setwd("G:/My Drive/BRB work/cyto_project/Cytogenetic software/R/")
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
  

if(require("DescTools")){
  print("DescTools is loaded correctly")
} else {
  print("trying to install DescTools")
  install.packages("DescTools")
  if(require("DescTools")){
    print("DescTools installed and loaded")
  } else {
    stop("could not install DescTools")
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


if(require("dplyr")){
  print("dplyr is loaded correctly")
} else {
  print("trying to install dplyr")
  install.packages("dplyr")
  if(require("dplyr")){
    print("dplyr installed and loaded")
  } else {
    stop("could not install dplyr")
  }  
}


if(require("hash")){
  print("hash is loaded correctly")
} else {
  print("trying to install hash")
  install.packages("hash", repos='https://mirrors.nics.utk.edu/cran')
  if(require("hash")){
    print("hash installed and loaded")
  } else {
    stop("could not install hash")
  }  
}
##takes first option for or 
## assumes order of abberations is from top to bottom 
##assumes all output goes top to bottom (eg no q20p23)
#' CytoConverter Function
#' 
#' This function accepts string or matrix input. Strings must be karyotypes and tables must have sample name in the first column and karyotype in the second column. The function outputs a table  
#' @param in_data
#' @param constitutional
#' @param build
#' @param guess
#' @param guess_q
#' @param guess_by_first_val
#' @param forMtn
#' @param orOption
#' @keyword
#' @export
#' @examples 
#' CytoConverter()
CytoConverter<-function(in_data,build="GRCh38",constitutional=T,guess=F,guess_q=F,guess_by_first_val=F,forMtn=T,orOption=T,sexstimate=F,complexSexEstimate=F,allow_Shorthand=F){

  
##if for dr. mitelman, overide parameters
  #if(forMtn==T)
 # {
#    constitutional=F
 #   quess_q=F
 #   guess=T
    ##dont make assumptions about sex or constitutional sex chromosomes if there is a XX? or X? or ? in sex chromosome field
#    sexstimate=F
    ##more complex handelintg pf sex choromsomes in clones and triploidy up 
 #   complexSexEstimate=T
 #   allow_Shorthand=F
#  }
  
##reference file to load
##specific loc
if(build =="GRCh38")
{
 
  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("Builds/cytoBand_GRCh38.txt", header = FALSE)
    ), as.character)
  
}else if(build =="hg19"){
 
  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("Builds/cytoBand_hg19.txt", header = FALSE)
    ), as.character)
  
}else if(build =="hg18"){
 
  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("Builds/cytoBand_hg18.txt", header = FALSE)
    ), as.character)
  
}else if(build=="hg17"){

  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("Builds/cytoBand_hg17.txt", header = FALSE)
    ), as.character)
  
}else if(is.null(build))
{
  ##default is grch38

  Cyto_ref_table <-
    sapply(as.data.frame(
      read.delim("Builds/cytoBand_GRCh38.txt", header = FALSE)
    ), as.character)
  
}else{
  return("Error : build incorrectly specified")
}
  
  
  ##reference object for end of the string 
  ref_table <-as.data.frame(Cyto_ref_table[sapply(unique(Cyto_ref_table[,1]),function(x){grep(x,Cyto_ref_table[,1])[length(grep(paste(x,"$",sep=""),Cyto_ref_table[,1]))]}),][,c(1,3)])
  ref_table<-apply(ref_table,2,as.character)
  
  
  
  
  ##table with desired output
  Final_table <- matrix(ncol = 5, nrow = 0)
  
  colnames(Final_table) <- c("Sample ID", "Chr", "Start", "End", "Type")
  
  ##convert any single string into a table
  if(is.vector(in_data))
  {
    in_data <- t(matrix(c("sample", in_data)))
  }
  
  ##dump table of stuff containing unprocessed reads
  ##Write this later
  ##get all fish recorded
  Dump_table <- matrix(ncol = 3, nrow = 0)
  ##double check that this does not delete later data potentially
  fish_table <- in_data[grep("ish.*$", in_data[, 2]), ]
  if (is.vector(fish_table) & length(fish_table)>0)
  {
    Dump_table <- rbind(Dump_table, c(fish_table, "Warning in fish reading"))
    ##now take out ish readings
    in_data[, 2] <- gsub("ish.*$", "", as.character(in_data[, 2]))
    
  }else {
    Dump_table <- rbind(Dump_table, cbind(fish_table, "Warning in fish reading"))
    ##now take out ish readings
    
    in_data[, 2] <- gsub("ish.*$", "", as.character(in_data[, 2]))
    
  }
  

  ##set everything to lowercase, then set x and y to uppercase
  in_data[, 2] <- tolower(in_data[, 2])
  in_data[, 2] <- chartr("x", "X", in_data[, 2])
  in_data[, 2] <- chartr("y", "Y", in_data[, 2])
  
  Con_data = matrix(nrow = 0, ncol = 2)
  ##spliting cell lines
  for (i in 1:nrow(in_data))
  {
    ##if activated, allow for shorthand outside of clones and outside of translocations eg. del(1), or if its a clone
    if(!is.na(in_data[i, 2]) && nchar(as.character(in_data[i, 2])) > 0 && (allow_Shorthand == T || grepl("/",in_data[i,2]))){
      ##string split this up
      data_split<-unlist(strsplit(in_data[i, 2], ","))
      data_split<-unlist(strsplit(gsub("/","/-;-;opj",data_split),"/"))
      mut_index<-grep("\\(.*\\)(?!\\()",data_split,perl=T)
      if(length(mut_index)>0)
      {
        mut_list<-sapply(data_split[mut_index],function(x){gsub("\\[|\\]","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",gsub("\\?","\\\\?",gsub("\\+","",substr(x,2,nchar(x))  )))))})
        index_match<-sapply(mut_list,function(x){grep(x,data_split)[1]})
        index_match<-index_match[intersect(which(!is.na(index_match)),which(index_match>0))]
        if(is.numeric(mut_index) && is.numeric(index_match)&& length(index_match)>0)
        {
          additional_to_add <- substring(names(index_match),0,regexpr("\\+",names(index_match)) )
          data_split[mut_index]<-paste(additional_to_add ,data_split[index_match],sep="")
        }
      }
      
      in_data[i,2]<-gsub(",-;-;opj","/",paste(data_split,collapse=","))
      
    }
    
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
        if(!is.vector(temp_data) && nrow(temp_data)>1)
        {
          for (j in 2:nrow(temp_data))
          {
                ##make sure temp data has more than 1 row
            
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
        }else{
          ##put it in the error table 
          Dump_table<-rbind(Dump_table, c(temp_data, "no proceeding cell line to refer to"))
          
        }
    }
    
    
    
  }
  
  rownames(Con_data)<-1:nrow(Con_data)
  
  
  Con_data[, 2] <- gsub(" ", "", Con_data[, 2])
  ##any rejoins will be ;; now
  ##actually lets depricate this
  ##Con_data[, 2] <- gsub(":", ";", Con_data[, 2])
  Con_data[, 2] <- gsub("crYp", "", Con_data[, 2])
  
  ## taking all reads with ? and ~ as well into dump table
  ##adding marker and add material as well
  unsure_table <- Con_data[grep("\\?|\\~|inc|mar|add", Con_data[, 2]), ]
  Dump_table <- if(is.vector(unsure_table)){
    rbind(Dump_table, c(unsure_table, "Warning in ?,~,marker, unknown additional material or incomplete karyotype detected"))
  }else{
    rbind(Dump_table, cbind(unsure_table, "Warning in ?,~, marker, unknown additional material or incomplete karyotype detected"))
    
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
    
    ##quit if karyotype returns false
    earlyReturn=F 
    
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
    
    # JP: if (any(grepl("::", temp[[lengthcount * 2]][i]) |
    # JP:         grepl("~>", temp[[lengthcount * 2]][i]) | grepl("->", temp[[lengthcount * 2]][i])  )) {
    if (any(grepl("::", temp[[lengthcount * 2]][o]) |
            grepl("~>", temp[[lengthcount * 2]][o]) | grepl("->", temp[[lengthcount * 2]][o])  )) {
      #parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##find p or q, take stuff before take stuff after, before is chromosomes after is positions, this will return 2 objects, must take into account
      ##parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
      ##ask tom about this one
      ##only splits first one
      longform_table <-
        strsplit(strsplit(temp[[lengthcount * 2]][o], "::")[[1]], "(~>)|(->)")
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
      
      
    
      
      
      ###########################################################################################################
      #####################for mitelman data only###############################################################
      #####check for more than one band per chromosome for translocations######################################
      ##############################################################################################################
      
      if(forMtn==T & grepl("t\\(",derMods[lengthcount] ) & length(positions) > 1  )
      {
        earlyReturn=T
        
      }else{
          
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
                                              -1, 4], chr_table[1, 4], sep = ''))
              ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])-1 
              ##                              , 4][1], chr_table[1, 4], sep = ''))
              ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4]) 
              ##                          , 4][1], chr_table[1, 4], sep = ''))
            }
          }else{
            
              ## restrict if positions are at the ends of the chromosome
              pos<-grep(paste(positions,sep="",collapse="|"), chr_table[, 4])
              if(pos[length(pos)] +1 >nrow(chr_table) && pos[1]==1)
              {
                currentvec<-c(
                  currentvec,
                  paste(chr_table[nrow(chr_table), 4], sep = ''),
                  paste(chr_table[1, 4], sep = '')
                )
                
              }else if(pos[length(pos)] +1 >nrow(chr_table)){
                ##if cyto is at the top
                currentvec<-   c(
                  currentvec,
                  paste(chr_table[nrow(chr_table), 4], sep = ''),
                  paste(chr_table[pos[1] -1, 4], chr_table[1, 4], sep = '')
                )
                
                
              }else if(pos[1]==1){
                ##if cyto is on the bottom
                
                currentvec<-    c(
                  currentvec,
                  paste(chr_table[pos[length(pos)] +
                                    1, 4], chr_table[nrow(chr_table), 4], sep = ''),
                  paste( chr_table[1, 4], sep = '')
                )
                
                
              }else{
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
    if(isiso)
    {
      currentvec<-sapply(currentvec,function(x){gsub(paste(unused,"[[:digit:]]+(\\.[[:digit:]])*",sep=''),temp[[(grep("ider",derMods))+1]],x)})
    }
    
    return(list(currentvec,earlyReturn))
  }
  
  
  
  ##make sure the order of the bands is in the order we like (p to q highers ps to low, lower qs to high )
  positionsorter<-function(positions){
    ##if(length(positions)==2)
    ##{
      ##if(grepl("q",positions[1])&grepl("p",positions[2]))
      ##{
      ##  tempstorage<-positions[2]
      ##  positions[2]<-positions[1]
      ##  positions[1]<-tempstorage
      ##}
      ##if(grepl("q",positions[1])&grepl("q",positions[2]))
      ##{
      ##  bands<-sapply(positions,function(x){as.numeric(gsub("q","",x))})
      ##  if(bands[1]>bands[2])
      ##  {
      ##    tempstorage<-positions[2]
      ##    positions[2]<-positions[1]
      ##    positions[1]<-tempstorage
      ##  }
      ##}
      ##if(grepl("p",positions[1])&grepl("p",positions[2]))
      ##{
        ##bands<-sapply(positions,function(x){as.numeric(gsub("p","",x))})
        ##if(bands[1]<bands[2])
        ##{
        ##  tempstorage<-positions[2]
        ##  positions[2]<-positions[1]
        ##  positions[1]<-tempstorage
        ##}
      ##}
    ##}
    
    ###########################333
    #######new
    ##########################
    ################################
    if(length(unlist(strsplit(positions,"-|~")))>1)
    {
      positions<-unlist(lapply(strsplit(positions,"-|~"),function(x){x[1]}))
    }
    
    return(positions)
  }
  
  
  
  ##function for separating normal data
  ##take into acc same chrom insestion
  ##add cen into here (pter qter analouge)
  
  parser <- function(coln, xmod, ymod, transloctable,addtot,Cyto)
  {
    Cyto_sample<-Cyto
    
    ##for or statements, take first statement
    Cyto_sample[coln]<-gsub("or.*$","",Cyto_sample[coln])
    
    
    ##if we are guessing ? marks
    if(guess_q == T & any(grepl("\\?",Cyto_sample[coln])))
    {
      Cyto_sample[coln] <- gsub("\\?","",Cyto_sample[coln])
    }
    
    ##if we are counting constitutional, substitute c
    if(constitutional==T)
    {
      Cyto_sample[coln]<-gsub("c$","",Cyto_sample[coln])
    }
    

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
    derModsBackup<-strsplit(derMods, "\\)")[[1]]
    
    ##skip second +t if first one is triggered
    plusT<-F
    
    ##regins
    regtranschrom=NULL
    
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
    } else if (grepl("der|rec|dic", Cyto_sample[coln]) &
               !grepl("ider|idic", Cyto_sample[coln]))
    {
      if(grepl("der|rec",Cyto_sample[coln]))
      {
        addBool <- paste(addBool, derMods[1], sep = "")
      }
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
      }else if((grepl("der|rec",Cyto_sample[coln]) && !is.na(derMods) && length(derMods)>1 && grepl("^r\\(",derMods[2]) )){
        derMods <-tail(derMods, length(derMods) - 1)
        temp <- tail(temp, length(temp) - 1)
         
      }else if(grepl("der\\([[:digit:]]+;[[:digit:]]+\\)t\\(|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\(",Cyto_sample)[coln])
      {
        
        ##handle der(11;13)t( and dic(11;13)t( esq cases here
        
        ##remember to search t(11;13 if its referring to previous call of translocation
        ##figure out how to do this 
        
        ##this will only work with two translocations and the translocation must immediantly follow 
        ##the der or dic
        
        if(grepl("der\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)\\(|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)\\(",Cyto_sample[coln])){
          ##translocation is fully described 
          if(grepl("der\\(",Cyto_sample[coln])){
            ##replace der with dic
            derMods[1]<-gsub("der","dic",derMods[1])  
            ##delete der so its treated like a pure dic
            addBool[1]<-gsub("der","dic",addBool[1])
          }
          
          ##delete translocation from both temp and derMods
          derMods<-derMods[-2]
          temp<-temp[-2]
        }else if(grepl("der\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)|dic\\([[:digit:]]+;[[:digit:]]+\\)t\\([[:digit:]]+;[[:digit:]]+\\)",Cyto_sample[coln])){
          ##if the translocation needs lookup 
          
          ##see if translocation is found
          ##else this thing crashes
          if(length(transloctable) > 0 && grepl(derMods[2],names(transloctable),perl=F,fixed=T))
          {
            ##last chunk of translocation indicating position
            addOnTrans<-names(transloctable)[grep(derMods[2],names(transloctable),perl=F,fixed=T)]
            addOnTrans<- strsplit(addOnTrans,"\\)")[[1]][2]
            addOnTransTemp<-strsplit(gsub("\\(|\\)","",addOnTrans),";")
            ##if the string is longer than just the translocation
            if(length(derMods)>2){
              derMods<-c(derMods[1:2],addOnTrans,derMods[3:length(derMods)])
              temp<-c(temp[1:2],addOnTransTemp,temp[3:length(temp)])
            }else{
              ##just the der and the translocation
              derMods<-c(derMods[1:2],addOnTrans)
              temp<-c(temp[1:2],addOnTransTemp)
            }
            
            if(grepl("der\\(",Cyto_sample[coln])){
              ##replace der with dic
              derMods[1]<-gsub("der","dic",derMods[1])  
              ##delete der so its treated like a pure dic
              addBool[1]<-gsub("der","dic",addBool[1])
              
            }
            
            ##delete translocation from both temp and derMods
            derMods<-derMods[-2]
            temp<-temp[-2]
          }else{
            print("translocation undefined")
            return(NULL)
          }
        }
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
    ##handle differently if we are not counting constitutional
    
    multi = 1 ##check if there is a X3 etc value , if there is pick that up, store, add to addtot, make it process through twice later
    if(constitutional==F)
    {
      temp_cyto<-gsub("(c$)|(c\\?$)","",Cyto_sample[coln])
      if (grepl("\\)X|\\)x", temp_cyto))
      {
        multi = unlist(strsplit(temp_cyto, ")X|)x"))[2]
        if (grepl("-|~", multi))
        {
          multi <- unlist(strsplit(multi, "~|-"))[1]
          
        }
        ##multi=gsub("[^0-9]","",multi)
        
        multi = as.numeric(multi)
        addBool <- paste(addBool, "multi", multi, sep = "")
      }
      
    }else{
      if (grepl("\\)X|\\)x", Cyto_sample[coln]))
      {
        multi = unlist(strsplit(Cyto_sample[coln], ")X|)x"))[2]
        if (grepl("-|~", multi))
        {
          multi <- unlist(strsplit(multi, "~|-"))[1]
          
        }
        ##multi=gsub("[^0-9]","",multi)
        
        multi = as.numeric(multi)
        addBool <- paste(addBool, "multi", multi, sep = "")
      }
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
    

 

    ##make sure you count + properly for ? marks or constitutional
    if((any(grepl("\\?|\\~",Cyto_sample[coln])) & any(grepl("\\+",Cyto_sample[coln])) )  |(constitutional==F & multi > 1 & grepl("(c$)|(c\\?$)",Cyto_sample[coln]) ))
    {
      if(constitutional==F & multi > 2 & !grepl("\\+",Cyto_sample[coln]) & grepl("(c$)|(c\\?$)",Cyto_sample[coln] )){
      addtot<-addtot+(1*(multi-2))
        
      }else{
      addtot<-addtot+(1*multi)
      }
    }
    
    #########
    ##if guess is true, try to process ? marks
    ##think about how this can affect counting  + and \\?
    ######new######
    #############################
    ##not processing things like t(9;22)(p?;q10)
    

    
    if(length(temp)==1){
      
      if (grepl("(t|ins)\\(", derMods[1] ))
      {
        
        
        transchrom <-
          str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[1]))),"\\)\\(.+?\\)",sep=''))
        ##if this is not labled
        if(is.na(transchrom))
        {
          transchrom <-
            str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[1]))),"\\)",sep=''))
        }
        
        regtranschrom<-gsub("\\+","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom)))
      }
    }
    
    
    if ((length(test) > 1 | any(grepl("p|q",temp))| grepl("(9;22)|(22;9)",paste(unlist(temp),collapse=';',sep=";")) |(!is.null(regtranschrom) && any(grepl(regtranschrom, names(transloctable))) ) )  & (!any(grepl("\\?|\\~", temp))  ) & !(guess_q==F & grepl("\\?",Cyto_sample[coln])) )
    {
      ##goes by steps of 2, odd indexes indicate chromosomes, even indicate positions
      ##length_temp<-(if((length(temp) / 2)==0.5){1}else{length(temp)/2})
      lengthcount=1
      repeat
      {
        if(lengthcount > (if(!is.integer((length(temp) / 2))){ceiling(length(temp)/2)}else{length(temp)/2})) break
        if(lengthcount> 60) {print("while loop not terminating");break}  
        
        
        Allchr <-
          c(Allchr, as.vector(paste(temp[[(lengthcount * 2-1)]], "$", sep = "")))
        
        ##handle those t(__;__) and ins (__;__) here with no follow up
        if ((((lengthcount * 2)-1)) == length(temp) || if(length(temp) > (lengthcount*2-1)){
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
            if(any(grepl("p",temp[[(lengthcount*2-1)]])))
            {
              arm="p10"
            }
            if(any(grepl("q",temp[[(lengthcount*2-1)]])))
            {
              arm="q10"
            }
            temp<-c(if((lengthcount*2-1)!=1){temp[1:(lengthcount*2-1)]}else{temp[(lengthcount*2-1)]},arm,if(length(temp)>=2){temp[[(lengthcount*2):length(temp)]]})
            temp[[(lengthcount*2-1)]]<-gsub("p|q","",temp[[(lengthcount*2-1)]])
            
        }
         
          
          if(any(grepl("^r\\(",derMods[lengthcount*2-1])) & forMtn==F){
           
             rindex<-lengthcount*2-1
             
             
               rChr<-paste("chr",temp[[rindex]],"$",sep="")
               tempRingPosition<-NULL
               tempRing<-NULL
               
                for(c in 1:length(rChr))
                {
                   
                   tempRingTable<-matrix(nrow=0,ncol=5)
                   chr<-Cyto_ref_table[grep(rChr[c],Cyto_ref_table[,1]),] 
                   tempRingTable<-rbind(tempRingTable,rbind(chr[1,],chr[nrow(chr),]))
                   tempRingTable[,4]
                   
                  tempRingPosition<- c(tempRingPosition,paste(tempRingTable[,4],collapse=""))
                   
                }
               
                 tempRing<-list(unlist(temp[[rindex]]), tempRingPosition)
                 derModsRing<-c(paste("r(",paste(unlist(temp),sep="",collapse=";"),sep=""),paste("(",paste(tempRingPosition,collapse = ";"),sep=""))
                 tempstorage<-NULL
                 tempdermods<-NULL
                 
                 if((length(temp))>(lengthcount*2))
                 {
                   tempstorage<-temp[((lengthcount*2)):length(temp)]
                   tempdermods<-derMods[((lengthcount*2)):length(temp)]
                   
                   
                 }
                 
                 temp[[(lengthcount*2-1)]]<-tempRing[[1]]
                 
                 if(length(temp)==1)
                 {
                   temp<-list(unlist(temp),tempRing[[2]])
                 }else{
                   
                   temp[[lengthcount*2]]<-tempRing[[2]]
                   
                 }
                 temp<-temp[1:(lengthcount*2)]
                 
                 derMods[(lengthcount*2)-1]<-derModsRing[1]
                 
                 if(length(derMods)==1)
                 {
                   derMods<-c(derMods,derModsRing[2])
                 }else{
                   derMods[lengthcount*2]<-derModsRing[2]
                   
                 }
                 
                 derMods<-derMods[1:(lengthcount*2)]
                 
                 if(!is.null(tempstorage))
                 {
                   temp<-c(temp,tempstorage)
                   derMods<-c(derMods,tempdermods)
                   
                 }
               
              
             
          } 
          
          
          if (grepl("(t|ins)\\(", derMods[(lengthcount * 2-1)]) )
          {
            ##consider adding t(num;num) option
            
            transchrom <-
              str_extract(Cyto_sample[coln], gsub("\\?","\\\\?",paste(gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)])),"\\)\\(.+?\\)",sep='')))
            ##if this is not labled
            if(is.na(transchrom))
            {
              transchrom <-
                str_extract(Cyto_sample[coln], gsub("\\?","\\\\?",paste(gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)])),"\\)",sep='')))
            }
            
            regtranschrom<-gsub("\\+","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom)))
            temptrans <- data.frame()
            
            
            
            ##make stuff now
            if( !any(grepl(regtranschrom, names(transloctable))) & ((lengthcount*2)> length(temp)|| all(!grepl("p|q",temp[[(lengthcount*2)]]))) & all(grepl("9|22",transchrom)))
            {
              temptrans<-
                rbind(temptrans,cbind(
                  paste("der(", "9", ")", sep = ''),
                  paste(
                    "der(",
                    "9",
                    ";",
                    "22",
                    ")(",
                    "q34.1p24.3",
                    ";",
                    "q13.33q11",
                    ")",
                    sep = '',
                    collapse = ";"
                  )
                  
                )
                )
              
              temptrans<-
                rbind(temptrans,cbind(
                  paste("der(", "22", ")", sep = ''),
                  paste(
                    "der(",
                    "9",
                    ";",
                    "22",
                    ")(",
                    "q34.3q34.1",
                    ";",
                    "q11.2p13",
                    ")",
                    sep = '',
                    collapse = ";"
                  )
                  
                )
                )
            
              
              transloctable <- c(transloctable, list(temptrans))
              names(transloctable)[length(transloctable)] <-
                gsub("\\+","",transchrom)
              }
            
            
            if(grepl("^\\++t\\(",derMods[(lengthcount*2-1)]) & !grepl("der",Cyto_sample[coln]))
            {
              plusT=T
              
              transderplus<-transloctable[grep(regtranschrom,names(transloctable))]
              dertransextract<-transderplus[[1]][grep(paste("der\\(",temp[[(lengthcount*2)-1]],"\\)",sep="",collapse="|"),transderplus[[1]][,1]),2]
              
              transtemp <-
                unlist(lapply(as.list(dertransextract),function(x){strsplit(gsub("[\\(\\)]", "", regmatches(
                  x, gregexpr("\\(.*?\\)",  x)
                )[[1]]), ";")}),recursive=F)
              
              transdermod<-unlist(lapply(dertransextract,function(x){y<-paste("+",x,sep="");strsplit(y, "\\)")}))
              ##make temp storage for rest of temp
              tempstorage<-NULL
              tempdermods<-NULL
              if((length(temp))>(lengthcount*2))
              {
                tempstorage<-temp[[((lengthcount*2)+1):length(temp)]]
                tempdermods<-derMods[((lengthcount*2)+1):length(temp)]
                
                
              }
              
              temp[[(lengthcount*2-1)]]<-transtemp[[1]]
              
              if(length(temp)==1)
              {
                temp<-list(unlist(temp),transtemp[[2]])
              }else{
                
                temp[[lengthcount*2]]<-transtemp[[2]]
                
              }
              temp<-c(temp[1:(lengthcount*2)],transtemp[3:length(transtemp)])
              
              derMods[(lengthcount*2)-1]<-transdermod[1]
              if(length(derMods)==1)
              {
                derMods<-c(derMods,transdermod[2])
              }else{
                derMods[lengthcount*2]<-transdermod[2]
                
              }
              derMods<-c(derMods[1:(lengthcount*2)],transdermod[3:length(transdermod)])
              
              
              if(!is.null(tempstorage))
              {
                temp<-c(temp,tempstorage)
                derMods<-c(derMods,tempdermods)
                
              }
              
            }else{
            
              transchrom <-
                str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]))),"\\)",sep=''))
              selectedTransTable<-as.matrix(transloctable[grep((gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom))),names(transloctable))][[1]])
              transDer<-selectedTransTable[grep(paste("der\\(",paste(sapply(Mainchr,function(x){substr(x,0,nchar(x)-1)}),sep='',collapse=';'),"\\)",sep=''),selectedTransTable[,1]),][2]
              transtemp <-
                strsplit(gsub("[\\(\\)]", "", regmatches(
                  transDer, gregexpr("\\(.*?\\)",  transDer)
                )[[1]]), ";")
              #####
              ##this isnt working 
              ###
              temp<-c(if((lengthcount*2-1)!=1){temp[1:((lengthcount*2-1))]}else{transtemp[1]},transtemp[2],if((length(temp)>=lengthcount*2) && (!grepl("p|q|->|:",derMods[length(temp)]) || (IsOdd(length(temp))))){temp[lengthcount*2:length(temp)]})
            }
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
        if (grepl("t\\(", derMods[(lengthcount * 2-1)]) ) {
          ##if translocation , do this , shouldnt be placed , should be placed one loop above
          ## if not in the table, add to table
          ##names are just t(1;12)
          
          ######################################################################
          ###################################################333
          ##############################################################################
          #######################################################################3
          ##if its t(9;22) or t(22;9) , add to table if not already found
          #########################################################################
          #########################################################################
          #############################################################################
          #######################################################################
          
          transchrom <-
            str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]))),"\\)\\(.+?\\)",sep=''))
          ##if this is not labled
          if(is.na(transchrom))
          {
            transchrom <-
              str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]))),"\\)",sep=''))
          }
          regtranschrom<-gsub("\\+","",gsub("\\)","\\\\)",gsub("\\(","\\\\(",transchrom)))
          if (!any(grepl(regtranschrom, names(transloctable))))
          {
            ###
            ########
            ##entire translocation only if nessesary
            ##str_extract(mem, "t\\([[:digit:]]+(;[[:digit:]])*?\\)\\(.+?\\)")
            temptrans <- data.frame()
            
            
            
            ##make stuff now
            if(  !any(grepl(regtranschrom, names(transloctable))) &((lengthcount*2)> length(temp)|| all(!grepl("p|q",temp[[(lengthcount*2)]]))) & all(grepl("9|22",transchrom)))
            {
              temptrans<-
                rbind(temptrans,cbind(
                  paste("der(", "9", ")", sep = ''),
                  paste(
                    "der(",
                    "9",
                    ";",
                    "22",
                    ")(",
                    "q34.1p24.3",
                    ";",
                    "q13.33q11",
                    ")",
                    sep = '',
                    collapse = ";"
                  )
                  
                )
                )
              
              temptrans<-
                rbind(temptrans,cbind(
                  paste("der(", "22", ")", sep = ''),
                  paste(
                    "der(",
                    "9",
                    ";",
                    "22",
                    ")(",
                    "q34.3q34.1",
                    ";",
                    "q11.2p13",
                    ")",
                    sep = '',
                    collapse = ";"
                  )
                  
                )
                )
              
            }else {
              for (o in 1:length(temp[[(lengthcount * 2-1)]]))
              {
                tempcurvec<- getCytoBands(lengthcount, o, temp,coln,derMods)
                currentvec <- tempcurvec[[1]]
                earlyReturn<-tempcurvec[[2]]
                
                if(earlyReturn==T){
                  return(NA)
                }else{
                  ##vector of o cytobands ##includes everything but cytoband on der o indicated, double check how your doing this (going one band before, is this right?)
                  
                  
                  ##above is all to find vectors of o that are excluded
                  if (o == 1)
                  {
                    temptrans <-
                      rbind(temptrans, cbind(
                        paste("der(", temp[[(lengthcount * 2-1)]][o], ")", sep = ''),
                        paste(
                          "der(",
                          temp[[(lengthcount * 2-1)]][o],
                          ";",
                          temp[[(lengthcount * 2-1)]][length(temp[[(lengthcount * 2-1)]])],
                          if (length(currentvec) > 1) {
                            ";"
                          },
                          rep(temp[[(lengthcount * 2-1)]][o], length(currentvec) - 1),
                          ")(",
                          currentvec[1],
                          ";",
                          temp[[lengthcount * 2]][length(temp[[(lengthcount * 2-1)]])],
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
                          temp[[(lengthcount * 2-1)]][o - 1],
                          if (length(currentvec) > 1) {
                            ";"
                          },
                          rep(temp[[(lengthcount * 2-1)]][o], length(currentvec) - 1),
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
                      ))
                    
                  }
                }
              }
            }
            transloctable <- c(transloctable, list(temptrans))
            names(transloctable)[length(transloctable)] <-
              gsub("\\+","",transchrom)
          }
          
          
          ##for extending der(X)t(X;y) type 
          if(!grepl("^\\++t\\(",derMods[(lengthcount*2-1)]) & grepl("der",Cyto_sample[coln]))
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
          
          ########
          ####new #####
          ## fix +t() here
          ################
          if(grepl("^\\++t\\(",derMods[(lengthcount*2-1)]) & !grepl("der",Cyto_sample[coln]) & plusT==F)
         {
            print(lengthcount)
            transderplus<-transloctable[grep(regtranschrom,names(transloctable))]
            dertransextract<-transderplus[[1]][grep(paste("der\\(",temp[[(lengthcount*2)-1]],"\\)",sep="",collapse="|"),transderplus[[1]][,1]),2]
            
           transtemp <-
              unlist(lapply(as.list(dertransextract),function(x){strsplit(gsub("[\\(\\)]", "", regmatches(
                x, gregexpr("\\(.*?\\)",  x)
              )[[1]]), ";")}),recursive=F)
            
            transdermod<-unlist(lapply(dertransextract,function(x){y<-paste("+",x,sep="");strsplit(y, "\\)")}))
            ##make temp storage for rest of temp
            tempstorage<-NULL
            tempdermods<-NULL
            if((length(temp))>(lengthcount*2))
            {
              tempstorage<-temp[[((lengthcount*2)+1):length(temp)]]
              tempdermods<-derMods[((lengthcount*2)+1):length(temp)]
            
              
            }
            
            temp[[(lengthcount*2-1)]]<-transtemp[[1]]
            
            if(length(temp)==1)
            {
              temp<-list(unlist(temp),transtemp[[2]])
            }else{
              
              temp[[lengthcount*2]]<-transtemp[[2]]

            }
            temp<-c(temp[1:(lengthcount*2)],transtemp[3:length(transtemp)])
            
            derMods[(lengthcount*2)-1]<-transdermod[1]
            if(length(derMods)==1)
            {
              derMods<-c(derMods,transdermod[2])
            }else{
              derMods[lengthcount*2]<-transdermod[2]
              
            }
            derMods<-c(derMods[1:(lengthcount*2)],transdermod[3:length(transdermod)])
            
            
            if(!is.null(tempstorage))
            {
              temp<-c(temp,tempstorage)
              derMods<-c(derMods,tempdermods)
              
            }
            
          }
          
          
        }
          
        if (grepl("ins\\(", derMods[(lengthcount * 2-1)]) ) {
          currentvec<-vector()
          inschrom <-str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]))),"\\)\\(.+?\\)",sep=''))
          
          if(is.na(inschrom))
          {
            inschrom <-
              str_extract(Cyto_sample[coln], paste(gsub("\\?","\\\\?",gsub("\\+","\\\\+",gsub("\\(","\\\\(",derMods[(lengthcount*2-1)]))),"\\)",sep=''))
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
          if(length(temp[[lengthcount*2-1]])==1 && any(length(temp[[(lengthcount * 2-1)]]) < str_count(temp[[lengthcount*2]],"p|q")))
          {
            temp[[(lengthcount * 2-1)]]<- rep(temp[[(lengthcount * 2-1)]],2)
            
            ##location to split string (2nd q or p)
            splitloc<-str_locate_all(temp[[lengthcount*2]],"p|q")[[1]][,1][2]
            
            temp[[lengthcount*2]]<-c(str_sub(temp[[lengthcount*2]],1,splitloc-1),str_sub(temp[[lengthcount*2]],splitloc,str_length(temp[[lengthcount*2]])))
            within=T
          }
          
          if (!any(grepl(reginschrom, names(transloctable))))
          {
            ##entire translocation only if nessesary
            tempins <- data.frame()
            
             
            currentvec <-getCytoBands(lengthcount, 1,temp,coln,derMods)[[1]]
            
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
            
            ##change this to be a del chromosome
            currentvec <- getCytoBands(lengthcount, 2,temp,coln,derMods)[[1]]
            
            ##tempins <-
            ##rbind(tempins, cbind(
            ##paste("der(", temp[[(lengthcount * 2-1)]][2], ")", sep = ''),
            ##paste(
            ##"del(",
            ##temp[[(lengthcount * 2-1)]][2],
            ##")(",
            ##temp[[(lengthcount*2)]][2],
            ##")",
            ##sep = '',
            ##collapse = ";"
            ##),
            ##""
            ##))
            
            tempins<-rbind(tempins, cbind(
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
          
          if(length(test) > 1 | plusT){
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
            if (grepl("~|(- &!(->))", temp[[lengthcount * 2]][i]) )
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
                    grepl("~>", temp[[lengthcount * 2]][i]) |   grepl("->", temp[[lengthcount * 2]][i]) )) {
              #parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
              ##ask tom about this one
              ##find p or q, take stuff before take stuff after, before is chromosomes after is positions, this will return 2 objects, must take into account
              ##parse data according to ::, in front of p and q are chromosomes, if qter or pter, do stuff, afterward is position, make table of things included, then make list of stuff excluded
              ##ask tom about this one
              ##only splits first one
              longform_table <-
                strsplit(strsplit(temp[[lengthcount * 2]][i], "::")[[1]], "(~>)|(->)")
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
                if(length(unlist(strsplit(positions,"-|~")))>1)
                {
                  positions<-unlist(lapply(strsplit(positions,"-|~"),function(x){x[1]}))
                }
                ##put stuff in table
                positions_table <-
                  matrix(chr_table_2[grep(paste(positions, collapse = "|", sep = "|"),
                                          chr_table_2[, 4]), ], ncol = 5)
                
                if (is.vector(positions_table)| nrow(positions_table)==1 |ncol(positions_table)==1)
                {
                  if(ncol(positions_table)==1)
                  {
                    positions_table<-t(positions_table)
                  }
                  in_table <-
                    rbind(
                      in_table,
                      cbind(
                        chr_table_2[1, 1],
                        positions_table[2],
                        positions_table[3],
                        derMods[(lengthcount * 2-1)]
                      )
                    )
                  
                }else{
                  
                  
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
                
              }
              
              
              
              
              ##detect multis, combine with coord
              ##do multi (X2) processing 
              ##for over X2 times, its a gain 
              ##if (nrow(in_table) > 0 && any(grepl("multi", in_table[, 4])))
              ##{
              
              ##n <- multi - 1
              ##multimasterin_table<-in_table
              ## multimasterexin_table<-exin_table
              ##for (f in 1:n)
              ##{
              ##if(f>1)
              ##{
              ##multitemp<-in_table[grep("multi", in_table[, 4]), ]
              ##multitemp[,4]<-paste("+",multitemp[,4],sep='')
              ##multimasterin_table<- rbind(multimasterin_table, multitemp)
              ##if (any(nrow(exin_table) > 0 && grepl("multi", exin_table[, 4])))
              ##{
              ##multitemp<-exin_table[grep("multi", exin_table[, 4]), ]
              ##  multitemp[,4]<-paste("+",multitemp[,4],sep='')
              ##  multimasterexin_table <- rbind(multimasterexin_table, multitemp)
              ##  }
              
              ##}else{
              ##multimasterin_table<- rbind(in_table, in_table[grep("multi", in_table[, 4]), ])
              ##  if (any(nrow(exin_table) > 0 && grepl("multi", exin_table[, 4])))
              ##{
              ##multimasterexin_table <- rbind(multimasterexin_table, exin_table[grep("multi", exin_table[, 4]), ])
              
              ##}
              ##}
              ##}
              ##in_table<-multimasterin_table
              ##excoord<-multimasterexcoord
              
              
              ##}
              
              
              #add to coord
              
              ##combine piecewise with total count    
              coord <- rbind(coord, in_table)
              ##make note where the breaks are
              ##excoord<-rbind(excoord,) ##constant exclusion
              
            } else {
              in_table<-data.frame()
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
              if(length(unlist(strsplit(positions,"-|~")))>1)
              {
                positions<-unlist(lapply(strsplit(positions,"-|~"),function(x){x[1]}))
              }
              
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
                    in_table <-
                      rbind(in_table,
                            cbind(chr_table[1, 1], "0", positions_table[3], derMods[lengthcount * 2 -
                                                                                      1]))
                    
                  } else{
                    in_table <-
                      rbind(in_table,
                            cbind(chr_table[1, 1], "0", positions_table[nrow(positions_table), 3], derMods[lengthcount *
                                                                                                             2 - 1]))
                  }
                }
                if (any(grepl("q", positions)))
                {
                  if (is.vector(positions_table))
                  {
                    in_table <-
                      rbind(
                        in_table,
                        cbind(
                          chr_table[1, 1],
                          positions_table[2],
                          ref_table[grep(paste(chr_table[1, 1], "$", sep = ""), ref_table), ][2],
                          derMods[(lengthcount * 2-1)]
                        )
                      )
                    
                  } else{
                    in_table <-
                      rbind(
                        in_table,
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
                in_table <-
                  rbind(
                    in_table,
                    cbind(
                      chr_table[1, 1],
                      positions_table[1, 2],
                      positions_table[nrow(positions_table), 3],
                      derMods[(lengthcount * 2-1)]
                    )
                  )
              }
              
              ##do multi (X2) processing 
              ##for over X2 times, its a gain 
              ##if (nrow(in_table) > 0 && any(grepl("multi", in_table[, 4])))
              ##{
              
              ##n <- multi - 1
              ##multimasterin_table<-in_table
              ## multimasterexin_table<-exin_table
              ##for (f in 1:n)
              ##{
              ##if(f>1)
              ##{
              ##multitemp<-in_table[grep("multi", in_table[, 4]), ]
              ##multitemp[,4]<-paste("+",multitemp[,4],sep='')
              ##multimasterin_table<- rbind(multimasterin_table, multitemp)
              ##if (any(nrow(exin_table) > 0 && grepl("multi", exin_table[, 4])))
              ##{
              ##multitemp<-exin_table[grep("multi", exin_table[, 4]), ]
              ##  multitemp[,4]<-paste("+",multitemp[,4],sep='')
              ##  multimasterexin_table <- rbind(multimasterexin_table, multitemp)
              ##  }
              
              ##}else{
              ##multimasterin_table<- rbind(in_table, in_table[grep("multi", in_table[, 4]), ])
              ##  if (any(nrow(exin_table) > 0 && grepl("multi", exin_table[, 4])))
              ##{
              ##multimasterexin_table <- rbind(multimasterexin_table, exin_table[grep("multi", exin_table[, 4]), ])
              
              ##}
              ##}
              ##}
              ##  in_table<-multimasterin_table
              ##excoord<-multimasterexcoord
              
              
              
              ##}
              
              ##combine piecewise with total count    
              coord <- rbind(coord, in_table)
              ##make note where the breaks are
              ##excoord<-rbind(excoord,) ##constant exclusion
              
            }
            
            
            
          } }
        }
      
        
        lengthcount <- lengthcount+1  
       
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
      ##coord[,2:3]<-t(apply(coord[,2:3],1,function(x){sort(x)}))
      
      ##if there are two translocations, fix coordinate overlap so only overlap counts,
      
      ##if(str_count(Cyto_sample[coln],"t\\(")>1)
      ##{
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
      
      ##make sure coords are in numerical order 
      coord<-coord[order(coord[,1],coord[,2],coord[,3]),]
      
      
      
      ##all chromosomes used
      Chr<-unique(c(Mainchr,Allchr))
      
      
      #make sure translocations are consistant
      for(p in 1:length(Chr))
      {
        curChr <- Chr[p]
        
        tempChr = coord[grep(paste(
          "chr",
          curChr,
          "$" ,
          sep = "",
          collapse = "|"
        ),
        coord[, 1]), ]
        tempChr = apply(tempChr, 2, as.character)
        
        
        
        if(!is.vector(tempChr) && grepl("t\\(",tempChr[,4]) && nchar(addBool)>0)
        {
          ##sub out translocations--delete areas of no overlap, mainchr only, translocatios only
          ##transloc<-tempChr[grep(paste("der\\(",gsub("\\$","",curChr),"t\\(",sep=""),tempChr[,4]),]
          transloc<-tempChr[grep("t\\(",tempChr[,4]),]
          if((is.data.frame(transloc) | is.matrix(transloc))&&nrow(transloc) > 1)
          {
            if(DescTools::Overlap(as.numeric(transloc[1,2:3]), as.numeric(transloc[2,2:3])))
            {  
              
              ##if nothing , handle
              if(as.numeric(transloc[1,3])==as.numeric(transloc[2,2]))
              {
                temptransloc<-c(tempChr[1,1],as.numeric(transloc[1,2]),as.numeric(transloc[2,3]),tempChr[1,4])
                tempChr<-tempChr[-1*grep("t\\(",tempChr[,4]),]
                
                
                tempChr<-rbind(tempChr,temptransloc)
                coord<-rbind(temptransloc,coord[-1*intersect(grep("t\\(",coord[,4]),grep(paste("chr",curChr,sep=""),coord[,1])),])
                
              }else{
              temptransloc<-c(tempChr[1,1],mergeIntOverlap(as.numeric(transloc[1,2:3]),as.numeric(transloc[2,2:3])),tempChr[1,4])
              tempChr<-tempChr[-1*grep("t\\(",tempChr[,4]),]
              
              
              tempChr<-rbind(tempChr,temptransloc)
              coord<-rbind(temptransloc,coord[-1*intersect(grep("t\\(",coord[,4]),grep(paste("chr",curChr,sep=""),coord[,1])),])
              }
            }
          }
        }
        if(is.data.frame(tempChr)&&nrow(tempChr)==1)
        {
          tempChr<-as.vector(tempChr)
        }
      }
      
      
      ##make sure coords are in numerical order again
      coord<-coord[order(coord[,1],coord[,2],coord[,3]),]
      
      
      ##make excoord table
      
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
        
        if(length(tempChr)==0)
        {
          ###go back to here
          
          excoord = rbind(excoord,cbind(paste("chr",gsub("\\$","",curChr),sep=""),0,ref_table[grep(paste("chr",curChr, sep =""), ref_table), ][2],derModsBackup[1]))
          
        }else if (is.vector(tempChr))
        {
          if (!identical(gsub(" ","",as.character(tempChr[2])), "0"))
          {
            excoord = rbind(excoord, cbind(tempChr[1], 0, tempChr[2], tempChr[4]))
          }
          
          if (!identical(gsub(" ","",as.character(tempChr[3])),
                         gsub(" ","",as.character(ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table), ][2]))))
          {
              excoord = rbind(excoord,
                            cbind(tempChr[1], tempChr[3], ref_table[grep(paste(tempChr[1], "$", sep =
                                                                                ""), ref_table), ][2], tempChr[4]))
          }
        }
        else{
          
          ##make sure none are inside intervals, if it is, delete out of tempChr 
          for(z in 1:nrow(tempChr))
          {
            inbetween<-F
            for(x in 1:(nrow(tempChr))){

              if((as.numeric(tempChr[x,2])<as.numeric(tempChr[z,2]) & as.numeric(tempChr[x,3]>tempChr[z,3])))
              {
                inbetween<-T
              }
              
              
            }
            
            if(inbetween==T)
            {
              tempChr[z,]<-c(0,0,0,0)
            }
        }
           
          tempChr<-tempChr[which(tempChr[,1]!="0"),]
            
          
          ##if its a vector, do top part again
                if (is.vector(tempChr))
                {
                  if (!identical(gsub(" ","",as.character(tempChr[2])), "0"))
                  {
                    excoord = rbind(excoord, cbind(tempChr[1], 0, tempChr[2], tempChr[4]))
                  }
                  
                  if (!identical(gsub(" ","",as.character(tempChr[3])),
                                 gsub(" ","",as.character(ref_table[grep(paste(tempChr[1], "$", sep = ""), ref_table), ][2]))))
                  {
                    excoord = rbind(excoord,
                                    cbind(tempChr[1], tempChr[3], ref_table[grep(paste(tempChr[1], "$", sep =
                                                                                         ""), ref_table), ][2], tempChr[4]))
                  }
                }else{
                          
                    #else
                    for (z in 1:nrow(tempChr))
                    {
                      
                      
                      if (z == 1)
                      {
                        if( !identical(gsub(" ","",as.character(tempChr[1, 2])), "0"))
                        {
                          excoord = rbind(excoord,
                                          cbind(tempChr[z, 1], 0, tempChr[z, 2], tempChr[z, 4]))
                        }
                        
                        if (nrow(tempChr) > 1 && !identical(gsub(" ","",as.character(tempChr[z,3])),gsub(" ","",as.character(tempChr[z+1,2]))))
                        {
                          excoord = rbind(excoord,
                                          cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4]))
                        }
                        
                        
                        
                      }else if (z == nrow(tempChr) &
                                 !identical(gsub(" ","",as.character(tempChr[z, 3])),
                                            gsub(" ","",as.character(ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep =
                                                                                                       "")), ref_table), ][2])))) {
                        excoord = rbind(excoord,
                                        cbind(tempChr[z, 1], tempChr[z, 3], ref_table[grep(as.character(paste(tempChr[z, 1], "$", sep =
                                                                                                                "")), ref_table), ][2], tempChr[z, 4]))
                        
                      } else if (z < nrow(tempChr) & z > 1 && !identical(gsub(" ","",as.character(tempChr[z,3])),gsub(" ","",as.character(tempChr[z+1,2])))) {
                        
                        ##if this interval is not inside another 
                          excoord = rbind(excoord,
                                        cbind(tempChr[z, 1], tempChr[z, 3], tempChr[z + 1, 2], tempChr[z, 4]))
                        
                        
                      }
            
            
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
    if (length(excoord) > 0 && !is.vector(excoord)&& ncol(excoord)>1)
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
        if(nrow(excoord[-bads,])>0){
          excoord<-excoord[-bads,]
        }else
        {
          excoord<-data.frame()
        }
      }
      if(is.vector(excoord) && length(excoord)>0)
      {
        excoord[4]<-paste(addBool, excoord[ 4], sep = "")
      }else if(ncol(excoord)==1){
        excoord<-t(excoord)
        excoord[4]<-paste(addBool, excoord[ 4], sep = "")
        
      }else if(length(excoord)>0)
      {
        excoord[, 4] <- paste(addBool, excoord[, 4], sep = "")
      }
    }else if(ncol(excoord)==1)
    {
      excoord<-t(excoord)
      ecoord[2:3]<-as.numeric(as.character(excoord[2:3]))
      excoord[4]<-paste(addBool, excoord[ 4], sep = "")
    }else if(is.vector(excoord) && length(excoord) > 0){
      
      ecoord[2:3]<-as.numeric(as.character(excoord[2:3]))
      excoord[4]<-paste(addBool, excoord[ 4], sep = "")
    }
    
    
    
    ##make sure this does what i want it to
    ##excoord<-unique(excoord)
    ##coord<-unique(coord)
    
    ##}
    
    
    ##have to fix derivative chromosomes  with two translocations , if they overlap in a chromosome , its the overlapping increment that maters
    
    
    ##do multi (X2) processing right now
    ##for over X2 times, its a gain 
    if ((nrow(coord) > 0 && any(grepl("multi", coord[, 4])) ) && (constitutional==T| (constitutional ==F & !grepl("(c$)|(c\\?$)",coord[,4]) )))
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
  colparse <- function(j, xmod, ymod, transloctable,addtot,Cyto)
  {
    out <- try(parser(j, xmod, ymod, transloctable,addtot,Cyto_sample))
    
    if (inherits(out, "try-error")) {
      return(gsub("\n"," ",paste(geterrmessage(), "in", j, "field")))
    }
    return(out)
  }


  
  
  ### Now the whole deal
  
  ##for merginign loss gain pairs (returns index )
  ##consolidatesimple<-function(gTable,lTable,i){ 
  ##think i should do this via recursion some how
  ##return(which(lTable[,1] == gTable[i,1] & lTable[,2] == gTable[i,2] & lTable[,3] == gTable[i,3]))
  
  ##}
  
  
  ### Now the whole deal
  
  mergeTable_jp <- function(M,keep_extras=F) {
     ##store non gains and losses for intermediate steps
     M_temp<-M[grep("Gain|Loss",M[,4],invert=T),]
      
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
      
      for (h_key in keys(h_)) {
        if (start > as.numeric(h_key) & start < h_[[h_key]][['End']] ) {
          # add a new section starting at "start"
          h_[start] <- hash(
            End=h_[[h_key]][['End']],
            Gain=h_[[h_key]][['Gain']],
            Loss=h_[[h_key]][['Loss']],
            "+Loss"=h_[[h_key]][['+Loss']]
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
      
      for (h_key in keys(h_)) {
        if (end > as.numeric(h_key) & end < h_[[h_key]][['End']] ) {
          # add a new section starting at "end"
          h_[end] <- hash(
            End=h_[[h_key]][['End']],
            Gain=h_[[h_key]][['Gain']],
            Loss=h_[[h_key]][['Loss']],
            "+Loss"=h_[[h_key]][['+Loss']]
          )
          # modify the old section to end at "end"
          h_[[h_key]][['End']] <- end
          
          break
        }
      }
      
      ### Now add new section by splitting it accordingly
      start_vals <- sort(as.numeric(keys(h_)))
      prev_end <- start
      for (start_val in start_vals) {
        
        ## Check if existing section is completely within new region
        ## Update existing section accordingly by adding to gain/loss
        ## EXAMPLE: Adjust gain/loss of ** sections
        ## h_     |-----|*****|       |**|------|
        ## new   start->|----------------|<-end
        
        if (start_val >= start & h_[[as.character(start_val)]][['End']] <= end) {
          # keep track of original section order
          h_[[as.character(start_val)]][[type]][[
            length(h_[[as.character(start_val)]][[type]])+1
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
            
            h_[prev_end] <- hash(
              End=start_val, Gain=list(), Loss=list(),"+Loss"=list()
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
        
        h_[start] <- hash(
          End=end, Gain=list(), Loss=list(),"+Loss"=list()
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
          
          h_[prev_end] <- hash(
            End=end, Gain=list(), Loss=list(),"+Loss"=list()
          )
          h_[[as.character(prev_end)]][[type]] <- list(end)
          
        }
      }
      
      return(h_)
    }
    
    deleteIntersections <- function(h_) {
      start_vals <- sort(as.numeric(keys(h_)))
      for (start_val in start_vals) {

        ## If both gain and loss are > 0, then there is overlap
        ## If min_val is > 1, then there is duplicate overlap
        ## Only delete one gain for each loss or one loss for each gain
        gain <- length(h_[[as.character(start_val)]][['Gain']])
        loss <- length(h_[[as.character(start_val)]][['Loss']])
        ploss <-length(h_[[as.character(start_val)]][['+Loss']])
        if(gain == (ploss + loss)){
          ###delete section
          del(start_val, h_)
        }else{
          
          ##delete ploss first
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
          ##delete loss afterwards 
          min_val <- min(gain, loss)
          # Delete leading elements from Gain and pLoss lists
          if (min_val > 0) {
            
            for (i in 1:min_val) {
              h_[[as.character(start_val)]][['Gain']][[1]] <- NULL
              h_[[as.character(start_val)]][['Loss']][[1]] <- NULL
            }
          }
          
        }
        
        ##do this again at the end after everything is cleared 
        ##gain <- length(h_[[as.character(start_val)]][['Gain']])
        ##loss <- length(h_[[as.character(start_val)]][['Loss']])
        ##ploss <-length(h_[[as.character(start_val)]][['+Loss']])
        ##if(gain == (ploss+ loss)){
          ###delete section
          ##del(start_val, h_)
        ##}
        
        
      }
      return(h_)
    }
    
    mergeAdjacentSections <- function(h_) {
      
      getContiguousSection <- function(h__, section, orig_end) {
        # This is a recursive function that crawls the hash to build
        # contiguous sections
        
        sect_end <- section[['End']]
        if (has.key(as.character(sect_end), h__)) {
          
          # extend section only if it matches the original ending of the previous section
          # and also matches Type (e.g., Gain/Loss)
          if (orig_end %in% h__[[sect_end]][[section[['Type']]]]) {
            # delete item from list
            orig_end_index <- match(orig_end, h__[[sect_end]][[section[['Type']]]])
            h__[[sect_end]][[section[['Type']]]][[orig_end_index]] <- NULL
            
            # extend section
            section[['End']] <- h__[[sect_end]][['End']]
            
            if (length(h__[[sect_end]][['Gain']]) == 0
                & length(h__[[sect_end]][['Loss']]) == 0
                & length(h__[[sect_end]][['+Loss']]) == 0)
            {
              # delete entire section if no more gain or loss
              delete(sect_end, h__)
            }
            
            # continue searching recursively
            section <- getContiguousSection(h__, section, orig_end)
            
          } 
          
        }
        return(section)
      }
      
      coord_table <- data.frame(matrix(ncol=4, nrow=0))
      colnames(coord_table) <- c('Chr', 'Start', 'End', 'Type')
      
      for (chr in keys(h_)) {
        # while there are still values in the section hash
        start_vals <- sort(as.numeric(keys(h_[[chr]])))
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
          } else if (length(h_[[chr]][[ch_start]][['+Loss']]) > 0)
          {
            type<-'+Loss'
          }else{
            # this should never happen
            print('Error, should have deleted this key')
          }
          
          # get original end value and delete from list
          orig_end <- h_[[chr]][[ch_start]][[type]][[1]]
          h_[[chr]][[ch_start]][[type]][[1]] <- NULL
          
          # delete section if no more gain or ploss
          if (length(h_[[chr]][[ch_start]][['Gain']]) == 0
              & length(h_[[chr]][[ch_start]][['+Loss']]) == 0 
              & length(h_[[chr]][[ch_start]][['Loss']]) == 0) {
            delete(ch_start, h_[[chr]])
          }
          
          
          # first part of contiguous section
          section <- c(Chr=chr, Start=start, End=end, Type=type)
          # build the contiguous section recursively
          section <- getContiguousSection(h_[[chr]], section, orig_end)
          
          ### add contiguous section to the table
          # find existing sections that fit at the end of the new section
          end_index <- which(
            as.numeric(coord_table[['Start']])==end & coord_table[['Type']]==type
          )
          # or existing sections that fit at the beginning of the new section
          start_index <- which(
            as.numeric(coord_table[['End']])==start & coord_table[['Type']]==type
          )
          if (length(start_index) > 0) {
            if (length(end_index) > 0) {
              # modify existing row to combine with new section and
              # existing end section
              coord_table[start_index[1],][['End']] <- coord_table[end_index[1],][['End']]
              # delete existing end section
              coord_table <- coord_table[-end_index[1],]
            } else {
              # modify existing row to combine with new section at the end
              coord_table[start_index[1],][['End']] <- section[['End']]
            }
          } else {
            if (length(end_index) > 0) {
              # modify existing row to combine with new section at the beginning
              coord_table[end_index[1],][['Start']] <- section[['Start']]
            } else {
              # add new section to table
              coord_table <- rbind(coord_table, as.data.frame(t(section)))
            }
          }
          
          # update list of start values for loop's stopping condition
          start_vals <- sort(as.numeric(keys(h_[[chr]])))
        }
      }
      return(coord_table)
    }
    
    ## Create an empty hash map for each chromosome
    chrs <- unique(M[,1])
    h <- hash()
    for (chr in chrs) {
      h[chr] <- hash()
    }
    
    ## Populate the hash map with sections
    for (row in 1:nrow(M)) {
      if (M[row, 'Type'] == 'Gain' | M[row, 'Type'] == 'Loss' | M[row, 'Type'] == '+Loss' ) {
        h[M[row, 'Chr']] <- insertSection(
          h[[M[row, 'Chr']]],
          as.numeric(M[row, 'Start']),
          as.numeric(M[row, 'End']),
          M[row, 'Type']
        )
      }
    }
    
    ## Delete intersections
    for (key in keys(h)) {
      h[key] <- deleteIntersections(h[[key]])
    }
    
    ## Merge adjacent sections and create final table
    final_coord_table <- mergeAdjacentSections(h)
    
    ##if want other info 
    if(keep_extras){
      final_coord_table<-rbind(final_coord_table,M_temp)
    }
    
    return(final_coord_table)
  }
  
  
  
  ##section to cancel out loss gains
  ### First code that takes two intervals (start, loss), where gain is first, that 
  ### returns a 3-column data frame with
  ### the intersection missing. 
  ### Can have 0, 1, or two rowsAlso has a column for type.
  
  mergeGL<-function(v1,v2){
    labs<-c("Gain", "Loss")
    v1<-as.vector(as.integer(as.numeric(as.character(v1[1:2]))))
    v2<-as.vector(as.integer(as.numeric(as.character(v2[1:2]))))
    out<-as.data.frame(matrix(nrow=0,ncol=3))
    colnames(out)<-c("Start","End","Type")
    out<-cbind(as.integer(as.numeric(out[,1])),as.integer(as.numeric(out[,2])),out[,3])
    if(v1[1]>v2[1]){
      labs<-labs[2:1]
      tmp<-v1
      v1<-v2
      v2<-tmp}
    
    if(v1[2]>=v2[2]){
      out<-data.frame(Start=c(v1[1],v2[2]), End=c(v2[1], v1[2]), 
                      Type=rep(labs[1],2))} else{if(v1[2]>=v2[1]){
                        out<-data.frame(Start=c(v1[1],v1[2]), End=c(v2[1], v2[2]), 
                                        Type=labs)} else{out<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                                                                         Type=labs)}}
    
    bads<-which(out[,1]>=out[,2])
    if(length(bads)>0){out<-out[-bads,]}
    origMat<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                        Type=labs)
    if(nrow(out)==nrow(origMat) && identical(out,origMat))
    {
      out=list(out,TRUE)
    }else{
      out=list(out,FALSE)
    }
    return(out)}
  
  ### Now one that takes two two-column matrices of gains and losses
  ### (gains first) and returns a 3-column data frame of merged
  
  mergeGLmat<-function(G, L){
    i<-1
    ##G[,1:2]<-apply(G[,1:2],2,as.numeric)
    ##L[,1:2]<-apply(L[,1:2],2,as.numeric)
    
    ##sort stuff
    L<-L[order(L[,1],(-L[,2])),]
    G<-G[order(G[,1],(-G[,2])),]
    
    origL<-L
    while(nrow(L)>0 & i<=nrow(G)){
      newL<-matrix(ncol=3, nrow=0)
      j<-1
      modified<-TRUE
      print(list(i,modified))
      while(j<=nrow(L) & modified){
        nxt<-mergeGL(G[i,], L[j,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        w<-which(nxt[,3]=="Loss")
        if(length(w)>0){newL<-rbind(newL, nxt[w,],
                                    if(nrow(L) >1 && !is.vector(L[-j,]) && nrow(L[-j,])>=j && !modified){L[-j,][j:nrow(L[-j,]),]}
        )
        }else if(nrow(nxt)>0){newL<-rbind(newL,if(nrow(L)>1 && !is.vector(L[-j,]) && nrow(L[-j,])>=j && !modified){L[-j,][j:nrow(L[-j,]),]})}else if(!is.vector(L[-j,])){newL<-rbind(newL,apply(L[-j,],2,as.character))}else{newL<-rbind(newL,sapply(L[-j,],as.character))}
        print(list(j,modified,newL))
        j<-j+1
      }
      L<-newL
      i<-i+1
    }
    
    
    i<-1
    modL<-origL
    while(nrow(G)>0 & i<=nrow(modL)){
      newG<-matrix(ncol=3, nrow=0)
      j<-1
      modified<-TRUE
      ##print(list(i,G))
      
      while(j<=nrow(G) & modified){
        nxt<-mergeGL(G[j,], modL[i,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        
        ##if returns as loss, collapse
        if(length(nxt[,3])> 0 && nxt[,3]=="Loss"){
          
          modL[i,]<-nxt
          modified=T
          ##newG<-rbind(newG,G[j,])
        }
        
        w<-which(nxt[,3]=="Gain")
        if(length(w)>0){newG<-rbind(newG, nxt[w,],
                                    if(nrow(G)>1 && !is.vector(G[-j,]) && nrow(G[-j,])>=j &&!modified ){G[-j,][j:nrow(G[-j,]),]}
        )
        }else if(nrow(nxt)>0){newG<-rbind(newG,if(nrow(G)>1 && !is.vector(G[-j,]) && nrow(G[-j,])>=j && !modified){G[-j,][j:nrow(G[-j,]),]})}
        else if(!is.vector(G[-j,])){newG<-rbind(newG,apply(G[-j,],2,as.character))}
        else{newG<-rbind(newG,sapply(G[-j,],as.character))}
        
        print(list(j,modified,newG))
        
        j<-j+1
      }
      G<-newG
      i<-i+1
    }
    
    G[,1]<-as.numeric(as.character(G[,1]))
    G[,2]<-as.numeric(as.character(G[,2]))
    L[,1]<-as.numeric(as.character(L[,1]))
    L[,2]<-as.numeric(as.character(L[,2]))
    
    
    
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
      if(is.vector(Msub))
      {
        out<-rbind(out, Msub)
        
      }else{
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
    }
    
    return(out)}
  
  ### This function takes a table M (a data frame) with columns:
  ### chromosome, start, stop, and type (Gain or Loss)
  ## and does proper merging
  
  ### First we need a function that takes two intervals (start, loss) that 
  ### returns NA if they don't overlap and merges them if they do.
  ### v1 and v2 are each two-vectors giving the interval  
  
  mergeInt<-function(v1,v2){
    v1<-as.integer(as.character(v1))
    v2<-as.integer(as.character(v2))
    if(is.na(prod(c(v1,v2)))){return(NA)}
    if(v1[1]>v2[1]){
      tmp<-v1
      v1<-v2
      v2<-tmp}
    
    if(v1[2]+1<v2[1]){return(NA)}
    
    return(c(v1[1], max(c(v1[2], v2[2]))))}
  
  ##returns overlaped region, deletes unoverlaps
  mergeIntOverlap<-function(v1,v2){
    v3=NA
    if(is.na(prod(c(v1,v2)))){return(NA)}
    if(v1[1] > v2[1]| v2[2] < v1[2] | v1[2] > v2[1] | v2[2] < v1[1]){
      
      v3<-sort(c(v1,v2))
      v3<-v3[c(2,3)]
      
    }
    
    ##if(v1[2]+1<v2[1]){return(NA)}
    
    return(v3)
  }
  
  
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
  
  
  ###use machinery for deletion intervals
  
  mergeDeletions<-function(M,Mainchr){
    OldM<-M
    startDel<-grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",M[,4])[1]
    
    
    wdel<-union(intersect(grep("del", M[, 4]),grep("^del", M[, 4],invert=T)),intersect(grep("add", M[, 4]),grep("^add", M[, 4],invert=T)))
    
    if(!is.na(startDel))
    {
      startDel<-1:startDel  
      wdel<-setdiff(union(intersect(grep("del", M[, 4]),grep("^del", M[, 4],invert=T)),intersect(grep("add", M[, 4]),grep("^add", M[, 4],invert=T))),startDel)
      
    }
    
    resttosee<-setdiff(grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",M[,4]),wdel)
    rest<-(1:nrow(M))[-1*c(wdel,resttosee)]
    
    
    
    M<-M[c(wdel,resttosee),]
    if(length(resttosee)>0 & length(wdel)>0)
    {
      M<- bigDelMerge(M)
    }
    
    startDel<-grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",M[,4])[1]
    
    
    wdel<-union(intersect(grep("del", M[, 4]),grep("^del", M[, 4],invert=T)),intersect(grep("add", M[, 4]),grep("^add", M[, 4],invert=T)))
    
    if(!is.na(startDel))
    {
      startDel<-1:startDel  
      wdel<-setdiff(union(intersect(grep("del", M[, 4]),grep("^del", M[, 4],invert=T)),intersect(grep("add", M[, 4]),grep("^add", M[, 4],invert=T))),startDel)
      
    }
    wg<-setdiff(grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",M[,4]),wdel)
    Mg<-M[wg,1:3]
    Ml<-M[wdel,1:3]
    
    if(length(wg)>1){
      ##Mg<-M[wg,1:3]
      ##for(j in 1:(nrow(Mg)-1)){
      ##Mg[j:nrow(Mg),]<-mergeBelow(Mg[j:nrow(Mg),])}
      M[wg,1:3]<-Mg}
    
    if(length(wdel)>1){
      ##Ml<-M[wl,1:3]
      ##for(j in 1:(nrow(Ml)-1)){
      ##Ml[j:nrow(Ml),]<-mergeBelowL(Ml[j:nrow(Ml),])}
      M[wdel,1:3]<-Ml}
    
    ##bads<-which(!M[,4]%in%c("Loss", "Gain"))
    
    ##deletions that are to be removed
    
    
    w<-which(!is.na(M[,2]))
    if(length(rest)>0)
    {
      out<-rbind(M[w,],OldM[rest,])
      w<-c(1:2)
    }else{
      
      out<-M[w,]
      
    }
    
    if(length(w)==1 ){out<-t(out)}
    
    return(out)}
  
  
  
  ##variant for deletions only
  
  
  mergeDelmat<-function(G, L){
    i<-1
    ##G[,1:2]<-apply(G[,1:2],2,as.numeric)
    ##L[,1:2]<-apply(L[,1:2],2,as.numeric)
    origL<-L
    while(nrow(L)>0 & i<=nrow(G)){
      newL<-matrix(ncol=3, nrow=0)
      j<-1
      modified<-TRUE
      print(list(i,modified))
      while(j<=nrow(L) & modified){
        nxt<-mergeDel(G[i,], L[j,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        w<-union(intersect(grep("del", nxt[, 3]),grep("^del", nxt[, 3],invert=T)),intersect(grep("add", nxt[, 3]),grep("^add", nxt[,3],invert=T)))
        if(length(w)>0){newL<-rbind(newL, nxt[w,],
                                    if(nrow(L) >1 && !is.vector(L[-j,]) && nrow(L[-j,])>=j && !modified){L[-j,][j:nrow(L[-j,]),]}
        )
        }else if(nrow(nxt)>0){newL<-rbind(newL,if(nrow(L)>1 && !is.vector(L[-j,]) && nrow(L[-j,])>=j && !modified){L[-j,][j:nrow(L[-j,]),]})}else if(!is.vector(L[-j,])){newL<-rbind(newL,apply(L[-j,],2,as.character))}else{newL<-rbind(newL,sapply(L[-j,],as.character))}
        print(list(j,modified,newL))
        j<-j+1
      }
      L<-newL
      i<-i+1
    }
    
    
    i<-1
    while(nrow(G)>0 & i<=nrow(origL)){
      newG<-matrix(ncol=3, nrow=0)
      j<-1
      modified<-TRUE
      print(list(i,G))
      
      while(j<=nrow(G) & modified){
        nxt<-mergeDel(G[j,], origL[i,])
        modified<-nxt[[2]]
        nxt<-nxt[[1]]
        w<-setdiff(grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",nxt[,3]),union(intersect(grep("del", nxt[, 3]),grep("^del", nxt[, 3],invert=T)),intersect(grep("add", nxt[, 3]),grep("^add", nxt[,3],invert=T))))
        if(length(w)>0){newG<-rbind(newG, nxt[w,],
                                    if(nrow(G)>1 && !is.vector(G[-j,]) && nrow(G[-j,])>=j &&!modified ){G[-j,][j:nrow(G[-j,]),]}
        )
        }else if(nrow(nxt)>0){newG<-rbind(newG,if(nrow(G)>1 && !is.vector(G[-j,]) && nrow(G[-j,])>=j && !modified){G[-j,][j:nrow(G[-j,]),]})}else if(!is.vector(G[-j,])){newG<-rbind(newG,apply(G[-j,],2,as.character))}else{newG<-rbind(newG,sapply(G[-j,],as.character))}
        
        print(list(j,modified,newG))
        j<-j+1
      }
      G<-newG
      i<-i+1
    }
    
    G[,1]<-as.numeric(as.character(G[,1]))
    G[,2]<-as.numeric(as.character(G[,2]))
    L[,1]<-as.numeric(as.character(L[,1]))
    L[,2]<-as.numeric(as.character(L[,2]))
    
    
    
    out<-data.frame(Start=c(G[,1], L[,1]), End=c(G[,2],L[,2]), 
                    Type=c(G[,3], L[,3]))
    return(out)}
  
  ### Now the whole thing:
  bigDelMerge<-function(M){
    
    ##ask what does this do 
    out<-M[-(1:nrow(M)),]
    
    chrs<-unique(M[,1])
    
    for(i in 1:length(chrs)){
      w<-which(M[,1]==chrs[i])
      Msub<-M[w,]
      
      startDel<-grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\())",Msub[,4])[1]
      
      
      wl<-union(intersect(grep("del", Msub[, 4]),grep("^del", Msub[, 4],invert=T)),intersect(grep("add", Msub[, 4]),grep("^add", Msub[, 4],invert=T)))
      
      if(!is.na(startDel))
      {
        startDel<-1:startDel  
        wl<-setdiff(union(intersect(grep("del", Msub[, 4]),grep("^del", Msub[, 4],invert=T)),intersect(grep("add", Msub[, 4]),grep("^add", Msub[, 4],invert=T))),startDel)
        
      }
      wg<-setdiff(grep("((t\\()|(idic\\()|(rob\\()|(trc\\()|(dic\\()|(Gain))",Msub[,4]),wl)
      
      if(length(wg)==0 | length(wl)==0){
        out<-rbind(out, Msub)} else{
          nxt<-mergeDelmat(Msub[wg,-1],Msub[wl,-1])
          if(nrow(nxt)>0)
          {
            nxt<-cbind(chrs[i],nxt)
            colnames(nxt)[1]<-"Chr"
          }
          out<-rbind(out, nxt)
        }
    }
    
    return(out)}
  
  
  mergeDel<-function(v1,v2){
    labs<-c(v1[3], v2[3])
    v1<-as.vector(as.integer(as.numeric(as.character(v1[1:2]))))
    v2<-as.vector(as.integer(as.numeric(as.character(v2[1:2]))))
    out<-as.data.frame(matrix(nrow=0,ncol=3))
    colnames(out)<-c("Start","End","Type")
    out<-cbind(as.integer(as.numeric(out[,1])),as.integer(as.numeric(out[,2])),out[,3])
    if(v1[1]>v2[1]){
      labs<-labs[2:1]
      tmp<-v1
      v1<-v2
      v2<-tmp}
    
    if(v1[2]>=v2[2]){
      out<-data.frame(Start=c(v1[1],v2[2]), End=c(v2[1], v1[2]), 
                      Type=rep(labs[1],2))} else{if(v1[2]>=v2[1]){
                        out<-data.frame(Start=c(v1[1],v1[2]), End=c(v2[1], v2[2]), 
                                        Type=labs)} else{out<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                                                                         Type=labs)}}
    
    bads<-which(out[,1]>=out[,2])
    if(length(bads)>0){out<-out[-bads,]}
    origMat<-data.frame(Start=c(v1[1],v2[1]), End=c(v1[2], v2[2]), 
                        Type=labs)
    if(nrow(out)==nrow(origMat) && identical(out,origMat))
    {
      out=list(out,TRUE)
    }else{
      out=list(out,FALSE)
    }
    return(out)}
  
  
  
  
  ##parseing per sample
  
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
    xdel = 0 ##counts if -X occures as constitutional
    ydel = 0 ##counts if -Y occures as consitutional 
    
    xconstitutional=0 ##shift counts for xc indications (kind of a correction factor)
    yconstitutional=0  ##shift counts for yc indications (kind of a correction factor)
    
    idealx=0 ##estimate of what the x value should be
    idealy=0 ##estimate of what the y value should be
    
  
    addtot = 0 ##counts total "new chromosomes"
    deltot = 0 ##counts total complete chrom deletions
    modtot = 0 ##for idems only, counts modification chromosomes
    
    n = 1 ##ploidy count
    ploidy=2 ## ploidy non additive ##default 2 for diploid
    
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
        if(grepl("(c$)|(c\\?$)",Cyto_sample[2]))
        {
          
          xconstitutional=str_count(Cyto_sample[2], "X") 
          yconstitutional=str_count(Cyto_sample[2], "Y") 
          
        }
        
        #if(constitutional==F & grepl("(c$)|(c\\?$)",Cyto_sample[2])){
        #  
        #  normX = str_count(Cyto_sample[2], "X")
        #  normY = str_count(Cyto_sample[2], "Y")
          
        #}else 
        if (any(grepl("Y", Cyto_sample)) &
            (sum(grepl("Y", Cyto_sample), na.rm = TRUE) > sum(grepl("\\+Y", Cyto_sample), na.rm =
                                                              TRUE)))
        {
          normX = 1
          normY = 1
        } ##if ? is a chromosome and sexstimate is off, dont make guesses on constitutional change
        else if(( sexstimate==F && any(grepl("\\?",Cyto_sample[2])) )){
          normX=2
          normY=0
        }else
        {
          normX = 2
        }
        
        ##count number of XY in 2nd slot ##make sure 2nd slot is sex chromosomes, change this code later
        if (!grepl("^[XY]+", Cyto_sample[2]))
        {
          #make sure reg ecpression means what you want it to
        }  else {
          xcount = str_count(Cyto_sample[2], "X")
          xcount = xcount + str_count(Cyto_sample[2], "\\?")
          ycount = str_count(Cyto_sample[2], "Y")
          ##x,y calculation
        }
      }
      
      ##check for ploidy levels and digit is over 2
      if (grepl("<[[:digit:]]n>", Cyto_sample[1]))
      {
        ##extract range of stuff before n
        n = as.numeric(strsplit(Cyto_sample[1], "<|n>")[[1]][2])
        n=ploidy
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
        
        ##for if deletions get completely eliminated in t() type abberations
        deletions=F
        
        #########
        ##if guess is true, try to process ? marks
        ##think about how this can affect counting  + and \\?
        ######new######
        #############################
        if(guess_q == T)
        {
          Cyto_sample[j] <- gsub("\\?","",Cyto_sample[j])
        }
        
        ##if or, delete second option based on booleen
        if(orOption==T){
          Cyto_sample[j]<-gsub("or.*$","",Cyto_sample[j])
        }
        
        #addition and deletions of entire chromosomes can be handeled easily and separate from the rest
        if (grepl(
          "mar|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*r\\(*[[:digit:]]*\\)*$|^\\+*([[:digit:]]((~|-)[[:digit:]])*)*neo[[:digit:]]*$",
          Cyto_sample[j]
        ) || (constitutional==F && grepl("(c$)|(c\\?$)",Cyto_sample[j])))
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
              if(constitutional==F & grepl("\\+[[:digit:]]+c\\?*$", Cyto_sample[j]))
              {
                ##just one addition
                tem<-1 
              }else{
                if(!grepl("\\+[[:digit:]]+c\\?*$", Cyto_sample[j]))
                  {
                  tem <-
                    as.numeric((strsplit(
                      strsplit(Cyto_sample[j], "mar|r\\(*[[:digit:]]*\\)*$|neo|c$|c\\?$")[[1]][1],
                      "\\+"
                    )[[1]][2]))
                }
              }
            }
            
            ##only if tem is a number 
            if(is.numeric(tem)||is.integer(tem)||is.double(tem))
            {
              addtot <- addtot + tem
            }else{
              ##output an error
              Dump_table <- rbind(Dump_table, c(Con_data[i,], "Error in markers and other ambiguous objects not accounted for"))
              
            }
            
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
        } else if (grepl("^-[[:digit:]]+c*$", Cyto_sample[j]) | grepl("-X", Cyto_sample[j]) | grepl("-Y", Cyto_sample[j]))
        {
          cytoName <- gsub("c", "", substring(Cyto_sample[j], first = 2))
          chr_name <-
            ref_table[grep(paste("chr", as.character(cytoName), "$", sep = ""), ref_table), ]
         
           ##exclude x and y deletions until the end
           if(!grepl("Y",chr_name[1]) & !grepl("X",chr_name[1])){
            temp_table[1, 1] = chr_name[1]
            temp_table[1, 2] = "0"
            temp_table[1, 3] = chr_name[2]
            temp_table[1, 4] = "Loss"
          }
          
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
          inc_table <- colparse(j, xmod, ymod, transloctable,addtot,Cyto_sample)
          if (is.character(inc_table))
          {
            Dump_table <- rbind(Dump_table, c(Con_data[i, ], inc_table))
          }
          else if (is.null(inc_table)) {
            
          }
          else if(length(inc_table)==1 && is.na(inc_table)){
            Dump_table <- rbind(Dump_table, c(Con_data[i,], "Error in more than one band associated with a chromosome in a translocation"))
            
          }
          else{
            temp_table <- inc_table[[1]]
            original_temp_table<-inc_table[[1]]
            ex_table <- inc_table[[2]]
            original_ex_table<-inc_table[[2]]
            xmod <- inc_table[[3]]
            ymod <- inc_table[[4]]
            Mainchr <- inc_table[[5]]
            multi <- inc_table[[6]]
            transloctable<-inc_table[[7]]
            addtot<-inc_table[[8]]
            
            print(multi)
            ##for future implementation: switch to searching at ex table if temp table is empty, default is temp_table
            pointerConditional<-temp_table
            
            if (nrow(temp_table) > 0)
            {
              ##if its a translocated vector, fix now
              if(nrow(temp_table)==4 & ncol(temp_table)==1)
              {
                temp_table<-t(temp_table)
                colnames(temp_table)<-c("Chr","Start","End","Type")
              }
              
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
              }else if(multi>2){
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
                    tempindexchromtracker<-sapply(gsub("\\$", "", Mainchr), function(x) {if(!grepl("X|Y",x))
                      which.max(clone_chrom_tracker == as.numeric(x))
                    })
                    
                    ##fix up if its a list
                    if(is.list(tempindexchromtracker))
                    {
                      tempindexchromtracker<-unlist(tempindexchromtracker)
                    }
                    
                    if(!is.null(tempindexchromtracker) && length(tempindexchromtracker)>0 )
                    {
                      clone_chrom_tracker <-
                        clone_chrom_tracker[-1 *tempindexchromtracker ]
                    }
                  }
                }
                
                ##check how this handles X chromosomes
                ## if more than 2 chromosomes involved and its not a translocation an insertion or an inversion, something is "deleted"
                if (any(grepl("-[:alpha:]", temp_table[, 4])) &
                        any(!grepl("X\\$|Y\\$", Mainchr)))
                {
                  ##deltot <- deltot + (1 * multi)
                }else
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
              if(any(grepl("((t)|(ins))\\(",temp_table[,4])&!grepl("^((t)|(ins))\\(",temp_table[,4])))
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
                  )), grep("^((t)|(ins))\\(",temp_table[,4],invert=T)
                ), 4] <- "Gain"
                if(nrow(ex_table)>0)
                {
                  ex_table[intersect(intersect(
                    grep(paste("chr", Mainchr, sep = "",collapse="|"), ex_table[, 1]),
                    grep(
                      "((t)|(ins))\\(",
                      ex_table[, 4]
                    )), grep("^((t)|(ins))\\(",ex_table[,4],invert=T)
                  ), 4] <- "Loss"
                }
                ##temp_table[, 4] <-
                ##  paste(additionTable[[1]], temp_table[, 4], sep = "")
                ##ex_table[, 4] <-
                ##  paste(additionTable[[2]], ex_table[, 4], sep = "")
                ##temp_table <- rbind(temp_table, ex_table)
              }
              
              
              ##take deletion areas out of translocation adjacent things
              if ( ((any(grepl("del", temp_table[, 4])&!grepl("^del", temp_table[, 4])&!grepl("multi[[:digit:]]*del",temp_table[,4]))) | (any(grepl("add", temp_table[, 4])&!grepl("^add", temp_table[, 4])&!grepl("multi[[:digit:]]*add",temp_table[,4])))))
              {
                if(nrow(temp_table)>length(union(((intersect(grep("del", temp_table[, 4]), intersect(grep("^del", temp_table[, 4],invert=T),grep("multi[[:digit:]]*del",temp_table[,4],invert=T))))) ,(intersect(grep("add", temp_table[, 4]), intersect(grep("^add", temp_table[, 4],invert=T),grep("multi[[:digit:]]*add",temp_table[,4],invert=T)))))) )
                {
                  temp_table<-mergeDeletions(temp_table,Mainchr)
                  deletions=T
                  ##print(temp_table)
                  ##if it is a vector, convert 
                  if(is.vector(temp_table))
                  {
                    temp_table<-as.data.frame(as.list(temp_table))
                    temp_table[, 2:3] <-as.numeric(temp_table[,2:3])
                    colnames(temp_table)<-c("Chr","Start","End","Type")
                  }else if(nrow(temp_table)==4 & ncol(temp_table)==1)
                  {
                    ##if its a translocated vector, fix now
                    temp_table<-t(temp_table)
                    temp_table<-as.data.frame(as.list(temp_table))
                    temp_table[, 2:3] <- as.numeric(temp_table[, 2:3])
                    colnames(temp_table)<-c("Chr","Start","End","Type")
                  }
                  
                  additionTable <- detectAdd(temp_table[, 4], if(nrow(ex_table)>0){ex_table[, 4]}else{NULL})
                  
                }
                 
         
                
              }
              
              
              ##if temp table is empty, skip all this
              if(nrow(temp_table)>0)
              {
                
              }
                    if (  (any(grepl(
                      "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                      temp_table[, 4]
                    )) & (length(temp_table[,4])>length(c(which(grepl("del",temp_table[,4])),which(grepl("add",temp_table[,4])))  )) & (length(temp_table[,4])>1 & ((length(temp_table[,4])-length(c(which(grepl("del",temp_table[,4])),which(grepl("add",temp_table[,4])))))>1 ) ))   || (any(grepl(
                      "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                      temp_table[, 4]
                    ))
                      && deletions)
                    )
                    {
                      mod = TRUE
                      if(nrow(temp_table)>0)
                      {
                        temp_table[intersect(
                          grep(
                            paste(paste("chr", Mainchr, sep = ""),collapse = "|"),
                            temp_table[, 1],
                            invert = T
                          ),
                          grep(
                            "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                            temp_table[, 4]
                          )
                        ), 4] <- "Gain"
                      }
                      
                      if(nrow(ex_table)>0)
                      {                
                        ex_table[intersect(
                          grep(paste(paste("chr", Mainchr, sep = ""),collapse="|"), ex_table[, 1]),
                          grep(
                            "((^\\++((der)|(rec))\\(.*)|(^((der)|(rec))\\(.*))",
                            ex_table[, 4]
                          )
                        ), 4] <- "Loss"
                      
                      }
                      
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      
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
                      if(nrow(temp_table)>0)
                      {
                        original_temp_table<-temp_table
                        temp_table[grep("i\\(.*", temp_table[, 4]), 4] <- "Gain"
                        additionTable[[1]] <- c(additionTable[[1]],additionTable[[1]][grep("i\\(.*", original_temp_table[, 4])]) 
                        temp_table<-rbind(temp_table,original_temp_table[grep("i\\(.*", original_temp_table[, 4]), ])
                        
                        ##add to additionTable in light of new addition
                      }
                      if(nrow(ex_table)>0)
                      {
                        ex_table[grep("i\\(.*", ex_table[, 4]), 4] <- "Loss"
                      }
                      if(nrow(temp_table)>0)
                      {
                       temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
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
                      if(nrow(temp_table)>0)
                      {
                        original_temp_table<-temp_table
                        temp_table[intersect(grep("ider\\(.*", temp_table[, 4]),grep("del\\(.*|add\\(.*",temp_table[,4],invert=T)), 4] <-
                        "Gain"
                        additionTable[[1]] <- c(additionTable[[1]],additionTable[[1]][intersect(grep("ider\\(.*", temp_table[, 4]),grep("del\\(.*|add\\(.*",temp_table[,4],invert=T))])
                        temp_table<-rbind(temp_table,original_temp_table[grep("ider\\(.*", original_temp_table[, 4]), ])
                        
                        
                        ##activate deletions ? Think about this 
                        ##JAN Thinking time 
                        ######
                        #######
                      }
                      if(nrow(ex_table)>0)
                      {
                        ex_table[grep("ider\\(.*", ex_table[, 4]), 4] <- "Loss"
                      }
                      
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      if(nrow(ex_table)>0)
                      {
                        ex_table[, 4] <-
                          paste(additionTable[[2]], ex_table[, 4], sep = "")
                      }
                      temp_table <- rbind(temp_table, ex_table)
                      
                      ##mark deletions now
                      ##temp_table[intersect(grep("ider\\(.*", temp_table[,4]),grep("del\\(.*|add\\(.*",temp_table[,4])),4] <-"Loss"
                      
                      
                    }
                    
                    if (any(grepl("idic\\(.*", temp_table[, 4])))
                    {
                      mod = TRUE
                      ##need to either break it up or replace,
                      ##think about this one too
                      ##deleted <-
                      ##ex_table[grep("idic\\(.*&del\\(", ex_table[, 4]), 4]
                      if(nrow(temp_table)>0)
                      {
                        
                        original_temp_table<-temp_table
                        temp_table[intersect(grep("idic\\(.*", temp_table[, 4]),grep("dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", temp_table[, 4],invert=T)  ), 4] <- "Loss"
                        additionTable[[1]] <- c(additionTable[[1]],additionTable[[1]][intersect(grep("del\\(.*|add\\(.*", original_temp_table[, 4],invert=T),grep("idic\\(.*",original_temp_table[,4]))])
                        temp_table<-rbind(temp_table,original_temp_table[intersect(grep("del\\(.*|add\\(.*", original_temp_table[, 4],invert=T),grep("idic\\(.*",original_temp_table[,4])), ])
                      }
                      if(nrow(ex_table)>0)
                      {
                        original_ex_table<-ex_table
                        ex_table[grep("idic\\(.*", ex_table[, 4]), 4] <- "Gain"
                        additionTable[[2]] <- c(additionTable[[2]],additionTable[[2]][grep("idic\\(.*", original_ex_table[, 4])]) 
                        ex_table<-rbind(ex_table,original_ex_table[grep("idic\\(.*", original_ex_table[, 4]), ])
                        
                      }
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      if(nrow(ex_table)>0)
                      {
                        ex_table[, 4] <-
                          paste(additionTable[[2]], ex_table[, 4], sep = "")
                      }
                      temp_table <- rbind(temp_table, ex_table)
                      ##mark deletions now
                      ##temp_table[intersect(grep("idic\\(.*", temp_table[,4]),grep("del\\(.*|add\\(.*",temp_table[,4])),4] <-"Loss"
                      
                    }
                    
                    ##figour out why you did this 
                    ##make sure this is reversed for long form
                    if (any(which(grepl("dic\\(.*", temp_table[, 4]) & !grepl("idic\\(.*", temp_table[, 4]) )))
                    {
                      mod = TRUE
                      if(any(grepl("long",temp_table[,4])))
                      {
                        if(nrow(temp_table)>0)
                        {
                          temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                        }
                        if(nrow(ex_table)>0)
                        {
                          ex_table[grep("dic\\(.*", ex_table[, 4]) ,4]<-"Loss"
                          ex_table[, 4] <-
                            paste(additionTable[[2]], ex_table[, 4], sep = "")
                        }
                        
                      }else
                      {
                        if(nrow(temp_table)>0)
                        {
                          temp_table[intersect(grep("dic\\(.*", temp_table[, 4]),grep("idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", temp_table[, 4],invert=T)  ) , 4] <- "Loss"
                          
                          temp_table[, 4] <-
                            paste(additionTable[[1]], temp_table[, 4], sep = "")
                        }
                        if(nrow(ex_table)>0)
                        {
                          
                          ex_table[intersect(intersect(grep("dic\\(.*", ex_table[, 4]),grep("idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", ex_table[, 4],invert=T)  ),
                                             grep("del\\(.*|add\\(.*",ex_table[,4])) , 4] <- gsub("del\\(.*|add\\(.*","",
                                                                                                  ex_table[intersect(intersect(grep("dic\\(.*", 
                                                                                                  ex_table[, 4]),grep("idic\\(.*|dup\\(.*|tan\\(.*|trp\\(*.|qdq\\(.*", ex_table[, 4],invert=T)  ),
                                                                                                  grep("del\\(.*|add\\(.*",ex_table[,4])) , 4]) 

                          

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
                        if(nrow(ex_table)>0)
                        {
                          ex_table[intersect(grep("trc\\(.*", ex_table[, 4]),grep(paste("chr",Mainchr[1],"chr",Mainchr[3],sep='|'),ex_table[,1])),4]<-"Loss"
                        }
                        
                      }else
                      {
                        if(nrow(temp_table)>0)
                        {
                          temp_table[intersect(grep("trc\\(.*", temp_table[, 4]),grep(paste("chr",Mainchr[1],"chr",Mainchr[3],sep='|'),temp_table[,1])),4]<-"Loss"
                        }
                      }
                      
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      
                      if(nrow(ex_table)>0)
                      {
                        ex_table[intersect(grep("trc\\(.*", ex_table[, 4]),grep("chr",Mainchr[2],ex_table[,1])), 4] <- "Loss"
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
                      
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
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
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
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
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      
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
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      
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
                      if(nrow(temp_table)>0)
                      {
                        temp_table[grep("del", temp_table[, 4]), 4] <-
                          "Loss"
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
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
                    
                    if (any(grepl("add", temp_table[, 4])))
                    {
                      mod = TRUE
                      if(nrow(temp_table)>0)
                      {
                        temp_table[grep("add", temp_table[, 4]), 4] <-
                          "Loss"
                        temp_table[, 4] <-
                          paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
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
                      if(nrow(temp_table)>0)
                      {
                        temp_table[, 4] <-
                        paste(additionTable[[1]], temp_table[, 4], sep = "")
                      }
                      if(nrow(ex_table)>0)
                      {
                        ex_table[, 4] <-
                          paste(additionTable[[2]], ex_table[, 4], sep = "")
                      }
                      temp_table <- rbind(temp_table, ex_table)
                    }
                    
                    
                    
                    ##handel stuff for minus as well
                    ##if +Gain 
                if(nrow(temp_table)>0)
                {
                    temp_table[grep("^-Gain|^-$", temp_table[, 4]), 4] <- "Loss"
                    temp_table[grep("\\+Loss",temp_table[,4]),4] <- "+Loss"
                    
                    temp_table[grep("\\+Gain", temp_table[, 4]), 4] <- "Gain"
                    temp_table[intersect(grep("\\+", temp_table[, 4]),grep("\\+Loss", temp_table[, 4],invert=T)), 4] <- "Gain"
                    ##check other permutations of this/ dont think insertions belong here
                    temp_table[grep("dup|qdp|tan|trp|\\+$", temp_table[, 4]), 4] <-
                      "Gain"
                    
                    Plus_Loss<-temp_table[grep("\\+Loss",temp_table[,4]),]
                    
                    ##cut off intersections early
                    ###Think aboiut this 
                    if(length(Plus_Loss)>0 && nrow(Plus_Loss)>0 && any(grepl("Gain",temp_table[,4]) & any(grepl("Loss",temp_table[,4]))))
                    {
                     ## temp_table[grep("\\+Loss", temp_table[, 4]), 4] <- "Loss"
                      
                      temp_table<-mergeTable_jp(temp_table)
                      ##indexLoss<--1
                      ##loop Plus Loss here
                      ##ask for more efficient methods
                      ##for(i in 1:length(Plus_Loss)){
                      ##  indexLoss<-c(indexLoss,match(Plus_Loss,temp_table[,1:3]))
                      ##}
                      ##indexLoss<-indexLoss[-1]
                      ##if(length(indexLoss)>0){
                      ##  temp_table[indexLoss,4]<-"+Loss"
                      ##}
                    
                    }
                    
                    temp_table[grep("\\+Loss", temp_table[, 4]), 4] <- ""
                    
                }
              
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
        val<-as.numeric(val)
        
        if(as.numeric(val[2])>=as.numeric(val[1]))
        {
          val = seq(from = val[1],
                    to = val[2],
                    by = 1)
        }else{
          val = seq(from = val[2],
                    to = val[1],
                    by = 1)
        }
      }
      else{
        val<-as.numeric(val)
        val = val[2]
      }
      
    }
    
    val = as.numeric(val)
    
    ##prelim ploidy gueeses for xideal and yideal from val
    
    
    
    
    ##take x and y completely out of this 
    val = val -xcount -ycount-(xadd+yadd-xdel-ydel)
    if(normX+normY != 2)
    {
      if(normX+normY>2)
      {
        val=val+(normX+normY-2)
      }else{
        val=val+(2-normX+normY)
      }
    }
    
    ##original val
    orgval = strsplit(Cyto_sample[1], "<")[[1]] ##numer of chromosomes indicated in first value
    orgval = paste(strsplit(orgval, "[[:alpha:]]+")[[1]],
                sep = "",
                collapse = "")
    if (grepl("~|-", orgval))
    {
      orgval = unlist(strsplit(orgval, "~|-"))
      if (all(!is.na(as.numeric(orgval))))
      {
       
          if (all(!is.na(as.numeric(orgval))))
          {
            orgval<-as.numeric(orgval)
            if(as.numeric(orgval[2])>=as.numeric(orgval[1])){
             orgval = seq(from = orgval[1],
                           to = orgval[2],
                           by = 1)
            }else{
              orgval = seq(from = orgval[2],
                           to = orgval[1],
                           by = 1)
            }
          }
      }
      else{
        orgval<-as.numeric(orgval)
        orgval = orgval[2]
      }
      
    }
    
    if (any(is.na(val)))
    {
      Dump_table <-
        rbind(Dump_table, c(as.vector(Con_data[i, ]), "Error in unclear chrom number"))
    } else
    {
      idealval = 46 + addtot - deltot- xcount - ycount -(xadd+yadd-xdel-ydel)
      if(normX+normY != 2)
      {
        if(normX+normY>2)
        {
          idealval=idealval+(normX+normY-2)
        }else{
          idealval=idealval+(2-normX+normY)
        }
      }
      
      diffval = val - idealval
      val_remainder = diffval %% 22
      val_remainder_2=diffval %% 23
      val_divider = diffval / 22 ## for stuff over diploidy
      val_divider_2=diffval/23
      val_temp=0
      
      if (!is.null(val_remainder))
      {
        if (any(val_remainder == 0))
        {
          val_divider <- val_divider[grep(TRUE, val_remainder == 0)]
          val_remainder <- val_remainder[grep(TRUE, val_remainder == 0)]
          
        }else if(any(val_remainder_2==0)){
          val_divider <- val_divider_2[grep(TRUE, val_remainder_2 == 0)]
          val_remainder<- val_remainder_2[grep(TRUE, val_remainder_2 == 0)]
        }else{
          val_remainder = val_remainder[1]
          val_divider = val_divider[1]
        }
        
        print(c(val, idealval, val_remainder, val_divider))
        print(c(diffval,orgval))
        
        
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
            
            ploidy=val_divider+2

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
              colnames(temp_table)<-colnames(sample_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              
              ploidy=1
              
              
              
            }
            
            print(c("val_div<0", Cyto_sample))
          }
       
        }else if (grepl("^ids$", Cyto_sample[length(Cyto_sample)]) && any(diffval / length(clone_chrom_tracker) > 0 &
                                                                          diffval %% length(clone_chrom_tracker) == 0) )
        {
          ##if unaccounted chromosom number == diff val, add them , let this estimate for uncertainty

          
            
            diffval<-diffval[intersect(grep(TRUE,diffval / length(clone_chrom_tracker) > 0),grep(TRUE,diffval %% length(clone_chrom_tracker) == 0))]
            ploidy=(diffval/length(clone_chrom_tracker))+2
            
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
              colnames(temp_table)<-colnames(sample_table)
              
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            }
            

          
        }else if(guess==T & (any(diffval< (-14))| (any((diffval<31&diffval>14)))|any((diffval<53&diffval>38))|any((diffval>56 &(floor(val_divider)*22+3<diffval) & (ceiling(val_divider)*23-3>diffval) )))){
          ##clone_chrom_tracker <- rep(1:22, 2)
          if(length(diffval) > 1)
          {
            ## take the one that fits and first one that is true
            new_diffval<-diffval[which(diffval< (-14)| (diffval<31&diffval>14)|(diffval<53&diffval>38)|(diffval>56 &(floor(val_divider)*22+3<diffval) & (ceiling(val_divider)*23-3>diffval) ))][1]
            
          }else{
            new_diffval<-diffval
          }
          ##think about how to implement this beyond teraploidy
          ##doesnt work
          print(new_diffval)
          if(new_diffval < (-14))
           {
     
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Loss", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
                ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                  ##              temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            ploidy=1
            
           }
           
           
          if(new_diffval<31&new_diffval>14)
          {
            
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
                ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                  ##              temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              ploidy=3
            
          }
          
          if(new_diffval<53&new_diffval>38)
          {
            
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
                ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                  ##              temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
                ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                  ##              temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            ploidy=4
          }
            
          if(new_diffval>56 & floor(val_divider)*22+3<new_diffval& ceiling(val_divider)*22-3>new_diffval)
          {
            end<-round(val_divider)
            ploidy=end+2
            for (k in 1:end)
            {
              temp_table <-
                data.frame(ref_table[, 1],
                            rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
                ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
                  ##                temp_table[, 1]), ]
                
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            }
              
          }  
          

          ##time to guess according to ranges if guessing is true
          
          ##put in dump table
          Dump_table <-
            rbind(Dump_table,
                  c(
                    as.vector(Con_data[i, ]),
                    "Warning in some chromosomes unaccounted for"
                  ))
          print(c("unaccounted", Cyto_sample))
          
        
        }else if((guess_by_first_val==T | forMtn==T) &
                 (any(orgval > 0 & orgval <= 34)|any(orgval >= 58 & orgval <= 80)
                  |any(orgval >= 81 & orgval <= 103)|any(orgval >= 104 & orgval <= 126)
                  |any(orgval >= 127 & orgval <= 149)|any(orgval >= 150 & orgval <= 172 )
                  |any(orgval >= 173 & orgval <= 195))){
          ###################################################
          ###################################################
          ##################################################
          ###THIS IS NEW
          ### DONT FORGET TO COMMIT TO DR PHANS SCRIPT
          ###DONT FORGET
          ##PAY ATTENTION
          ##ETC
          #ETC 
          ##ETC
          ################################################
          ##################################################
          ################################################
          
          ##if the initial exact match or the estimate match does not work, and the options are turned on, estimate from the initial value as defined in the iscn up until the octaploid level 
          if(length(orgval) > 1)
          {
            ## take the one that fits and first one that is true
            new_val<-orgval[which((orgval >= 0 & orgval <= 34)| (orgval>=58&orgval<=80)|
                                 (orgval>=81 & orgval <= 103)|(orgval>=104 & orgval<=126)|(orgval >= 127 & orgval <= 149)
                               |(orgval>= 150 & orgval <= 172)|(orgval >= 173 & orgval <= 195) )][1]
            
          }else{
            new_val<-orgval
          }
          print(c(new_val))
          
          ##assign ploidy here 
          if(new_val >= 0 & new_val <= 34){
            ploidy=1
          }
          if(new_val>=58 & new_val<=80){
            ploidy=3
          }
          if(new_val>=81 & new_val <= 103){
            ploidy=4
          }
          if(new_val>=104 & new_val<=126){
            ploidy=5
          }
          if(new_val >= 127 & new_val <= 149){
            ploidy=6
          }
          if(new_val>= 150 & new_val <= 172){
            ploidy=7
          }
          if(new_val >= 173 & new_val <= 195) {
            ploidy=8
          }
          ##loop to get values here
          if(ploidy == 1 )
          {
            
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Loss", nrow(ref_table)))
            temp_table<-temp_table[1:22,]
            ##temp_table <-
            ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
            ##              temp_table[, 1]), ]
            
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])

          }else if(ploidy > 2 )
          {
            end<-(ploidy-2)
            for (k in 1:end)
            {
              temp_table <-
                data.frame(ref_table[, 1],
                           rep(0, nrow(ref_table)),
                           ref_table[, 2],
                           rep("Gain", nrow(ref_table)))
              temp_table<-temp_table[1:22,]
              ##temp_table <-
              ##temp_table[grep(paste("chr", clone_chrom_tracker, collapse = '|'),
              ##                temp_table[, 1]), ]
              
              temp_table[, 4] <- as.character(temp_table[, 4])
              temp_table <- as.matrix(temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
              sample_table <- rbind(sample_table, temp_table)
              sample_table[, 4] <- as.character(sample_table[, 4])
            }
            
          }  
          
          print(c(new_val,ploidy))
          
          ##time to guess according to ranges if guessing is true
          
          ##put in dump table
          Dump_table <-
            rbind(Dump_table,
                  c(
                    as.vector(Con_data[i, ]),
                    "Warning in some chromosomes unaccounted for"
                  ))
          print(c("unaccounted", Cyto_sample))
          
        }else if(val_divider>1){
          ##if current ploidy cannot be estimated, see if the next decrease level in ploidy can be calculated
          val_temp=floor(val_divider)
          ploidy=val_temp+2
          if (val_temp > 0)
          {
            for (f in 1:(val_temp))
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
           
              }
            
              print(c("val_div>0", Cyto_sample))
              ##add this to uncertain 
              
              ##put in dump table
              Dump_table <-
                rbind(Dump_table,
                      c(
                        as.vector(Con_data[i, ]),
                        "Warning in some chromosomes unaccounted for"
                      ))
              print(c("unaccounted", Cyto_sample))
            
            }

          
          

      }
    }
    ##Deal with sex chromosomes here
      #remember to deal with 46,xxx
      
      count_before_extras = 0 ##theoretically original
      ploidy_count = ploidy ##count_before_mods * plody
      count_after_mods = 0 ## counts with deletions
      constitutionalcount=xconstitutional+yconstitutional      

      count_before_extras = xcount+ycount+xmod+ymod

      idealTotal=normX+normY
      sexDev_from_norm=0
      difference=0
      count_after_mods=count_before_extras
      
      ##if mitelman specifications are true and there is no sex count, assume constitutionality, just count -X and - Y straight
      if(forMtn==T & (xcount+ycount)==0){
        constitutional=T
        if((ydel) >= 1)
        {
          for(f in 1:ydel){
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Loss", nrow(ref_table)))
            ##temp_table <-
            temp_table <-
              temp_table[grep("chrY", temp_table[, 1]), ]
            
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            colnames(temp_table)<-colnames(sample_table)
            
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
          }
        }
        if((xdel) >= 1){
          for(f in 1:xdel){
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Loss", nrow(ref_table)))
            ##temp_table <-
            temp_table <-
              temp_table[grep("chrX", temp_table[, 1]), ]
            
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            colnames(temp_table)<-colnames(sample_table)
            
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
          }
        } 
      }else if(constitutionalcount>0 && constitutional==F){
        ##if the x count is consitutional and we dont want to count the constitutional state, just print -xs as is 
        
        if((ydel) >= 1)
        {
          for(f in 1:ydel){
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Loss", nrow(ref_table)))
            ##temp_table <-
            temp_table <-
              temp_table[grep("chrY", temp_table[, 1]), ]
            
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            colnames(temp_table)<-colnames(sample_table)
            
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
          }
        }
        if((xdel) >= 1){
          for(f in 1:xdel){
            temp_table <-
              data.frame(ref_table[, 1],
                         rep(0, nrow(ref_table)),
                         ref_table[, 2],
                         rep("Loss", nrow(ref_table)))
            ##temp_table <-
            temp_table <-
              temp_table[grep("chrX", temp_table[, 1]), ]
            
            temp_table[, 4] <- as.character(temp_table[, 4])
            temp_table <- as.matrix(temp_table)
            colnames(temp_table)<-colnames(sample_table)
            
            sample_table[, 4] <- as.character(sample_table[, 4])
            sample_table <- rbind(sample_table, temp_table)
            sample_table[, 4] <- as.character(sample_table[, 4])
          }
        } 
      }else{
        ##do the complicated calculations
        ##if we need to take into account ploidy
        ##handle 69,xx,-y (would be one xgain, one y loss)
        if((idealTotal* (ploidy-2)/2 + idealTotal) != (count_before_extras) |( (count_before_extras)+xdel+ydel > (idealTotal* (ploidy-2)/2 + idealTotal) & (xdel+ydel) > 0 )){
            if(ploidy>2)
            {
              idealTotal=idealTotal*(ploidy-2)/2+idealTotal
              
            }else if(ploidy==1){
                idealTotal=1
            }
  
            ##count_after_mods=(ploidy-2)*(count_before_extras) 
          
            
              if(constitutionalcount>0){
                
                count_after_mods=count_before_extras - xdel - ydel
              }else if(idealTotal == (count_before_extras+xdel+ydel)){
                
                ##think about this , (46,xx,-x,-x)
                count_after_mods=count_before_extras
              }else if(count_before_extras +xdel+ydel > idealTotal )
              {
                ##if theres overflow, determine if there is an addition and subtraction needed to be calculated, or if raw counts should count as total
                 if((count_before_extras)+xdel+ydel > idealTotal & (xdel+ydel) > 0){
                  
                   ##if the raw counts are greater than idealTotal, assume raw counts are the total and deletions are subtractions from that value
                   if(count_before_extras >= idealTotal){
                     
                     count_after_mods=count_before_extras
                     
                     ##set equal to ploidy, then subtract difference
                     count_after_mods=count_after_mods-xdel-ydel
                   }else{
                     ##calculate difference
                      count_after_mods=count_before_extras
                      
                      ##set equal to ploidy, then subtract difference
                      difference= xdel+ydel+ count_after_mods - ploidy
                      count_after_mods=count_after_mods-difference
                   }
                }else if(count_before_extras < idealTotal)
                {
                  count_after_mods=(ploidy-2)*(count_before_extras) 
                  
                  ##set equal to ploidy, then subtract difference
                  difference= xdel+ydel-count_after_mods
                  count_after_mods=count_after_mods+difference
                }else{
                  count_after_mods=count_before_extras+xdel+ydel
                  
                }
                
              }
            
            if(constitutional==F){
              ##add a gain/loss to counter to negate loss
              #if haploidy, add a gain to counteract
                if(constitutionalcount>2){
                  ynew=0
                  
                  if(yconstitutional>0){
                    for(f in 1:yconstitutional-1)
                    {
                      temp_table <-
                        data.frame(ref_table[, 1],
                                   rep(0, nrow(ref_table)),
                                   ref_table[, 2],
                                   rep("Loss", nrow(ref_table)))
                      ##temp_table <-
                      temp_table <-
                        temp_table[grep("chrY", temp_table[, 1]), ]
                      
                      temp_table[, 4] <- as.character(temp_table[, 4])
                      temp_table <- as.matrix(temp_table)
                      colnames(temp_table)<-colnames(sample_table)
                      
                      sample_table[, 4] <- as.character(sample_table[, 4])
                      sample_table <- rbind(sample_table, temp_table)
                      sample_table[, 4] <- as.character(sample_table[, 4])
                    }
                  }
                  
                  if(xconstitutional>0){
                    for(f in 1:(xconstitutional-xcount))
                    {
                      temp_table <-
                        data.frame(ref_table[, 1],
                                   rep(0, nrow(ref_table)),
                                   ref_table[, 2],
                                   rep("Loss", nrow(ref_table)))
                      temp_table <-
                        temp_table[grep("chrX", temp_table[, 1]), ]
                      
                      temp_table[, 4] <- as.character(temp_table[, 4])
                      temp_table <- as.matrix(temp_table)
                      colnames(temp_table)<-colnames(sample_table)
                      
                      sample_table[, 4] <- as.character(sample_table[, 4])
                      sample_table <- rbind(sample_table, temp_table)
                      sample_table[, 4] <- as.character(sample_table[, 4])
                      
                    }
                  }
                }
              
              }
              ##if polyploidy, loss ploidy-2
              ##for loop for each constitutional value
            }
          
            
          
              sexDev_from_norm <- count_after_mods-2
              
              if(sexDev_from_norm>0){
                  ##add gain of sex chrom for loop
                ##copy paste old code
                ynew=0
                for(f in 1:sexDev_from_norm)
                {
                  ##include if y and x discordanceb 
                  if((ycount+ymod) > (ynew +normY))
                  {
                    temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Gain", nrow(ref_table)))
                    ##temp_table <-
                    temp_table <-
                      temp_table[grep("chrY", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    ynew=ynew+1
                  }else{
                    
                    
                    temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Gain", nrow(ref_table)))
                    temp_table <-
                      temp_table[grep("chrX", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    
                  }
                }
              
            }else if(sexDev_from_norm < 0){
              ynew=0
              for(f in 1:(-1*sexDev_from_norm))
              {
               
            
                  if((ycount+ymod) < (ynew +normY)  )
                  {
                    temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Loss", nrow(ref_table)))
                    ##temp_table <-
                    temp_table <-
                      temp_table[grep("chrY", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    ynew=ynew+1
                    
                  }else{
                    
                    
                    temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Loss", nrow(ref_table)))
                    temp_table <-
                      temp_table[grep("chrX", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    
                  }
              }
            
            }else{
              ##if sexDev == 0, make sure to take into account -y and XX canceling out 
                if((ydel) >= 1)
                {
                    for(f in 1:ydel){
                     temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Loss", nrow(ref_table)))
                    ##temp_table <-
                    temp_table <-
                      temp_table[grep("chrY", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                  }
  
                  if(xcount+xmod > 1){
                    temp_table <-
                      data.frame(ref_table[, 1],
                                 rep(0, nrow(ref_table)),
                                 ref_table[, 2],
                                 rep("Gain", nrow(ref_table)))
                    ##temp_table <-
                    temp_table <-
                      temp_table[grep("chrX", temp_table[, 1]), ]
                    
                    temp_table[, 4] <- as.character(temp_table[, 4])
                    temp_table <- as.matrix(temp_table)
                    colnames(temp_table)<-colnames(sample_table)
                    
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    sample_table <- rbind(sample_table, temp_table)
                    sample_table[, 4] <- as.character(sample_table[, 4])
                    
                  }
                }
            
              }
         }
      
            ###count_after_mods -2 gives gains/losses
            ##make this more complicated to handle cases like 46,xxxc,-x (- x is printed)
            
          
     
        
  
      
      ##&& ((constitutionalcount==0) ||((constitutionalcount*(ploidy-2) + constitutionalcount) > (count_before_extras)))
        
      print(c("count_before_extras",count_before_extras,"ploidy_count", ploidy_count,"count_after_mods",count_after_mods ,
              "constitutionalcount",constitutionalcount,"idealtotal",idealTotal,"sexDev_from_norm",sexDev_from_norm,"difference",difference,"normx",normX,"normY",normY,"ycount",ycount,"ymod",ymod))
      
      
      ##conditional if not the same after sex counts 
      ##print(c("val_div>0", Cyto_sample))
      ##add this to uncertain 
      
      ##put in dump table
      ##Dump_table <-
     ##   rbind(Dump_table,
      ##        c(
     ##           as.vector(Con_data[i, ]),
      ##          "Warning in some chromosomes unaccounted for"
      ##        ))
     ## print(c("unaccounted", Cyto_sample))
      
      
    ##something is wrong here 
    if(any(is.na(sample_table[,2]))|any(is.na(sample_table[,3]))){
      sample_table<-sample_table[-union(which(is.na(sample_table[,2])) ,which(is.na(sample_table[,3]))),]
      ##Dump_table<-rbind(Dump_table,c(as.vector(Con_data[i, ]),
                                        ## "Error in NA found"))
     }    
        
    if(!is.vector(sample_table))
    {
      sample_table <- sample_table[rowSums(!is.na(sample_table)) > 0, ]
    }

        
        
        
      
    
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
      sorted_sample_table<-mergeTable_jp(sample_table)
    }
    
    ##sorted_sample_table<-sample_table
    
    ##correct format for sample table
    ##something is  wrong here
    if (is.vector(sorted_sample_table))
    {
      sorted_sample_table[1]<-as.character(sorted_sample_table[1])
      sorted_sample_table[2]<-as.integer(as.character(sorted_sample_table[2]))
      sorted_sample_table[3]<-as.integer(as.character(sorted_sample_table[3]))
      sorted_sample_table[4]<-as.character(sorted_sample_table[4])
      sorted_sample_table <- t(sorted_sample_table)
    }else if(ncol(sorted_sample_table)==1){
      sorted_sample_table <- t(sorted_sample_table)
      sorted_sample_table<-as.data.frame(sorted_sample_table)
      sorted_sample_table[,1]<-as.character(sorted_sample_table[,1])
      sorted_sample_table[,2]<-as.integer(as.character(sorted_sample_table[,2]))
      sorted_sample_table[,3]<-as.integer(as.character(sorted_sample_table[,3]))
      sorted_sample_table[,4]<-as.character(sorted_sample_table[,4])
    }else if(nrow(sorted_sample_table)>1){
      sorted_sample_table<-as.data.frame(sorted_sample_table)
      sorted_sample_table[,1]<-as.character(sorted_sample_table[,1])
      sorted_sample_table[,2]<-as.integer(as.character(sorted_sample_table[,2]))
      sorted_sample_table[,3]<-as.integer(as.character(sorted_sample_table[,3]))
      sorted_sample_table[,4]<-as.character(sorted_sample_table[,4])
    }
    
    return(list(sorted_sample_table,Dump_table,transloctable))
      
    }
  }
 

  
  ##try catch for each sample
  samparse <- function(i)
  {
    out <- try(rowparse(i))
    
    if (inherits(out, "try-error")) {
      return(gsub("\n"," ",paste(geterrmessage(), "in", i, "sample")))
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
    
    
    
    
    ##if idem or sl, cancel out any - details , add in any shorthands between code
    if(any(grepl("ids|idem|sl|sd",Cyto_sample)) || allow_Shorthand )
    {
      
      ##index of anything in a clonal evolution step with no )(
      ##think about what to do if assymmetric (2 defined, one mystery, one defined, 2 mystery)
      ##takes first occurance if more than 2 of the same category -- f
     ## fix it to skip
      mut_index<-grep("\\(.*\\)(?!\\()",Cyto_sample,perl=T)
      if(length(mut_index)>0)
      {
        mut_list<-sapply(Cyto_sample[mut_index],function(x){gsub("\\)","\\\\)",gsub("\\(","\\\\(",gsub("\\?","\\\\?",gsub("\\+","",substr(x,2,nchar(x))  ))))})
        index_match<-sapply(mut_list,function(x){grep(x,Cyto_sample)[1]})
        index_match<-index_match[intersect(which(!is.na(index_match)),which(index_match>0))]
        if(is.numeric(mut_index) && is.numeric(index_match)&& length(index_match)>0)
        {
          additional_to_add <- substring(names(index_match),0,regexpr("\\+",names(index_match)) )
          Cyto_sample[mut_index]<-paste(additional_to_add ,Cyto_sample[index_match],sep="")
        }
      }
      ##index of anything thats - in a clonal evolution step
      mut_gone_index<-grep("-[[:alpha:]]+\\(",Cyto_sample)
      if(length(mut_gone_index)>0)
      {
        mut_list<-sapply(Cyto_sample[mut_gone_index],function(x){gsub("\\)","\\\\)",gsub("\\(","\\\\(",gsub("\\?","\\\\?",substr(x,2,nchar(x)))))})
        index_cancel<-sapply(mut_list,function(x){grep(x,Cyto_sample)[1]})
        index_cancel<-index_cancel[intersect(which(!is.na(index_cancel)),which(index_cancel>0))]
        if(is.numeric(mut_gone_index) && is.numeric(index_cancel)&& length(index_cancel)>0)
        {
          Cyto_sample<-Cyto_sample[-1*c(mut_gone_index,index_cancel)]
        }
      }
    }
    
    ##check to make sure the input is roughly a karyotype before processing
    if(grepl("[[:digit:]]+((~|-)[[:digit:]]+)*(<[[:digit:]]+n>)*,",paste(Cyto_sample,collapse=',',sep=''))&&grepl("[[:digit:]]",Cyto_sample[1]))
    {
      
      
      tottable<-samparse(i)
      if(is.character(tottable) & length(tottable)==1)
      {
        Dump_table <- rbind(Dump_table, c(Con_data[i, ], tottable))
        transloctable<-data.frame()
      }else if(!is.list(tottable))
      {
        Dump_table <-  tottable
        transloctable<-data.frame()
      }
      else{
        
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
    }else{
      if(grepl("[[:digit:]]",Cyto_sample[1]))
      {
        Dump_table<-rbind(Dump_table,c(Con_data[i, ], "Warning in karyotype number not specified"))
        
      }else{
        Dump_table<-rbind(Dump_table,c(Con_data[i, ], "Warning in karyotype is incorrect"))
      }
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
  
  IDs <- unlist(lapply(IDs, function(x){
    paste(x, collapse = "_")
  }))
  IDs <- unique(IDs)
  
  for (i in 1:length(IDs)) {
    ##get all ids with same base name
    ##gsub away all regular expressions 
    w <- grep(gsub("\\+","\\\\+",gsub("\\*","\\\\*",gsub("\\!","\\\\!",paste("^", IDs[i], "_", sep = "")))), out[, 1])
    ##if (length(w) > 1) {
    nxt <- unlist(lapply(strsplit(as.character(out[w, 2]), split = "\\["),
                         function(x){
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
      if(any(is_cp))
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
  colnames(Dump_table) <- c("Sample ID","Karyotype","Error Message")
  
  result<-list(Final_Final_table,Dump_table)
  names(result)<-list("Result","Error_log")
  
return(result)
}



##file_list_CytoConverter<-function(in_data_vector,result_name,dir_name_file=NULL,dir_name_results=NULL,build="GRCh38",constitutional=T,guess=F,guess_q=F){
 
##  if(is.null(dir_name_results))
##  {
##    dir_name<-getwd()
##  }
##  if(is.null(dir_name_file))
##  {
##    dir_name<-getwd()
##  }
  
##for(i in 1:length(in_data_vector))
  ##{
    ##sapply(as.data.frame(read.delim(input$file$datapath,header=FALSE,sep="\t")),as.character)
    ##in_data_z<-as.matrix(sapply(as.data.frame(read.delim(file = paste(dir_name_file,"/",in_data_vector[1],sep=""),sep='\t',header=T)),as.character))
    ##results<-CytoConverter(in_data=in_data_z)
    ##write.table(file=paste(dir_name_results,"/",result_name,"_",in_data_vector[i],"_results",".txt",sep=""),x=results[[1]],sep='\t',row.names = F,quote=F)
    ##write.table(file=paste(dir_name_results,"/",result_name,"_",in_data_vector[i],"_error",".txt",sep=""),x=results[[2]],sep='\t',row.names = F,quote = F)
    
  ##}
  
  
##}



##
##file_list_CytoConverter(result_name = "cyto",in_data_vector = list.files("C:/Users/Janet/Downloads/CytoConverter/")[c(-1,-length(list.files("C:/Users/Janet/Downloads/CytoConverter/")))],dir_name_file = "C:/Users/Janet/Downloads/CytoConverter",dir_name_results ="C:/Users/Janet/Downloads/CytoConverter/Results" )

##files <- list.files(path = "C:/Users/Janet/Downloads/CytoConverter/Results/", pattern="*error*",full.names = T)

##tbl <- sapply(files, read.delim, simplify=FALSE) %>% 
##  bind_rows(.id = "id")
##tb<-tbl[grep("Error",tbl$X),]
##library(readr)
##library(dplyr)
##read.table()
##tbl <- sapply(files, read.delim, simplify=FALSE) %>% 
##  bind_rows(.id = "id")
##tb<-tbl[grep("Error",tbl$X),]





##else if((idealTotal* (ploidy-2) + idealTotal) < (count_before_extras)){
  ##make this more complex
##  count_after_mods=count_before_extras
  
