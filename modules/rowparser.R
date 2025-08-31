#' CytoConverter Row Parser
#' 
#' @description
#' This function parses each row of a karyotype table, processing individual sample 
#' karyotypes into structured genomic coordinate data. It handles the complex logic 
#' for interpreting cytogenetic nomenclature and converting it to genomic intervals
#' representing gains, losses, and structural rearrangements.
#' 
#' @param cyto_ref_table Reference table containing cytoband information for coordinate mapping
#' @param ref_table Additional reference data for chromosome processing  
#' @param Cyto_sample Vector containing the parsed components of a single karyotype
#' @param Con_data Context data for the current sample being processed
#' @param transloctable Table for tracking translocation events and complex rearrangements
#' @param Dump_table Table for collecting error messages and warnings during parsing
#' @param constitutional Boolean flag indicating constitutional vs somatic analysis mode
#' @param guess Boolean flag to enable guessing of ambiguous chromosomal regions
#' @param guess_q Boolean flag for q-arm specific guessing logic
#' @param guess_by_first_val Boolean flag to guess coordinates based on first values
#' @param forMtn Boolean flag for Montreal nomenclature compatibility
#' @param orOption Boolean flag to enable OR logic in parsing ambiguous cases
#' @param sexstimate Boolean flag for sex chromosome estimation and normalization
#' @param count_fusions Boolean flag to enable fusion detection and specialized fusion processing
#' 
#' @return List containing:
#'   \item{sample_table}{Matrix with processed genomic intervals for standard aberrations}
#'   \item{sample_fusion_table}{Matrix with processed fusion events and structural rearrangements}
#'   \item{Dump_table}{Updated error/warning table with any parsing issues encountered}
#'   
#' @details
#' The row parser performs several key functions:
#' \itemize{
#'   \item Interprets cytogenetic nomenclature into genomic coordinates
#'   \item Distinguishes between constitutional and somatic karyotypes  
#'   \item Handles complex structural rearrangements when fusion counting is enabled
#'   \item Manages sex chromosome normalization and counting
#'   \item Processes deletion, duplication, and translocation events
#'   \item Applies various heuristics for ambiguous cases when guess flags are enabled
#' }
#' 
#' @note This function is called internally by the main CytoConverter function and typically
#' should not be called directly by end users. It expects pre-processed karyotype components
#' in the Cyto_sample vector.
#' 
#' @seealso \code{\link{colparse}} for column-level parsing logic

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
        guess_by_first_val,
        forMtn,
        orOption,
        sexstimate,
        count_fusions
) {

    sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sample_table) <- c("Chr", "Start", "End", "Type")

    sorted_sample_table <- matrix(ncol = 4, nrow = 0)
    colnames(sorted_sample_table) <- c("Chr", "Start", "End", "Type")
  
    
    sample_fusion_table <- matrix(ncol = 4, nrow = 0)
    colnames(sample_fusion_table) <- c("Chr", "Start", "End", "Type")
    
    sorted_sample_fusion_table <- matrix(ncol = 4, nrow = 0)
    colnames(sorted_sample_fusion_table) <- c("Chr", "Start", "End", "Type")
    
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
  
    xdel_Q<- 0 ##counts of question of mark xdel 
    ydel_Q <- 0 ##counts of question mark ydel
    
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



    if (length(Cyto_sample) > 0 &
        grepl("[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
              Con_data[2]))
    {
      ##initiate clonal evolution counter if ids is present
      if (grepl("^ids$", Cyto_sample[length(Cyto_sample)]))
      {
        ##dont count sex chromosomes here
        ##this has 46 chromosomes, as they get deleted or modified (without \\+), knocks em out of this tracker, stuff that remains is either gained or lost (likely gained)
        clone_chrom_tracker <- rep(1:22, 2)
      }
      
      if (grepl("[^[:alpha:]*][[:digit:]+],|[^[:alpha:]*][[:digit:]+$]",
                Con_data[2]))
      {
        ##set normal count for XY chromosomes
        ##think about what 46,X,+Y would mean
        if(grepl("(c$)|(c\\?$)",Cyto_sample[2]))
        {
          
          xconstitutional= stringr::str_count(Cyto_sample[2], "X") 
          yconstitutional= stringr::str_count(Cyto_sample[2], "Y") 
          
        }
        
        #if(constitutional==F & grepl("(c$)|(c\\?$)",Cyto_sample[2])){
        #  
        #  normX = stringr::str_count(Cyto_sample[2], "X")
        #  normY = stringr::str_count(Cyto_sample[2], "Y")
        
        #}else 
        if (any(grepl("Y", Cyto_sample)) &
            (sum(grepl("Y", Cyto_sample), na.rm = TRUE) > sum(grepl("\\+Y", Cyto_sample), na.rm =
                                                              TRUE)))
        {
          normX = 1
          normY = 1
        } else if(( (sexstimate==F ) && any(grepl("\\?",Cyto_sample[2])) )){
          ##if ? is a chromosome and sexstimate is off, dont make guesses on constitutional change
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
          length_XY <-length(which(grepl("Y|X", Cyto_sample)==T))
          
          if(sexstimate==F & guess_q==F & grepl("\\?",Cyto_sample[2]) & length_XY<3 ){
            
            if(any(grepl("Y",Cyto_sample[2])))
            {
              if(length_XY==1){
                ycount=1
              }
            }else{
              if(length_XY==0){
                xcount=2
                
              }else if(length_XY==1){
                xcount=1
              }
            }
            Dump_table <- rbind(Dump_table, c(Con_data[1], "Warning in ambiguous sex chromosome count. Default of XX set"))
            
          }
        }  else {
          xcount = stringr::str_count(Cyto_sample[2], "X")
          
          ##do not count ? for mitelman or no sex estimates
          if(sexstimate == T)
          {
            xcount = xcount + stringr::str_count(Cyto_sample[2], "\\?")
            
          }
          
          ycount = stringr::str_count(Cyto_sample[2], "Y")
          ##x,y calculation
          
          ##if sexstimate is false and guess_q is false and there is only one sex chromosome but several ? marks, set value =2
          if(sexstimate==F & guess_q==F & (xcount+ycount) < 2 & grepl("\\?",Cyto_sample[2]) & !any(grepl("Y",Cyto_sample[-2]))){
            if(ycount > 0)
            {
              xcount=1
              ycount=1
              Dump_table <- rbind(Dump_table, c(Con_data[1], "Warning in ambiguous sex chromosome count. Default of XY set"))
              
            }else{
              xcount=2
              Dump_table <- rbind(Dump_table, c(Con_data[1], "Warning in ambiguous sex chromosome count. Default of XX set"))
              
            }
          }
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
      
    
    
        for (j in startcol:length(Cyto_sample)) {
          
         miniverter_data<-tryCatch({
                  mod_colparser$miniverter(
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
                    normX ,  
                    normY ,   
                    xcount , 
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
                    
                    startcol,
                    count_fusions=count_fusions
                    
                  )
                }, error = function(e) {
                  return(gsub("\n", " ", paste(e, "in", j, "field")))
                }, finally = {
                  # print(paste("  Parsed field: ", Cyto_sample[j]))
                })
         
         if (is.character(miniverter_data)) {
           Dump_table <- rbind(Dump_table, c(Con_data[1], miniverter_data))
           
         }else if (is.null(miniverter_data)) {
           
         }else if (length(miniverter_data) == 1 && is.na(miniverter_data)) {
           Dump_table <- rbind(
             Dump_table,
             c(
               Con_data[1],
               "Error in more than one band associated with a chromosome in a translocation"
             )
           )
           
         }else{
           temp_table <- miniverter_data[[1]]
           Con_data <- miniverter_data[[2]]
           transloctable <- miniverter_data[[3]]
           Dump_table <- miniverter_data[[4]]
           normX <- miniverter_data[[5]]
           normY <- miniverter_data[[6]]
           xcount <- miniverter_data[[7]]
           ycount <- miniverter_data[[8]]
           xadd <- miniverter_data[[9]]
           yadd <- miniverter_data[[10]]
           xmod <- miniverter_data[[11]]
           ymod <- miniverter_data[[12]]
           xdel <- miniverter_data[[13]]
           ydel <- miniverter_data[[14]]
           xconstitutional<- miniverter_data[[15]]
           yconstitutional <- miniverter_data[[16]]
           idealx <- miniverter_data[[17]]
           idealy <- miniverter_data[[18]]
           addtot <- miniverter_data[[19]]
           deltot <- miniverter_data[[20]]
           modtot <- miniverter_data[[21]]
           n  <- miniverter_data[[22]]
           ploidy <- miniverter_data[[23]]
           startcol <- miniverter_data[[24]]
           temp_fusion_table<-miniverter_data[[25]]
           
          # combine running file and temp new file together
          if (length(temp_table) != 0) {
            temp_table <- apply(temp_table, 2, as.character)
            sample_table <- rbind(sample_table, temp_table)
            
          }
           
           if(count_fusions==T){
             if(length(temp_fusion_table) != 0){
               temp_fusion_table <- apply(temp_fusion_table, 2, as.character)
               sample_fusion_table <- rbind(sample_fusion_table, temp_fusion_table)
               
             }
           }
         }
         

        } # for (j in startcol:length(Cyto_sample))
        
    
    
      
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
          rbind(Dump_table, c(as.vector(Con_data[1]), "Error in unclear chrom number"))
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
                      as.vector(Con_data[1]),
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
                      as.vector(Con_data[1]),
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
                      as.vector(Con_data[1]),
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
        if((ydel-ydel_Q) >= 1)
        {
          for(f in 1:(ydel-ydel_Q)){
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
        if((xdel-xdel_Q) >= 1){
          for(f in 1:(xdel-xdel_Q)){
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
        
        if((ydel-ydel_Q) >= 1)
        {
          for(f in 1:(ydel-ydel_Q)){
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
        if((xdel-xdel_Q) >= 1){
          for(f in 1:(xdel-xdel_Q)){
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
          
        }else if((sexDev_from_norm ) < 0){
          ynew=0
          if((sexDev_from_norm + xdel_Q +ydel_Q) < 0)
          {
            for(f in 1:(-1*(sexDev_from_norm+xdel_Q+ydel_Q)))
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
          }
        }else{
          ##if sexDev == 0, make sure to take into account -y and XX canceling out 
          if((ydel-ydel_Q) >= 1)
          {
            for(f in 1:(ydel-ydel_Q)){
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
        sorted_sample_table<-mod_merge$mergeTable(sample_table)
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
        
        ##fusion sample table formatting
        # Something is wrong here
        if(count_fusions==T)
        {
              if (
                any(is.na(sample_fusion_table[, 2]))
                | any(is.na(sample_fusion_table[, 3]))
              ) {
                sample_fusion_table <- sample_fusion_table[
                  -union(which(is.na(sample_fusion_table[, 2])),
                         which(is.na(sample_fusion_table[, 3]))),
                ]
              }
              
              if (!is.vector(sample_fusion_table)) {
                sample_fusion_table <- sample_fusion_table[rowSums(!is.na(sample_fusion_table)) > 0, ]
              }
              
              sample_fusion_table <- as.data.frame(sample_fusion_table, row.names = FALSE)
              if (is.vector(sample_fusion_table)) {
                sample_fusion_table[1] <- as.character(sample_fusion_table[1])
                sample_fusion_table[2] <- as.integer(as.numeric(as.character(sample_fusion_table[2])))
                sample_fusion_table[3] <- as.integer(as.numeric(as.character(sample_fusion_table[3])))
                sample_fusion_table[4] <- as.character(sample_fusion_table[4])
                sample_fusion_table <- t(sample_fusion_table)
                sorted_sample_fusion_table <- sample_fusion_table
                
              } else if (ncol(sample_fusion_table) == 1) {
                sample_fusion_table <- t(sample_fusion_table)
                sample_fusion_table[, 1] <- as.character(sample_fusion_table[, 1])
                sample_fusion_table[, 2] <- as.integer(as.numeric(as.character(sample_fusion_table[, 2])))
                sample_fusion_table[, 3] <- as.integer(as.numeric(as.character(sample_fusion_table[, 3])))
                sample_fusion_table[, 4] <- as.character(sample_fusion_table[, 4])
                sorted_sample_fusion_table <- sample_fusion_table
                
              } else if (nrow(sample_fusion_table) > 1) {
                sample_fusion_table[, 1] <- as.character(sample_fusion_table[, 1])
                sample_fusion_table[, 2] <- as.integer(as.numeric(as.character(sample_fusion_table[, 2])))
                sample_fusion_table[, 3] <- as.integer(as.numeric(as.character(sample_fusion_table[, 3])))
                sample_fusion_table[, 4] <- as.character(sample_fusion_table[, 4])
                # Eliminate duplicates and gain/loss with same coordinates
                ##sorted_sample_fusion_table <- mod_merge$mergeTable(sample_fusion_table)
                sorted_sample_fusion_table <- sample_fusion_table
              }
              
              # Correct format for sample table
              # Something is  wrong here
              if (is.vector(sorted_sample_fusion_table)) {
                sorted_sample_fusion_table[1] <- as.character(sorted_sample_fusion_table[1])
                sorted_sample_fusion_table[2] <- as.integer(as.character(sorted_sample_fusion_table[2]))
                sorted_sample_fusion_table[3] <- as.integer(as.character(sorted_sample_fusion_table[3]))
                sorted_sample_fusion_table[4] <- as.character(sorted_sample_fusion_table[4])
                sorted_sample_fusion_table <- t(sorted_sample_fusion_table)
                
              } else if (ncol(sorted_sample_fusion_table) == 1) {
                sorted_sample_fusion_table <- t(sorted_sample_fusion_table)
                sorted_sample_fusion_table <- as.data.frame(sorted_sample_fusion_table)
                sorted_sample_fusion_table[, 1] <- as.character(sorted_sample_fusion_table[, 1])
                sorted_sample_fusion_table[, 2] <- as.integer(as.character(sorted_sample_fusion_table[, 2]))
                sorted_sample_fusion_table[, 3] <- as.integer(as.character(sorted_sample_fusion_table[, 3]))
                sorted_sample_fusion_table[, 4] <- as.character(sorted_sample_fusion_table[, 4])
                
              } else if (nrow(sorted_sample_fusion_table) > 1) {
                sorted_sample_fusion_table <- as.data.frame(sorted_sample_fusion_table)
                sorted_sample_fusion_table[, 1] <- as.character(sorted_sample_fusion_table[, 1])
                sorted_sample_fusion_table[, 2] <- as.integer(as.character(sorted_sample_fusion_table[, 2]))
                sorted_sample_fusion_table[, 3] <- as.integer(as.character(sorted_sample_fusion_table[, 3]))
                sorted_sample_fusion_table[, 4] <- as.character(sorted_sample_fusion_table[, 4])
                
              }
        }else{
          sorted_sample_fusion_table<-NULL
        }
        
        
        if(count_fusions==F)
        {
          return(list(sorted_sample_table, Dump_table, transloctable,NULL))
        }else{
          return(list(sorted_sample_table, Dump_table, transloctable,sorted_sample_fusion_table))
          
        }
    }

}

