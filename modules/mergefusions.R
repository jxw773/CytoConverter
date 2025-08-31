#' Fusion-Specific Merge Functions
#' 
#' @description
#' This module contains specialized merge functions for handling fusion data and 
#' structural rearrangements. These functions extend the standard merge functionality
#' to properly handle fusion tags and complex chromosomal rearrangements.
#' 
#' @details
#' The fusion merge system maintains separate tracking for:
#' \itemize{
#'   \item Standard gains and losses
#'   \item Fusion events (marked with # tags)
#'   \item Plus-loss events (+Loss)
#'   \item Complex structural rearrangements
#' }
#' 
#' Key differences from standard merge functions:
#' \itemize{
#'   \item Preserves fusion tags during merge operations
#'   \item Handles overlapping fusion events appropriately
#'   \item Maintains breakpoint precision for structural aberrations
#'   \item Supports hierarchical fusion classification
#' }

##this is the the merge function for fusions 

#' Insert Fusion Section Function
#' 
#' @description
#' Inserts fusion regions into the data structure while maintaining proper tracking
#' of overlapping regions and fusion types. This function extends insertSection() 
#' to handle fusion-specific data.
#' 
#' @param h_ Hash data structure tracking non-overlapping sections within regions
#' @param start Starting coordinate of the fusion region
#' @param end Ending coordinate of the fusion region  
#' @param type Fusion type identifier (e.g., "#translocation", "#derivative_chrom")
#' 
#' @details
#' The function maintains the following structure in h_:
#' \itemize{
#'   \item Start coord -> End: end coordinate for this section
#'   \item Start coord -> Gain: list of end coords of "Gain" regions
#'   \item Start coord -> Loss: list of end coords of "Loss" regions  
#'   \item Start coord -> Fusion: list of end coords of fusion regions
#' }
#' 
#' For fusion types (those starting with "#"), the function:
#' \itemize{
#'   \item Stores fusion data in the "Fusion" category
#'   \item Preserves fusion type information for downstream processing
#'   \item Handles overlaps with existing gains/losses appropriately
#' }
#' 
#' @return Modified hash data structure with inserted fusion section
#' 
#' @examples
#' \dontrun{
#' # Insert a translocation region
#' h_ <- insertSection_Fus(h_, 1000000, 2000000, "#translocation_balanced")
#' 
#' # Insert a derivative chromosome region  
#' h_ <- insertSection_Fus(h_, 500000, 1500000, "#derivative_chrom")
#' }

insertSection_Fus <- function(h_, start, end, type) {
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
        Loss = h_[[h_key]][['Loss']],
        "+Loss" = h_[[h_key]][['+Loss']],
        Fusion=h_[[h_key]][[type]]
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
        Loss = h_[[h_key]][['Loss']],
        "+Loss" = h_[[h_key]][['+Loss']],
        Fusion=h_[[h_key]][[type]]
        
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
      if(grepl("#",type))
      {
        h_[[as.character(start_val)]][["Fusion"]][[
          length(h_[[as.character(start_val)]][["Fusion"]]) + 1
        ]] <- end  
      }else{
        h_[[as.character(start_val)]][[type]][[
          length(h_[[as.character(start_val)]][[type]]) + 1
        ]] <- end
      }
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
          Loss = list(),
          "+Loss" = list(),
          Fusion = list()
        )
        
        if(grepl("#",type)){
          h_[[as.character(prev_end)]][['Fusion']] <- list(end)
          
        }else{
          h_[[as.character(prev_end)]][[type]] <- list(end)
        }
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
      Loss = list(),
      "+Loss" = list(),
      Fusion = list()
    )
    
    if(grepl("#",type)){
      h_[[as.character(start)]][['Fusion']] <- list(end)
    }else{
      h_[[as.character(start)]][[type]] <- list(end)
    }
    
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
        Loss = list(),
        "+Loss" = list(),
        Fusion = list()
      )
      
      if(grepl("#",type)){
        h_[[as.character(prev_end)]][['Fusion']] <- list(end)
        
      }else{
        h_[[as.character(prev_end)]][[type]] <- list(end)
      }
    }
    
  }
  
  return(h_)
  
} # insertSection

deleteIntersections_Fus <- function(h_) {
  
  start_vals <- sort(as.numeric(hash::keys(h_)))
  
  for (start_val in start_vals) {
    ## If both gain and loss are > 0, then there is overlap
    ## If min_val is > 1, then there is duplicate overlap
    ## Only delete one gain for each loss or one loss for each gain
    fusion <- length(h_[[as.character(start_val)]][['Fusion']])
    loss <- length(h_[[as.character(start_val)]][['Loss']])
    ploss <- length(h_[[as.character(start_val)]][['+Loss']])
    
    if (fusion == (ploss + loss)) {
      # Delete section
      hash::delete(start_val, h_)
      
    } else {
      # Delete ploss first
      min_val <- min(fusion, ploss)
      
      # Delete leading elements from Gain and Loss lists
      if (min_val > 0) {
        for (i in 1:min_val) {
          h_[[as.character(start_val)]][['Fusion']][[1]] <- NULL
          h_[[as.character(start_val)]][['+Loss']][[1]] <- NULL
          
        }
        
      }
      
      ## If both gain and loss are > 0, then there is overlap
      ## If min_val is > 1, then there is duplicate overlap
      ## Only delete one gain for each loss or one loss for each gain
      gain <- length(h_[[as.character(start_val)]][['Fusion']])
      loss <- length(h_[[as.character(start_val)]][['Loss']])
      
      # Delete loss afterwards
      min_val <- min(fusion, loss)
      
      # Delete leading elements from Gain and pLoss lists
      if (min_val > 0) {
        for (i in 1:min_val) {
          h_[[as.character(start_val)]][['Fusion']][[1]] <- NULL
          h_[[as.character(start_val)]][['Loss']][[1]] <- NULL
          
        }
        
      }
      
    }
    
  }
  
  return(h_)
  
}

getContiguousSection_Fus <- function(h__, section, orig_end) {
  
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
        length(h__[[sect_end]][['Fusion']]) == 0
        && length(h__[[sect_end]][['Loss']]) == 0
        && length(h__[[sect_end]][['+Loss']]) == 0
      ) {
        # delete entire section if no more gain or loss
        hash::delete(sect_end, h__)
        
      }
      
      # continue searching recursively
      section <- getContiguousSection_Fus(h__, section, orig_end)
      
    }
    
  }
  
  return(section)
  
}

mergeAdjacentSections_Fus <- function(h_) {
  
  
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
      if (length(h_[[chr]][[ch_start]][['Fusion']]) > 0) {
        type <- 'Fusion'
        
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
        length(h_[[chr]][[ch_start]][['Fusion']]) == 0
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
      section <- getContiguousSection_Fus(h_[[chr]], section, orig_end)
      
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
      
      if (length(start_index) > 0 && type != "Fusion") {
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

mergeTable_Fus <- function(M, keep_extras = F) {
  
  # Store non gains and losses for intermediate steps
  M_temp <- M[grep("#|Loss", M[, 4], invert = T), ]
  
  # Create an empty hash map for each chromosome
  chrs <- unique(M[, 1])
  
  h <- hash::hash()
  for (chr in chrs) {
    h[chr] <- hash::hash()
  }
  
  # Populate the hash map with sections
  for (row in 1:nrow(M)) {
    if (
      grepl("#",M[row,'Type'] )
      || M[row, 'Type'] == 'Loss'
      || M[row, 'Type'] == '+Loss'
    ) {
      h[M[row, 'Chr']] <- insertSection_Fus(
        h[[M[row, 'Chr']]],
        as.numeric(M[row, 'Start']),
        as.numeric(M[row, 'End']),
        M[row, 'Type']
      )
      
    }
    
  }
  
  # Delete intersections
  for (key in hash::keys(h)) {
    h[key] <- deleteIntersections_Fus(h[[key]])
    
  }
  
  
  ##fix this to not merge fusions 
  # Merge adjacent sections and create final table
  final_coord_table <- mergeAdjacentSections_Fus(h)
  
  # If want other info
  if (keep_extras) {
    final_coord_table <- rbind(final_coord_table, M_temp)
    
  }
  
  return(final_coord_table)
  
}
