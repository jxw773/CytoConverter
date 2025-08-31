#' Cytoband Processing Module for CytoConverter
#' 
#' @description
#' This module contains functions for processing cytogenetic band information in the
#' context of chromosomal translocations, insertions, and complex structural
#' rearrangements. It handles the interpretation of cytogenetic nomenclature and
#' converts it to corresponding genomic regions while accounting for special cases
#' like isoderivative chromosomes and complex breakpoint descriptions.
#' 
#' @details
#' Key functionality includes:
#' \itemize{
#'   \item Processing breakpoint descriptions for translocations and insertions
#'   \item Handling complex cytogenetic notation including :: and -> separators
#'   \item Managing special chromosome features (centromeres, telomeres, heterochromatin)
#'   \item Supporting isoderivative chromosome analysis
#'   \item Integrating with Mitelman database conventions
#' }
#' 
#' The module works closely with the utils module for position sorting and
#' coordinate processing.

## function for getting cytobands for translocations and insertions
## (taking compliment of stuff not included in discription)

mod_utils <- modules::use('modules/utils.R')

#' Get Cytogenetic Bands for Structural Rearrangements
#' 
#' @description
#' This function processes cytogenetic band information for translocations and
#' insertions by identifying the complement of genomic regions not included in
#' the structural rearrangement description. It handles complex cytogenetic
#' nomenclature and converts it to genomic coordinates for downstream analysis.
#' 
#' @param Cyto_ref_table Reference table containing cytoband information with columns:
#'   chromosome, start, end, band, stain
#' @param Cyto_sample Vector containing the parsed components of the karyotype being analyzed
#' @param lengthcount Integer indicating the current position in the parsing sequence
#' @param o Integer index for accessing specific elements within the parsing structure
#' @param temp List containing parsed cytogenetic components and breakpoint information
#' @param coln Integer column index for accessing specific karyotype components
#' @param derMods Character vector containing derivative modification information
#' @param forMtn Logical flag indicating whether to apply Mitelman database-specific rules
#' 
#' @return List containing:
#'   \item{bands}{Character vector of cytogenetic bands representing regions involved}
#'   \item{earlyReturn}{Logical flag indicating whether processing should terminate early
#'     (used for Mitelman compatibility when multiple bands per chromosome are detected)}
#' 
#' @details
#' The function handles several complex scenarios:
#' 
#' **Isoderivative chromosomes**: When "ider" is detected, the function modifies
#' processing to account for the specific arm (p or q) that forms the isochromosome.
#' 
#' **Complex breakpoint notation**: Supports advanced cytogenetic notation including:
#' \itemize{
#'   \item :: separators for complex rearrangements
#'   \item -> and ~> for directional information
#'   \item Multiple chromosome involvement
#'   \item Terminal end specifications (pter, qter)
#'   \item Centromere references (cen, p10, q10)
#' }
#' 
#' **Mitelman compatibility**: When forMtn=TRUE, applies special rules for
#' translocations with multiple bands per chromosome to maintain database compatibility.
#' 
#' **Position processing**: Uses positionSorter from utils module to ensure
#' proper ordering of cytogenetic positions from p-arm to q-arm.
#' 
#' **Edge case handling**: Manages chromosome terminal positions and centromeric
#' regions with appropriate boundary conditions.
#' 
#' @examples
#' \dontrun{
#' # Process a simple translocation breakpoint
#' result <- getCytoBands(cyto_ref_table, cyto_sample, 1, 1, temp_data, 2, der_mods, FALSE)
#' 
#' # Handle complex breakpoint with multiple chromosomes
#' result <- getCytoBands(cyto_ref_table, cyto_sample, 2, 1, complex_temp, 3, der_mods, TRUE)
#' }
getCytoBands <- function(
        Cyto_ref_table,
        Cyto_sample,
        lengthcount,
        o,
        temp,
        coln,
        derMods,
        forMtn
) {
    # Must take into account acen and qter pter and only one listing
    # (will have to relte to two), reuse later code for this

    # Extract chromosome-specific cytoband data from reference table
    chr_table <- Cyto_ref_table[
        grep(
            paste(
                paste("chr", temp[[(lengthcount * 2 - 1)]][o], sep = ""),
                "$",
                sep = ""
            ),
            Cyto_ref_table
        ),
    ]
    
    # For isoderivative chromosomes, end point is potentially different,
    # make boolean now
    
    isiso = FALSE
    
    # Quit if karyotype returns false (used for Mitelman compatibility)
    earlyReturn = F
    
    # Detect isoderivative chromosomes and determine which arm is involved
    if (
        any(
            grepl("ider", Cyto_sample[coln])
            && grepl("t\\(", Cyto_sample[coln])
        )
    ) {

        isiso = TRUE
        # Extract the arm (p or q) that forms the isochromosome
        arm <- gsub("[[:digit:]]", "", temp[[(grep("ider", derMods)) + 1]])
        if (arm == "q") {
            unused = "p"  # If q arm isochromosome, p arm is unused

        }
        if (arm == "p") {
            unused = "q"  # If p arm isochromosome, q arm is unused

        }

    }
    
    currentvec <- vector()
    
    # Handle complex breakpoint notation with :: and -> separators
    if (any(grepl("::", temp[[lengthcount * 2]][o]) |
            grepl("~>", temp[[lengthcount * 2]][o]) |
            grepl("->", temp[[lengthcount * 2]][o]))) {
      # Parse complex format: chr1p11::chr2q22->chr3p13::chr4q24
      # Split on :: first, then on -> or ~> to get chromosome:position pairs
      # This format indicates complex rearrangements involving multiple breakpoints
      longform_table <-
        strsplit(strsplit(temp[[lengthcount * 2]][o], "::")[[1]], "(~>)|(->)")
      # Remove any leading colons from parsing artifacts
      longform_table <-
        lapply(longform_table, function(x) {
          gsub(':', '', x)
        })
      
      in_table = data.frame()
      
      # Mark for dicentric and tricentric chromosomes (add "long" prefix)
      if (grepl("dic|trc", derMods[lengthcount]))
      {
        addBool <- paste("long", addBool, sep = '')
      }
      
      # Process each segment of the complex breakpoint description
      for (j in 1:length(longform_table))
      {
        # Find positions of p and q arms in the description
        stringdx <- str_locate_all(pattern = "p|q", longform_table[[j]])
        # Extract chromosome name (everything before the first p or q)
        chr_name_long <-
          substr(longform_table[[j]][1], 0, stringdx[[1]][1] - 1)
        # Extract position information (from p/q to end of each segment)
        positions <-
          as.vector(cbind(
            substr(longform_table[[j]][1],
                   stringdx[[1]][1],
                   nchar(longform_table[[j]][1])),
            substr(longform_table[[j]][2],
                   stringdx[[2]][1],
                   nchar(longform_table[[j]][2]))
          ))
        
        # Use chromosome-specific reference table if chromosome is specified
        if (nchar(chr_name_long) != 0)
        {
          chr_table_2 <-
            Cyto_ref_table[grep(paste(paste("chr", chr_name_long, sep = ""), "$", sep = ""), Cyto_ref_table),]
        } else {
          chr_table_2 <- chr_table  # Use default chromosome table
        }
        
        # Handle special terminal and centromeric position references
        if (any(grepl("pter", positions)))
        {
          positions[grep("pter", positions)] <- chr_table_2[1, 4]  # p-terminal
        }
        if (any(grepl("qter", positions)))
        {
          positions[grep("qter", positions)] <-
            chr_table_2[nrow(chr_table_2), 4]  # q-terminal
        }
        # Handle centromere references (be careful about centromere boundaries)
        if (any(grepl("cen", positions)))
        {
          # Use first acrocentric band as centromere reference
          positions[grep("cen", positions)] <-
            chr_table_2[grep("acen", chr_table_2[, 5]),][1, 4]
        }
        
        # Handle p10/q10 centromere notation (standard cytogenetic convention)
        positions[grep("p10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]),][1, 4]  # p-side of centromere
        positions[grep("q10", positions)] <-
          chr_table_2[grep("acen", chr_table_2[, 5]),][2, 4]  # q-side of centromere
        
        # Ensure positions are properly ordered (p-arm to q-arm)
        positions <- mod_utils$positionSorter(positions)
        
        # Create position table by matching bands in reference table
        positions_table <-
          matrix(chr_table_2[grep(paste(positions, collapse = "|", sep = "|"),
                                  chr_table_2[, 4]),], ncol = 5)
        
        # Handle single row case (convert vector to matrix)
        if (is.vector(positions_table))
        {
          positions_table <- t(positions_table)
        }
        # Store start and end positions for this segment
        currentvec <-
          c(currentvec,
            as.vector(positions_table[, 4])[1],  # First position
            as.vector(positions_table[, 4])[length(as.vector(positions_table[, 4]))])  # Last position
        
      }
      
      
      
      
    } else {
      # Handle standard cytogenetic notation (simpler format)
      # Extract positions by splitting on p and q arms
      positions <-
        strsplit(gsub("q", ",q", gsub("p", ",p", temp[[lengthcount * 2]][o])), ",")[[1]][2:length(strsplit(gsub("q", ",q", gsub(
          "p", ",p", temp[[lengthcount * 2]][o]
        )), ",")[[1]])]
      # Handle centromere position notation
      positions[grep("p10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]),][1, 4]
      positions[grep("q10", positions)] <-
        chr_table[grep("acen", chr_table[, 5]),][2, 4]
      positions <- mod_utils$positionSorter(positions)
      
      
      
      
      ###########################################################################################################
      #####################for mitelman data only###############################################################
      #####check for more than one band per chromosome for translocations######################################
      ##############################################################################################################
      
      # Mitelman database compatibility: restrict complex translocations
      if (forMtn == T &
          grepl("t\\(", derMods[lengthcount]) &
          length(positions) > 1)
      {
        earlyReturn = T  # Signal early termination for Mitelman compatibility
        
      } else{
        ################################################
        ##############################################
        #########################################
        # Handle single position cases (add appropriate endpoints)
        if (length(positions) == 1)
        {
          if (any(grepl("p", positions)))
          {
            # For p-arm breakpoint, include region from breakpoint to q-terminal
            currentvec <-
              c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))] +
                                              1, 4], chr_table[nrow(chr_table), 4], sep = ''))
            ##currentvec <-
            ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))]
            ##           , 4], chr_table[nrow(chr_table), 4], sep = ''))                                - 1, 4], chr_table[1, 4], sep = ''))
          }
          if (any(grepl("q", positions)))
          {
            # For q-arm breakpoint, include region from p-terminal to breakpoint
            currentvec <-
              c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[1]
                         - 1, 4], chr_table[1, 4], sep = ''))
           ##currentvec <-
           ##c(currentvec, paste(chr_table[grep(positions, chr_table[, 4])[length(grep(positions, chr_table[, 4]))]
           ##           , 4], chr_table[nrow(chr_table), 4], sep = ''))
          }
        } else{
          # Handle multiple positions (complex breakpoints)
          # restrict if positions are at the ends of the chromosome
          pos <-
            grep(paste(positions, sep = "", collapse = "|"), chr_table[, 4])
          if (pos[length(pos)] + 1 > nrow(chr_table) &&
              pos[1] == 1)
          {
            # Breakpoints span entire chromosome
            currentvec <- c(currentvec,
                            paste(chr_table[nrow(chr_table), 4], sep = ''),
                            paste(chr_table[1, 4], sep = ''))
            
          } else if (pos[length(pos)] + 1 > nrow(chr_table)) {
            # Breakpoint extends to q-terminal
            currentvec <-   c(
              currentvec,
              paste(chr_table[nrow(chr_table), 4], sep = ''),
              paste(chr_table[pos[1] - 1, 4], chr_table[1, 4], sep = '')
            )
            
            
          } else if (pos[1] == 1) {
            # Breakpoint starts at p-terminal
            
            currentvec <-    c(
              currentvec,
              paste(chr_table[pos[length(pos)] +
                                1, 4], chr_table[nrow(chr_table), 4], sep = ''),
              paste(chr_table[1, 4], sep = '')
            )
            
            
          } else{
            # Standard case: breakpoints in middle of chromosome
            currentvec <-
              c(
                currentvec,
                paste(chr_table[pos[length(pos)] +
                                  1, 4], chr_table[nrow(chr_table), 4], sep = ''),
                paste(chr_table[pos[1] -
                                  1, 4], chr_table[1, 4], sep = '')
              )
          }
        }
      }
      
    }
    
    # Handle isoderivative chromosome adjustments
    if (isiso)
    {
      # Replace unused arm positions with the specified arm for isoderivatives
      currentvec <-
        sapply(currentvec, function(x) {
          gsub(paste(unused, "[[:digit:]]+(\\.[[:digit:]])*", sep = ''),
               temp[[(grep("ider", derMods)) + 1]],
               x)
        })
    }
    
    return(list(currentvec, earlyReturn))
  }


