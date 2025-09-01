#' Fusion Color Assignment Module
#' 
#' @description
#' This module provides functions for automatic color assignment to fusion types
#' in CytoConverter visualizations. It extracts unique fusion types from data,
#' generates colorblind-friendly color palettes, and ensures consistent color
#' mapping across all visualizations.
#' 
#' @details
#' The module implements:
#' \itemize{
#'   \item Fusion type extraction (excluding |chrom portions)
#'   \item Colorblind-friendly palette generation
#'   \item Consistent color mapping for fusion types
#'   \item Legend support for color-coded fusion types
#' }

#' Extract Unique Fusion Types
#' 
#' @description
#' Extracts unique fusion types from a data frame containing fusion information.
#' Removes the |chrom portion from fusion tags to identify distinct fusion types.
#' 
#' @param fusion_data Data frame with fusion information containing a Type column
#' @return Vector of unique fusion types (without |chrom portions)
#' 
#' @details
#' Fusion tags typically follow patterns like:
#' \itemize{
#'   \item "#translocation_balanced|chrom_1" → "translocation_balanced"
#'   \item "#derivative_chrom::translocation" → "derivative_chrom::translocation"
#'   \item "#insertion_chrom::inserted_piece|chrom_1" → "insertion_chrom::inserted_piece"
#' }
#' 
#' @examples
#' \dontrun{
#' fusion_data <- data.frame(
#'   Type = c("#translocation_balanced|chrom_1", 
#'            "#translocation_balanced|chrom_2",
#'            "#derivative_chrom::translocation")
#' )
#' types <- extract_fusion_types(fusion_data)
#' # Returns: c("translocation_balanced", "derivative_chrom::translocation")
#' }
extract_fusion_types <- function(fusion_data) {
  if (is.null(fusion_data) || nrow(fusion_data) == 0) {
    return(character(0))
  }
  
  # Get all fusion tags (those starting with #)
  fusion_tags <- fusion_data[grepl("^#", fusion_data[, 5]), 5]
  
  if (length(fusion_tags) == 0) {
    return(character(0))
  }
  
  # Remove the # prefix and any |chrom_X portion
  cleaned_types <- gsub("^#", "", fusion_tags)
  cleaned_types <- gsub("\\|chrom_.*$", "", cleaned_types)
  
  # Return unique fusion types
  unique(cleaned_types)
}

#' Generate Colorblind-Friendly Color Palette
#' 
#' @description
#' Generates a colorblind-friendly color palette for the specified number of colors.
#' Uses scientifically validated color schemes that are distinguishable for most
#' types of color vision deficiencies.
#' 
#' @param n Number of colors needed
#' @return Vector of color codes (hex format)
#' 
#' @details
#' The function uses multiple strategies:
#' \itemize{
#'   \item For n <= 12: Uses a curated set of colorblind-friendly colors
#'   \item For n > 12: Supplements with additional distinguishable colors
#'   \item Colors are optimized for both normal and colorblind vision
#' }
#' 
#' Color palette based on Paul Tol's schemes and colorbrewer2.org recommendations.
#' 
#' @examples
#' \dontrun{
#' colors <- generate_colorblind_palette(5)
#' # Returns 5 distinct, colorblind-friendly colors
#' }
generate_colorblind_palette <- function(n) {
  if (n <= 0) {
    return(character(0))
  }
  
  # Base colorblind-friendly palette (Paul Tol's bright scheme + additions)
  base_colors <- c(
    "#4477AA", # Blue
    "#EE6677", # Red
    "#228833", # Green
    "#CCBB44", # Yellow
    "#66CCEE", # Cyan
    "#AA3377", # Purple
    "#BBBBBB", # Grey
    "#EE7733", # Orange
    "#009988", # Teal
    "#CC6677", # Pink
    "#DDCC77", # Olive
    "#117733"  # Dark Green
  )
  
  if (n <= length(base_colors)) {
    return(base_colors[1:n])
  }
  
  # For more than base colors, generate additional distinguishable colors
  additional_colors <- c(
    "#882255", # Wine
    "#44AA99", # Turquoise
    "#999933", # Olive drab
    "#332288", # Indigo
    "#88CCEE", # Light blue
    "#DDDDDD", # Light grey
    "#661100", # Dark red
    "#AA4499", # Magenta
    "#44BB99", # Sea green
    "#EEBB77"  # Tan
  )
  
  all_colors <- c(base_colors, additional_colors)
  
  if (n <= length(all_colors)) {
    return(all_colors[1:n])
  }
  
  # If we need even more colors, use interpolation
  # This is a fallback for very large numbers of fusion types
  extra_needed <- n - length(all_colors)
  interpolated <- colorRampPalette(all_colors)(extra_needed)
  
  return(c(all_colors, interpolated))
}

#' Create Fusion Type Color Mapping
#' 
#' @description
#' Creates a consistent mapping between fusion types and colors. This function
#' ensures that the same fusion type always gets the same color across different
#' visualizations and sessions.
#' 
#' @param fusion_types Vector of unique fusion types
#' @return Named vector where names are fusion types and values are color codes
#' 
#' @details
#' The mapping is deterministic - it will always assign the same color to the
#' same fusion type, ensuring consistency across:
#' \itemize{
#'   \item Multiple plots in the same session
#'   \item Different sessions with the same data
#'   \item Text input vs file input processing
#' }
#' 
#' @examples
#' \dontrun{
#' types <- c("translocation_balanced", "derivative_chrom::translocation")
#' mapping <- create_fusion_color_mapping(types)
#' # Returns: named vector with consistent color assignments
#' }
create_fusion_color_mapping <- function(fusion_types) {
  if (length(fusion_types) == 0) {
    return(character(0))
  }
  
  # Sort fusion types for consistent ordering
  sorted_types <- sort(fusion_types)
  
  # Generate colors
  colors <- generate_colorblind_palette(length(sorted_types))
  
  # Create named vector
  names(colors) <- sorted_types
  
  return(colors)
}

#' Assign Colors to Fusion Data
#' 
#' @description
#' Assigns colors to fusion data based on fusion types. This function adds
#' a color column to the fusion data frame for use in plotting functions.
#' 
#' @param fusion_data Data frame with fusion information
#' @param color_mapping Named vector of fusion type to color mappings
#' @return Data frame with added Color column
#' 
#' @details
#' The function:
#' \itemize{
#'   \item Extracts fusion type from each row's Type column
#'   \item Looks up the corresponding color from the mapping
#'   \item Adds a Color column to the data frame
#'   \item Handles non-fusion entries (assigns default color)
#' }
#' 
#' @examples
#' \dontrun{
#' fusion_data <- data.frame(
#'   Chr = c("chr9", "chr22"),
#'   Start = c(1000, 2000),
#'   End = c(1500, 2500),
#'   Type = c("#translocation_balanced|chrom_1", "#translocation_balanced|chrom_2")
#' )
#' color_mapping <- create_fusion_color_mapping(c("translocation_balanced"))
#' colored_data <- assign_fusion_colors(fusion_data, color_mapping)
#' }
assign_fusion_colors <- function(fusion_data, color_mapping) {
  if (is.null(fusion_data) || nrow(fusion_data) == 0) {
    return(fusion_data)
  }
  
  # Initialize color column with default color for non-fusion entries
  fusion_data$Color <- "#CCCCCC"  # Light grey for non-fusion
  
  # Find fusion entries
  fusion_rows <- grepl("^#", fusion_data[, 5])
  
  if (any(fusion_rows)) {
    # Extract fusion types for fusion rows
    fusion_tags <- fusion_data[fusion_rows, 5]
    cleaned_types <- gsub("^#", "", fusion_tags)
    cleaned_types <- gsub("\\|chrom_.*$", "", cleaned_types)
    
    # Assign colors based on mapping
    for (i in which(fusion_rows)) {
      fusion_type <- cleaned_types[sum(fusion_rows[1:i])]
      if (fusion_type %in% names(color_mapping)) {
        fusion_data$Color[i] <- color_mapping[fusion_type]
      }
    }
  }
  
  return(fusion_data)
}

#' Generate Fusion Type Legend Data
#' 
#' @description
#' Generates data for creating a legend that shows fusion types and their
#' corresponding colors. This is used for creating informative legends
#' in fusion visualizations.
#' 
#' @param color_mapping Named vector of fusion type to color mappings
#' @return Data frame with fusion types and colors for legend creation
#' 
#' @details
#' The returned data frame contains:
#' \itemize{
#'   \item fusion_type: Clean, human-readable fusion type names
#'   \item color: Corresponding color codes
#'   \item display_name: Formatted names for display in legend
#' }
#' 
#' @examples
#' \dontrun{
#' color_mapping <- create_fusion_color_mapping(c("translocation_balanced"))
#' legend_data <- generate_fusion_legend_data(color_mapping)
#' }
generate_fusion_legend_data <- function(color_mapping) {
  if (length(color_mapping) == 0) {
    return(data.frame(
      fusion_type = character(0),
      color = character(0),
      display_name = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create display names (convert underscores to spaces, capitalize)
  display_names <- names(color_mapping)
  display_names <- gsub("_", " ", display_names)
  display_names <- gsub("::", " - ", display_names)
  
  # Capitalize first letter of each word
  display_names <- sapply(strsplit(display_names, " "), function(x) {
    paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
  })
  
  return(data.frame(
    fusion_type = names(color_mapping),
    color = as.character(color_mapping),
    display_name = display_names,
    stringsAsFactors = FALSE
  ))
}