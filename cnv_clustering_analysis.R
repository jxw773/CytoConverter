#' CNV Hierarchical Clustering Analysis Script
#' 
#' This script performs hierarchical clustering analysis on CNV (Copy Number Variation) data
#' based on Gain/Loss patterns. It reads tab-delimited CNV data, converts it to a binary
#' matrix, performs hierarchical clustering, and generates visualizations.
#' 
#' Author: CytoConverter Team
#' Date: 2024

# Load required libraries
if(!require("cluster", quietly = TRUE)){
  stop("cluster package is required but not installed")
}

if(!require("graphics", quietly = TRUE)){
  stop("graphics package is required but not installed")  
}

if(!require("grDevices", quietly = TRUE)){
  stop("grDevices package is required but not installed")
}

if(!require("stats", quietly = TRUE)){
  stop("stats package is required but not installed")
}

#' Read and process CNV data from tab-delimited file
#' 
#' @param file_path Path to the tab-delimited CNV data file
#' @return Data frame with CNV data
#' @examples
#' cnv_data <- read_cnv_data("cyto_result.txt")
read_cnv_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Read the tab-delimited file
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Check if required columns exist
  required_cols <- c("Sample.ID", "Chr", "Start", "End", "Type")
  if (!all(required_cols %in% colnames(data))) {
    stop(paste("Required columns missing. Expected:", paste(required_cols, collapse = ", ")))
  }
  
  # Filter for Gain and Loss only, convert unknown to 0 later
  cat("Read", nrow(data), "CNV records from", file_path, "\n")
  
  return(data)
}

#' Create genomic bins for CNV analysis
#' 
#' @param cnv_data CNV data frame
#' @param bin_size Size of genomic bins in base pairs (default: 10MB)
#' @return Data frame with genomic bins
create_genomic_bins <- function(cnv_data, bin_size = 10000000) {
  # Get unique chromosomes
  chromosomes <- unique(cnv_data$Chr)
  chromosomes <- chromosomes[order(chromosomes)]
  
  # Create bins for each chromosome
  bins <- data.frame()
  
  for (chr in chromosomes) {
    chr_data <- cnv_data[cnv_data$Chr == chr, ]
    if (nrow(chr_data) == 0) next
    
    # Get chromosome range
    chr_start <- min(chr_data$Start, na.rm = TRUE)
    chr_end <- max(chr_data$End, na.rm = TRUE)
    
    # Create bins
    bin_starts <- seq(chr_start, chr_end, by = bin_size)
    bin_ends <- c(bin_starts[-1] - 1, chr_end)
    
    chr_bins <- data.frame(
      Chr = chr,
      Start = bin_starts,
      End = bin_ends,
      Bin_ID = paste(chr, bin_starts, bin_ends, sep = "_")
    )
    
    bins <- rbind(bins, chr_bins)
  }
  
  return(bins)
}

#' Convert CNV data to binary matrix
#' 
#' @param cnv_data CNV data frame
#' @param use_bins Whether to use genomic bins (TRUE) or actual CNV regions (FALSE)
#' @param bin_size Size of genomic bins if use_bins = TRUE
#' @return Binary matrix with samples as rows and genomic regions as columns
#' Values: 1 = Gain, -1 = Loss, 0 = unknown/no data
cnv_to_binary_matrix <- function(cnv_data, use_bins = FALSE, bin_size = 10000000) {
  
  # Get unique samples
  samples <- unique(cnv_data$Sample.ID)
  samples <- samples[order(samples)]
  
  if (use_bins) {
    # Create genomic bins
    bins <- create_genomic_bins(cnv_data, bin_size)
    regions <- paste(bins$Chr, bins$Start, bins$End, sep = "_")
    
    # Initialize binary matrix
    binary_matrix <- matrix(0, nrow = length(samples), ncol = nrow(bins))
    rownames(binary_matrix) <- samples
    colnames(binary_matrix) <- regions
    
    # Fill matrix based on overlap with bins
    for (i in 1:nrow(cnv_data)) {
      sample <- cnv_data$Sample.ID[i]
      chr <- cnv_data$Chr[i]
      start <- cnv_data$Start[i]
      end <- cnv_data$End[i]
      type <- cnv_data$Type[i]
      
      # Find overlapping bins
      chr_bins <- bins[bins$Chr == chr, ]
      overlapping <- which(chr_bins$Start <= end & chr_bins$End >= start)
      
      if (length(overlapping) > 0) {
        value <- ifelse(type == "Gain", 1, ifelse(type == "Loss", -1, 0))
        for (bin_idx in overlapping) {
          global_idx <- which(bins$Chr == chr)[bin_idx]
          binary_matrix[sample, global_idx] <- value
        }
      }
    }
    
  } else {
    # Use actual CNV regions
    regions <- paste(cnv_data$Chr, cnv_data$Start, cnv_data$End, sep = "_")
    unique_regions <- unique(regions)
    unique_regions <- unique_regions[order(unique_regions)]
    
    # Initialize binary matrix
    binary_matrix <- matrix(0, nrow = length(samples), ncol = length(unique_regions))
    rownames(binary_matrix) <- samples
    colnames(binary_matrix) <- unique_regions
    
    # Fill matrix
    for (i in 1:nrow(cnv_data)) {
      sample <- cnv_data$Sample.ID[i]
      region <- regions[i]
      type <- cnv_data$Type[i]
      
      value <- ifelse(type == "Gain", 1, ifelse(type == "Loss", -1, 0))
      binary_matrix[sample, region] <- value
    }
  }
  
  cat("Created binary matrix with", nrow(binary_matrix), "samples and", ncol(binary_matrix), "regions\n")
  return(binary_matrix)
}

#' Perform hierarchical clustering on CNV binary matrix
#' 
#' @param binary_matrix Binary matrix from cnv_to_binary_matrix
#' @param distance_method Distance method for clustering (default: "binary")
#' @param clustering_method Clustering method (default: "complete")
#' @return List containing distance matrix and hierarchical clustering object
perform_hierarchical_clustering <- function(binary_matrix, distance_method = "binary", clustering_method = "complete") {
  
  # Remove columns (regions) with all zeros
  non_zero_cols <- apply(binary_matrix, 2, function(x) any(x != 0))
  if (sum(non_zero_cols) == 0) {
    stop("No CNV regions found in the data")
  }
  
  filtered_matrix <- binary_matrix[, non_zero_cols, drop = FALSE]
  cat("Using", ncol(filtered_matrix), "non-zero CNV regions for clustering\n")
  
  # Calculate distance matrix
  if (distance_method == "binary") {
    # Binary distance (Jaccard distance)
    dist_matrix <- dist(filtered_matrix, method = "binary")
  } else {
    # Other distance methods
    dist_matrix <- dist(filtered_matrix, method = distance_method)
  }
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = clustering_method)
  
  cat("Performed hierarchical clustering using", distance_method, "distance and", clustering_method, "linkage\n")
  
  return(list(
    distance_matrix = dist_matrix,
    clustering = hc,
    binary_matrix = filtered_matrix
  ))
}

#' Create dendrogram visualization
#' 
#' @param clustering_result Result from perform_hierarchical_clustering
#' @param output_file Output file path for saving the plot (optional)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param main Plot title
create_dendrogram <- function(clustering_result, output_file = NULL, width = 10, height = 8, main = "CNV Hierarchical Clustering Dendrogram") {
  
  if (!is.null(output_file)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  }
  
  # Create dendrogram
  plot(clustering_result$clustering, 
       main = main,
       xlab = "Samples",
       ylab = "Distance",
       cex = 0.8,
       cex.main = 1.2)
  
  # Add rectangles to highlight clusters if there are enough samples
  if (length(clustering_result$clustering$labels) >= 4) {
    # Determine optimal number of clusters (simple heuristic)
    n_clusters <- min(4, ceiling(length(clustering_result$clustering$labels) / 2))
    rect.hclust(clustering_result$clustering, k = n_clusters, border = "red")
  }
  
  if (!is.null(output_file)) {
    dev.off()
    cat("Dendrogram saved to:", output_file, "\n")
  }
}

#' Create heatmap of CNV patterns
#' 
#' @param clustering_result Result from perform_hierarchical_clustering
#' @param output_file Output file path for saving the plot (optional)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param main Plot title
create_cnv_heatmap <- function(clustering_result, output_file = NULL, width = 12, height = 8, main = "CNV Gain/Loss Patterns Heatmap") {
  
  if (!is.null(output_file)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  }
  
  # Get the binary matrix and reorder according to clustering
  binary_matrix <- clustering_result$binary_matrix
  sample_order <- clustering_result$clustering$order
  ordered_matrix <- binary_matrix[sample_order, , drop = FALSE]
  
  # Create color palette: blue for loss (-1), white for no change (0), red for gain (1)
  colors <- c("blue", "white", "red")
  breaks <- c(-1.5, -0.5, 0.5, 1.5)
  
  # Set up plot
  par(mar = c(5, 8, 4, 2) + 0.1)
  
  # Create heatmap
  image(1:ncol(ordered_matrix), 1:nrow(ordered_matrix), t(ordered_matrix),
        col = colors, breaks = breaks,
        xlab = "Genomic Regions", ylab = "Samples",
        main = main, axes = FALSE)
  
  # Add sample labels on y-axis
  if (nrow(ordered_matrix) <= 50) {  # Only show labels if not too many samples
    axis(2, at = 1:nrow(ordered_matrix), labels = rownames(ordered_matrix), las = 2, cex.axis = 0.7)
  } else {
    axis(2, at = seq(1, nrow(ordered_matrix), by = 5), labels = seq(1, nrow(ordered_matrix), by = 5))
  }
  
  # Add x-axis
  axis(1, at = seq(1, ncol(ordered_matrix), by = max(1, ncol(ordered_matrix) %/% 10)), 
       labels = seq(1, ncol(ordered_matrix), by = max(1, ncol(ordered_matrix) %/% 10)))
  
  # Add legend
  legend("topright", legend = c("Loss", "No Change", "Gain"), 
         fill = colors, cex = 0.8, bg = "white")
  
  if (!is.null(output_file)) {
    dev.off()
    cat("Heatmap saved to:", output_file, "\n")
  }
}

#' Main function to perform complete CNV clustering analysis
#' 
#' @param input_file Path to input CNV data file
#' @param output_dir Output directory for results (default: current directory)
#' @param use_bins Whether to use genomic bins instead of actual CNV regions
#' @param bin_size Size of genomic bins if use_bins = TRUE
#' @param distance_method Distance method for clustering
#' @param clustering_method Clustering method
#' @return List containing all analysis results
cnv_clustering_analysis <- function(input_file, output_dir = ".", use_bins = FALSE, bin_size = 10000000,
                                  distance_method = "binary", clustering_method = "complete") {
  
  cat("Starting CNV clustering analysis...\n")
  cat("Input file:", input_file, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read and process data
  cnv_data <- read_cnv_data(input_file)
  
  # Convert to binary matrix
  binary_matrix <- cnv_to_binary_matrix(cnv_data, use_bins = use_bins, bin_size = bin_size)
  
  # Perform clustering
  clustering_result <- perform_hierarchical_clustering(binary_matrix, distance_method, clustering_method)
  
  # Generate visualizations
  dendrogram_file <- file.path(output_dir, "cnv_dendrogram.png")
  heatmap_file <- file.path(output_dir, "cnv_heatmap.png")
  
  create_dendrogram(clustering_result, dendrogram_file)
  create_cnv_heatmap(clustering_result, heatmap_file)
  
  # Save binary matrix and clustering results
  matrix_file <- file.path(output_dir, "cnv_binary_matrix.txt")
  write.table(clustering_result$binary_matrix, matrix_file, sep = "\t", quote = FALSE)
  cat("Binary matrix saved to:", matrix_file, "\n")
  
  # Save cluster assignments
  if (nrow(clustering_result$binary_matrix) >= 2) {
    n_clusters <- min(4, ceiling(nrow(clustering_result$binary_matrix) / 2))
    clusters <- cutree(clustering_result$clustering, k = n_clusters)
    cluster_file <- file.path(output_dir, "cnv_cluster_assignments.txt")
    write.table(data.frame(Sample = names(clusters), Cluster = clusters), 
                cluster_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Cluster assignments saved to:", cluster_file, "\n")
  }
  
  cat("CNV clustering analysis completed successfully!\n")
  cat("Results saved in:", output_dir, "\n")
  
  return(list(
    cnv_data = cnv_data,
    binary_matrix = clustering_result$binary_matrix,
    distance_matrix = clustering_result$distance_matrix,
    clustering = clustering_result$clustering
  ))
}

# Example usage (commented out)
# # Run analysis on example data
# if (file.exists("cyto_result.txt")) {
#   results <- cnv_clustering_analysis("cyto_result.txt", output_dir = "clustering_results")
# }