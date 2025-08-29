#' CNV Annotation Functions for CytoConverter
#' 
#' This script provides functionality to annotate CNV regions with gene information
#' and perform pathway analysis on the identified genes.
#' 
#' @author CytoConverter Development Team
#' @version 1.0

# Required libraries
if(!require(stringr, quietly = TRUE)) {
  stop("stringr package is required but not installed")
}
if(!require(dplyr, quietly = TRUE)) {
  stop("dplyr package is required but not installed")
}

#' Load Gene Annotation Data
#' 
#' This function loads gene annotation data for the specified genome build.
#' Note: In a real implementation, this would load actual gene annotation files
#' such as GTF/GFF files or data from Ensembl/UCSC databases.
#' 
#' @param build Genome build (GRCh38, hg19, hg18, hg17)
#' @return A data frame with gene annotations
#' @export
load_gene_annotations <- function(build = "GRCh38") {
  
  # Validate input
  if (!build %in% c("GRCh38", "hg19", "hg18", "hg17")) {
    stop("Invalid build. Must be one of: GRCh38, hg19, hg18, hg17")
  }
  
  # For demonstration purposes, create a mock gene annotation dataset
  # In a real implementation, this would load from actual gene annotation files
  sample_genes <- data.frame(
    gene_symbol = c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "APC", "PTEN", 
                    "RB1", "CDKN2A", "ATM", "ERBB2", "MLH1", "MSH2", "VHL", "NF1",
                    "BCL2", "CCND1", "MDM2", "FGFR2", "NOTCH1", "WNT1", "CTNNB1", "AKT1"),
    chromosome = c("chr17", "chr17", "chr7", "chr8", "chr12", "chr3", "chr5", "chr10",
                   "chr13", "chr9", "chr11", "chr17", "chr3", "chr2", "chr3", "chr17",
                   "chr18", "chr11", "chr12", "chr10", "chr9", "chr12", "chr3", "chr14"),
    start_pos = c(7565097, 43044295, 55019017, 127735434, 25205246, 179198076, 112043414, 87863113,
                  48365646, 21967751, 108098351, 39687914, 36993332, 47630108, 10146555, 31094927,
                  63123346, 69455855, 69201904, 123237848, 136872279, 50317555, 41240941, 104769349),
    end_pos = c(7590856, 43125364, 55207337, 127742951, 25250929, 179240093, 112179823, 87971930,
                48473282, 21995300, 108236235, 39730426, 37050918, 47710367, 10234179, 31327581,
                63320128, 69469242, 69238134, 123357972, 136872279, 50396537, 41281939, 104826359),
    strand = c("-", "-", "+", "+", "-", "+", "+", "+", "+", "-", "+", "+", "-", "+", "+", "-", 
               "-", "+", "+", "+", "+", "+", "+", "+"),
    gene_type = rep("protein_coding", 24),
    description = c("Tumor protein p53", "Breast cancer 1", "Epidermal growth factor receptor", 
                    "MYC proto-oncogene", "KRAS proto-oncogene", "Phosphoinositide-3-kinase",
                    "Adenomatous polyposis coli", "Phosphatase and tensin homolog",
                    "RB transcriptional corepressor 1", "Cyclin dependent kinase inhibitor 2A",
                    "ATM serine/threonine kinase", "erb-b2 receptor tyrosine kinase 2",
                    "mutL homolog 1", "mutS homolog 2", "von Hippel-Lindau tumor suppressor",
                    "Neurofibromin 1", "BCL2 apoptosis regulator", "Cyclin D1", "MDM2 proto-oncogene",
                    "Fibroblast growth factor receptor 2", "Notch receptor 1", "Wnt family member 1",
                    "Catenin beta 1", "AKT serine/threonine kinase 1"),
    stringsAsFactors = FALSE
  )
  
  # Add build-specific adjustments (coordinate differences between builds)
  if (build != "GRCh38") {
    message(paste("Note: Using sample gene annotations. In production, load actual", build, "gene annotations."))
  }
  
  return(sample_genes)
}

#' Annotate CNV Regions with Genes
#' 
#' This function takes CNV results from CytoConverter and annotates them with 
#' genes that overlap the CNV regions.
#' 
#' @param cnv_results A data frame with CNV results from CytoConverter
#' @param build Genome build to use for gene annotations (default: "GRCh38")
#' @param overlap_threshold Minimum overlap percentage to consider a gene affected (default: 0.1)
#' @return A data frame with CNV results annotated with gene information
#' @export
annotate_cnv_with_genes <- function(cnv_results, build = "GRCh38", overlap_threshold = 0.1) {
  
  # Input validation
  if (!is.data.frame(cnv_results)) {
    stop("cnv_results must be a data frame")
  }
  
  required_cols <- c("Sample ID", "Chr", "Start", "End", "Type")
  missing_cols <- setdiff(required_cols, names(cnv_results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  if (!is.numeric(overlap_threshold) || overlap_threshold < 0 || overlap_threshold > 1) {
    stop("overlap_threshold must be a number between 0 and 1")
  }
  
  # Load gene annotations
  tryCatch({
    gene_annotations <- load_gene_annotations(build)
  }, error = function(e) {
    stop(paste("Failed to load gene annotations:", e$message))
  })
  
  # Initialize results list
  annotated_results <- list()
  
  # Process each CNV region
  for (i in seq_len(nrow(cnv_results))) {
    cnv_row <- cnv_results[i, ]
    cnv_chr <- cnv_row[["Chr"]]
    cnv_start <- as.numeric(cnv_row[["Start"]])
    cnv_end <- as.numeric(cnv_row[["End"]])
    
    # Find overlapping genes
    overlapping_genes <- gene_annotations[
      gene_annotations$chromosome == cnv_chr &
      !(gene_annotations$end_pos < cnv_start | gene_annotations$start_pos > cnv_end), 
    ]
    
    if (nrow(overlapping_genes) > 0) {
      # Calculate overlap for each gene
      overlapping_genes$overlap_start <- pmax(overlapping_genes$start_pos, cnv_start)
      overlapping_genes$overlap_end <- pmin(overlapping_genes$end_pos, cnv_end)
      overlapping_genes$overlap_length <- overlapping_genes$overlap_end - overlapping_genes$overlap_start
      overlapping_genes$gene_length <- overlapping_genes$end_pos - overlapping_genes$start_pos
      overlapping_genes$overlap_percentage <- overlapping_genes$overlap_length / overlapping_genes$gene_length
      
      # Filter by overlap threshold
      significant_genes <- overlapping_genes[overlapping_genes$overlap_percentage >= overlap_threshold, ]
      
      if (nrow(significant_genes) > 0) {
        # Create annotated results for each overlapping gene
        for (j in seq_len(nrow(significant_genes))) {
          gene_row <- significant_genes[j, ]
          annotated_row <- data.frame(
            cnv_row,
            Gene_Symbol = gene_row$gene_symbol,
            Gene_Chr = gene_row$chromosome,
            Gene_Start = gene_row$start_pos,
            Gene_End = gene_row$end_pos,
            Gene_Strand = gene_row$strand,
            Gene_Type = gene_row$gene_type,
            Gene_Description = gene_row$description,
            Overlap_Percentage = round(gene_row$overlap_percentage * 100, 2),
            stringsAsFactors = FALSE
          )
          annotated_results[[length(annotated_results) + 1]] <- annotated_row
        }
      } else {
        # No significant gene overlaps
        annotated_row <- data.frame(
          cnv_row,
          Gene_Symbol = "No significant genes",
          Gene_Chr = NA,
          Gene_Start = NA,
          Gene_End = NA,
          Gene_Strand = NA,
          Gene_Type = NA,
          Gene_Description = "No genes with sufficient overlap",
          Overlap_Percentage = 0,
          stringsAsFactors = FALSE
        )
        annotated_results[[length(annotated_results) + 1]] <- annotated_row
      }
    } else {
      # No overlapping genes found
      annotated_row <- data.frame(
        cnv_row,
        Gene_Symbol = "No genes found",
        Gene_Chr = NA,
        Gene_Start = NA,
        Gene_End = NA,
        Gene_Strand = NA,
        Gene_Type = NA,
        Gene_Description = "No genes in this region",
        Overlap_Percentage = 0,
        stringsAsFactors = FALSE
      )
      annotated_results[[length(annotated_results) + 1]] <- annotated_row
    }
  }
  
  # Combine all results
  if (length(annotated_results) > 0) {
    final_results <- do.call(rbind, annotated_results)
    rownames(final_results) <- NULL
    return(final_results)
  } else {
    return(data.frame())
  }
}

#' Extract Genes from Annotated CNV Results
#' 
#' This function extracts a list of unique genes from annotated CNV results
#' for pathway analysis.
#' 
#' @param annotated_cnv A data frame with annotated CNV results
#' @param cnv_type Filter by CNV type ("Gain", "Loss", or "All"). Default: "All"
#' @return A character vector of unique gene symbols
#' @export
extract_genes_from_cnv <- function(annotated_cnv, cnv_type = "All") {
  
  # Input validation
  if (!is.data.frame(annotated_cnv)) {
    stop("annotated_cnv must be a data frame")
  }
  
  if (!"Gene_Symbol" %in% names(annotated_cnv)) {
    stop("Gene_Symbol column not found. Use annotate_cnv_with_genes() first.")
  }
  
  if (!cnv_type %in% c("Gain", "Loss", "All")) {
    stop("cnv_type must be 'Gain', 'Loss', or 'All'")
  }
  
  # Filter by CNV type if specified
  if (cnv_type != "All") {
    if (!"Type" %in% names(annotated_cnv)) {
      stop("Type column not found for filtering")
    }
    annotated_cnv <- annotated_cnv[annotated_cnv$Type == cnv_type, ]
  }
  
  # Extract and clean gene symbols
  genes <- annotated_cnv$Gene_Symbol
  genes <- genes[!is.na(genes)]
  genes <- genes[!genes %in% c("No genes found", "No significant genes")]
  genes <- unique(genes)
  
  return(genes)
}

# Message about successful loading
message("CNV annotation functions loaded successfully.")
message("Available functions:")
message("  - load_gene_annotations()")
message("  - annotate_cnv_with_genes()")
message("  - extract_genes_from_cnv()")