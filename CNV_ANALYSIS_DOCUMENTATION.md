# CNV Annotation and Pathway Analysis for CytoConverter

This documentation provides comprehensive guidance for using the CNV annotation and pathway analysis functionality integrated with CytoConverter.

## Overview

The CNV annotation and pathway analysis extension adds powerful gene-level analysis capabilities to CytoConverter, enabling users to:

1. **Annotate CNV regions** with overlapping genes
2. **Perform pathway enrichment analysis** on affected genes
3. **Visualize results** in a clear, interpretable format
4. **Integrate seamlessly** with existing CytoConverter workflows

## Files and Components

### Core Function Files

- **`cnv_annotation.R`** - Functions for annotating CNV regions with gene information
- **`pathway_analysis.R`** - Functions for pathway enrichment analysis
- **`cnv_analysis_workflow.R`** - Integrated workflow combining all functionality

### Example and Test Files

- **`example_cnv_workflow.R`** - Comprehensive examples demonstrating usage
- **`test_cnv_functions.R`** - Unit tests ensuring code quality
- **`test_cnv_pathway_analysis.R`** - Basic functionality demonstration

## Quick Start

### Basic Usage

```r
# Load required libraries
library(stringr)
library(dplyr)

# Source the integrated workflow
source("cnv_analysis_workflow.R")

# Example 1: Single karyotype analysis
result <- cyto_cnv_analysis(
  input_data = "47,XX,+8",
  build = "GRCh38",
  perform_pathways = TRUE
)

# Display results summary
display_cnv_summary(result)

# Example 2: Multiple sample analysis with CNV data
cnv_data <- data.frame(
  "Sample ID" = c("Patient_1", "Patient_2"),
  "Chr" = c("chr17", "chr8"),
  "Start" = c(7000000, 127000000),
  "End" = c(8000000, 128000000),
  "Type" = c("Loss", "Gain"),
  "Cells Present" = c("45/50", "30/50"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

result <- cyto_cnv_analysis(cnv_data, save_results = TRUE)
```

## Detailed Function Reference

### CNV Annotation Functions (`cnv_annotation.R`)

#### `load_gene_annotations(build = "GRCh38")`

Loads gene annotation data for the specified genome build.

**Parameters:**
- `build`: Genome build ("GRCh38", "hg19", "hg18", "hg17")

**Returns:** Data frame with gene annotations

**Example:**
```r
genes <- load_gene_annotations("GRCh38")
head(genes)
```

#### `annotate_cnv_with_genes(cnv_results, build = "GRCh38", overlap_threshold = 0.1)`

Annotates CNV regions with overlapping genes.

**Parameters:**
- `cnv_results`: Data frame with CNV results from CytoConverter
- `build`: Genome build to use for annotations
- `overlap_threshold`: Minimum overlap percentage to consider a gene affected (0-1)

**Returns:** Data frame with CNV results annotated with gene information

**Example:**
```r
# Assume cnv_data is from CytoConverter
annotated <- annotate_cnv_with_genes(cnv_data, build = "GRCh38")
```

#### `extract_genes_from_cnv(annotated_cnv, cnv_type = "All")`

Extracts gene lists from annotated CNV results.

**Parameters:**
- `annotated_cnv`: Data frame from `annotate_cnv_with_genes()`
- `cnv_type`: Filter by CNV type ("Gain", "Loss", or "All")

**Returns:** Character vector of unique gene symbols

**Example:**
```r
all_genes <- extract_genes_from_cnv(annotated, "All")
gained_genes <- extract_genes_from_cnv(annotated, "Gain")
lost_genes <- extract_genes_from_cnv(annotated, "Loss")
```

### Pathway Analysis Functions (`pathway_analysis.R`)

#### `load_pathway_database(database = "KEGG")`

Loads pathway annotation data from various databases.

**Parameters:**
- `database`: Pathway database ("KEGG", "Reactome", "GO_BP", "GO_MF", "GO_CC")

**Returns:** List containing pathway information

#### `perform_pathway_enrichment(gene_list, database = "KEGG", ...)`

Performs pathway enrichment analysis using hypergeometric test.

**Parameters:**
- `gene_list`: Character vector of gene symbols
- `database`: Pathway database to use
- `background_size`: Total genes in genome (default: 20000)
- `min_pathway_size`: Minimum pathway size (default: 5)
- `max_pathway_size`: Maximum pathway size (default: 500)
- `p_cutoff`: P-value cutoff for significance (default: 0.05)

**Returns:** Data frame with enrichment results

**Example:**
```r
# Perform KEGG pathway enrichment
kegg_results <- perform_pathway_enrichment(
  gene_list = c("TP53", "BRCA1", "EGFR"),
  database = "KEGG",
  p_cutoff = 0.05
)
```

#### `visualize_pathway_results(enrichment_results, top_n = 10)`

Creates formatted visualization of pathway enrichment results.

**Parameters:**
- `enrichment_results`: Data frame from `perform_pathway_enrichment()`
- `top_n`: Number of top pathways to display

**Example:**
```r
visualize_pathway_results(kegg_results, top_n = 5)
```

#### `pathway_enrichment_summary(gene_list, databases = c("KEGG", "Reactome", "GO_BP"))`

Performs comprehensive pathway analysis across multiple databases.

**Parameters:**
- `gene_list`: Character vector of gene symbols
- `databases`: Vector of database names to analyze
- `p_cutoff`: P-value cutoff for significance

**Returns:** List containing results for each database

### Integrated Workflow Functions (`cnv_analysis_workflow.R`)

#### `cyto_cnv_analysis(input_data, build = "GRCh38", ...)`

Complete workflow from karyotype input to pathway analysis.

**Parameters:**
- `input_data`: Karyotype string or CNV data frame
- `build`: Genome build to use
- `perform_pathways`: Whether to perform pathway analysis (default: TRUE)
- `pathway_databases`: Databases to analyze (default: c("KEGG", "Reactome", "GO_BP"))
- `save_results`: Whether to save results to files (default: FALSE)
- `output_prefix`: Prefix for output files (default: "cnv_analysis")

**Returns:** List containing all analysis results

#### `display_cnv_summary(analysis_results, show_pathways = TRUE)`

Displays formatted summary of CNV analysis results.

**Parameters:**
- `analysis_results`: Results from `cyto_cnv_analysis()`
- `show_pathways`: Whether to display pathway results
- `max_pathways`: Maximum pathways to display per database

## Use Cases and Examples

### Use Case 1: Basic CNV Gene Annotation

```r
# Load functions
source("cnv_annotation.R")

# Create or load CNV data
cnv_data <- data.frame(
  "Sample ID" = "Patient_001",
  "Chr" = "chr17",
  "Start" = 7000000,
  "End" = 8000000,
  "Type" = "Loss",
  "Cells Present" = "45/50",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Annotate with genes
annotated <- annotate_cnv_with_genes(cnv_data)

# View results
print(annotated[, c("Sample.ID", "Chr", "Type", "Gene_Symbol", "Gene_Description")])
```

### Use Case 2: Pathway Enrichment Analysis

```r
# Load functions
source("pathway_analysis.R")

# Gene list from previous annotation
genes <- c("TP53", "BRCA1", "EGFR", "MYC")

# Analyze different pathway databases
kegg_pathways <- perform_pathway_enrichment(genes, "KEGG")
reactome_pathways <- perform_pathway_enrichment(genes, "Reactome")
go_pathways <- perform_pathway_enrichment(genes, "GO_BP")

# Visualize results
visualize_pathway_results(kegg_pathways, top_n = 5)
```

### Use Case 3: Complete Integrated Analysis

```r
# Load integrated workflow
source("cnv_analysis_workflow.R")

# Comprehensive analysis
results <- cyto_cnv_analysis(
  input_data = "46,XX,del(17)(p13.1)",
  build = "GRCh38",
  perform_pathways = TRUE,
  pathway_databases = c("KEGG", "Reactome", "GO_BP"),
  save_results = TRUE,
  output_prefix = "patient_analysis"
)

# Display summary
display_cnv_summary(results)

# Access specific results
affected_genes <- results$gene_lists$all_genes
kegg_results <- results$pathway_results$KEGG
```

### Use Case 4: Batch Analysis

```r
# Multiple samples
sample_data <- data.frame(
  "Sample ID" = rep(c("Sample_A", "Sample_B", "Sample_C"), each = 2),
  "Chr" = c("chr17", "chr8", "chr7", "chr9", "chr3", "chr12"),
  "Start" = c(7000000, 127000000, 55000000, 21900000, 179000000, 25000000),
  "End" = c(8000000, 128000000, 56000000, 22000000, 180000000, 26000000),
  "Type" = c("Loss", "Gain", "Gain", "Loss", "Gain", "Gain"),
  "Cells Present" = rep("unknown", 6),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Analyze all samples together
batch_results <- cyto_cnv_analysis(sample_data)

# Or analyze by sample
for (sample_id in unique(sample_data$`Sample ID`)) {
  sample_subset <- sample_data[sample_data$`Sample ID` == sample_id, ]
  sample_results <- cyto_cnv_analysis(sample_subset)
  cat(sprintf("\nResults for %s:\n", sample_id))
  display_cnv_summary(sample_results, show_pathways = FALSE)
}
```

## Output Files

When `save_results = TRUE`, the following files are created:

- **`{prefix}_cnv_results.txt`** - Original CNV results
- **`{prefix}_annotated_cnv.txt`** - CNV results with gene annotations
- **`{prefix}_all_genes.txt`** - All affected genes
- **`{prefix}_gained_genes.txt`** - Genes in gained regions
- **`{prefix}_lost_genes.txt`** - Genes in lost regions
- **`{prefix}_{database}_pathways.txt`** - Pathway enrichment results for each database

## Data Requirements

### Input CNV Data Format

The CNV data should be a data frame with these columns:

| Column | Description | Required |
|--------|-------------|----------|
| Sample ID | Sample identifier | Yes |
| Chr | Chromosome (e.g., "chr1", "chrX") | Yes |
| Start | Start position (numeric) | Yes |
| End | End position (numeric) | Yes |
| Type | CNV type ("Gain" or "Loss") | Yes |
| Cells Present | Cell count information | No |

### Supported Genome Builds

- **GRCh38** (default) - Latest human genome reference
- **hg19** - Previous standard human genome reference
- **hg18** - Older human genome reference
- **hg17** - Historical human genome reference

## Performance Considerations

- **Gene annotation**: Fast for typical CNV analyses (< 1000 regions)
- **Pathway analysis**: Scales with number of genes (typically fast for < 100 genes)
- **Memory usage**: Minimal for typical use cases
- **File I/O**: Optional file output for large analyses

## Troubleshooting

### Common Issues

1. **"Column not found" errors**
   - Ensure CNV data has required columns
   - Check for spaces in column names (use `check.names = FALSE`)

2. **"No genes found" results**
   - Verify chromosome naming convention (should include "chr" prefix)
   - Check coordinate system (1-based genomic coordinates)
   - Consider lowering `overlap_threshold` parameter

3. **"No significant pathways" results**
   - Increase `p_cutoff` parameter
   - Try different pathway databases
   - Ensure gene symbols are correctly formatted

4. **Memory issues**
   - Reduce `background_size` parameter
   - Process samples in smaller batches
   - Use `save_results = FALSE` to reduce memory usage

### Getting Help

1. **Check function documentation**: All functions include detailed parameter descriptions
2. **Run unit tests**: Execute `test_cnv_functions.R` to verify installation
3. **Review examples**: See `example_cnv_workflow.R` for comprehensive usage examples
4. **Validate input data**: Ensure data format matches requirements

## Best Practices

1. **Always validate input data** before analysis
2. **Use appropriate genome build** matching your coordinate system
3. **Save results for large analyses** to avoid re-computation
4. **Consider multiple pathway databases** for comprehensive analysis
5. **Adjust significance thresholds** based on study design
6. **Document analysis parameters** for reproducibility

## Integration with CytoConverter

This extension is designed to work seamlessly with the existing CytoConverter workflow:

1. **Input compatibility**: Accepts CytoConverter output format directly
2. **Coordinate system**: Uses same genomic coordinate system
3. **Genome builds**: Supports all CytoConverter genome builds
4. **Error handling**: Graceful fallback when CytoConverter is unavailable

The integration allows for complete karyotype-to-pathway analysis in a single workflow, making it easy to go from cytogenetic descriptions to biological interpretation.