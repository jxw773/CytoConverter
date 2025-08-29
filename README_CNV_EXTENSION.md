# CytoConverter CNV Analysis Extension

## Quick Start Guide

This extension adds gene annotation and pathway analysis capabilities to CytoConverter.

### Installation

Ensure you have the required R packages:
```r
# Required packages
install.packages(c("stringr", "stringi", "dplyr"))
```

### Basic Usage

```r
# Load the integrated workflow
source("cnv_analysis_workflow.R")

# Analyze a karyotype
result <- cyto_cnv_analysis("47,XX,+8")

# Display results
display_cnv_summary(result)
```

### Key Features

- ✅ **Gene Annotation**: Annotate CNV regions with overlapping genes
- ✅ **Pathway Analysis**: Perform enrichment analysis using KEGG, Reactome, GO
- ✅ **Multiple Genome Builds**: Support for GRCh38, hg19, hg18, hg17
- ✅ **Comprehensive Documentation**: Detailed examples and error handling
- ✅ **Unit Tests**: Robust testing ensures reliability
- ✅ **Integration**: Seamless integration with existing CytoConverter

### Files Included

| File | Description |
|------|-------------|
| `cnv_annotation.R` | Core gene annotation functions |
| `pathway_analysis.R` | Pathway enrichment analysis functions |
| `cnv_analysis_workflow.R` | Integrated workflow |
| `example_cnv_workflow.R` | Comprehensive usage examples |
| `test_cnv_functions.R` | Unit tests |
| `CNV_ANALYSIS_DOCUMENTATION.md` | Complete documentation |

### Example Output

```
===============================================
       CNV ANALYSIS RESULTS SUMMARY
===============================================

BASIC STATISTICS:
- Total CNV regions: 5
- Total gene annotations: 5
- Unique genes affected: 5
- Chromosomes affected: chr17, chr8, chr7, chr3, chr9
- Genome build: GRCh38

AFFECTED GENES:
- All genes (5): TP53, MYC, EGFR, PIK3CA, CDKN2A
- Gained genes (3): MYC, EGFR, PIK3CA
- Lost genes (2): TP53, CDKN2A

PATHWAY ENRICHMENT RESULTS:

KEGG Database (4 significant pathways):
  1. Pathways in cancer (p=3.78e-17, 5 genes)
  2. Cell cycle (p=5.00e-07, 2 genes)
  3. PI3K-Akt signaling pathway (p=5.00e-07, 2 genes)
  4. Colorectal cancer (p=5.00e-07, 2 genes)
```

See `CNV_ANALYSIS_DOCUMENTATION.md` for complete documentation.