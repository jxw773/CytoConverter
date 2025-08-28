# CytoConverter

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

CytoConverter is a powerful R tool that converts cytogenetic nomenclature (karyotypes) into genomic coordinates. This enables researchers to bridge the gap between traditional cytogenetic analysis and modern genomic approaches.

## Overview

**ISB-CGC-CytoConverter** is modified from a fork of the [CytoConverter project](https://github.com/jxw773/CytoConverter).

**Research Paper**: [CytoConverter: a web-based tool to convert karyotypes to genomic coordinates](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3062-4)

### What is CytoConverter?

Cytogenetic nomenclature describes chromosomal aberrations using a system of cytogenetic bands - regions on chromosomes that are microscopically visible after staining. Modern genomic analyses use precise genomic coordinates that specify chromosomal locations by distance from chromosome ends.

CytoConverter fills this critical gap by:
- **Converting** cytogenetic band notation to precise genomic coordinates
- **Identifying** regions of chromosomal gain and loss from karyotype data
- **Facilitating** integration with modern genomic feature databases
- **Supporting** multiple genome builds (GRCh38, hg19, hg18, hg17)

## Requirements

**Prerequisites:**
- R 4.0 or higher
- Required R packages (automatically installed via init script)

**Installation:**

1. Clone or download the CytoConverter repository
2. Navigate to the CytoConverter directory 
3. Install required dependencies:

```bash
./init.R
```

## Usage

### Command Line Interface

Run CytoConverter with the wrapper script:

```bash
./cytoconverter \
  --input input-file.txt \
  --threads 4 \
  --output output-file.txt \
  --log log-file.txt
```

**Parameters:**
- `input`: Input file of sample names and associated karyotypes, one per line, tab delimited
- `threads`: Number of parallel threads to run (input file will be split accordingly)
- `output`: Output file containing genomic coordinates and gain/loss indicators for all samples
- `log`: Log file containing warnings or errors encountered during processing

### R API Usage

```r
# Load the CytoConverter function
source("modules/cytoconverter.R")

# Basic usage with karyotype string
result <- CytoConverter("46,XY,+21")

# Usage with data table
karyotype_data <- data.frame(
  Sample = c("Sample1", "Sample2"),
  Karyotype = c("46,XY,+21", "46,XX,del(5q)")
)
result <- CytoConverter(karyotype_data)

# Access results
coordinates_table <- result$Results
error_log <- result$Error_log
```

**Function Parameters:**
- `in_data`: Input karyotype string or data frame
- `build`: Genome build ("GRCh38", "hg19", "hg18", "hg17") - default: "GRCh38"
- `constitutional`: Include constitutional variations - default: TRUE
- `guess`: Attempt to guess ambiguous notations - default: FALSE
- `guess_q`: Guess '?' marks in karyotypes - default: FALSE
- `forMtn`: Include Mountain regions - default: TRUE
- `orOption`: Handle 'or' statements - default: TRUE
- `sexstimate`: Estimate sex from karyotype - default: FALSE
- `allow_Shorthand`: Allow shorthand notation - default: FALSE

### Input Format

#### Karyotype String Format
```
46,XY,+21
46,XX,del(5q13q33)
47,XY,+8,der(16)t(1;16)(q23;q24)
```

#### Table Format
Tab-delimited file with sample names and karyotypes:
```
Sample1	46,XY,+21
Sample2	46,XX,del(5q13q33)
Sample3	47,XY,+8,der(16)t(1;16)(q23;q24)
```

### Output Format

The results table contains:
- **Sample**: Sample identifier
- **Clone**: Clone line number
- **Start**: Start genomic coordinate of gain/loss
- **End**: End genomic coordinate of gain/loss  
- **Type**: Indicator for "Gain" or "Loss"
- **CellCount**: Number of cells in clone out of total cells in sample

## Code Structure

The code has been modularized for better maintainability:

```
   ┌────────────────────┐      ┌───────────────────────────┐
   │ cytoconverter.R    │      │ merge.R                   │
   │                    │      │                           │
   │    CytoConverter() │      │    insertSection()        │
   └───┬────────────────┘      │    deleteIntersections()  │
       │                       │    getContiguousSection() │
       │                       │    mergeAdjacentSections()│
       │                       │    mergeTable()           │
   ┌───▼───────────┐           │    mergeDel()             │
   │ rowparser.R   ├───────────►    mergeDelmat()          │
   │               │           │    bigDelMerge()          │
   │    rowparse() │           │    mergeDeletions()       │
   └───┬────────┬──┘           └───────────────────────────┘
       │        │              ┌──────────────────────┐
       │        └──────────────► utils.R              │
       │                       │                      │
       │                       │    positionSorter()  │
   ┌───▼────────────┐          │    mergeIntOverlap() │
   │ colparser.R    ├──────────►    detectAdd()       │
   │                │          └──────────────────────┘
   │     colparse() │
   └────────────┬───┘          ┌───────────────────┐
                │              │ cytobands.R       │
                └──────────────►                   │
                               │    getCytoBands() │
                               └───────────────────┘
```

**Module Descriptions:**
- `cytoconverter.R`: Main entry point and core conversion logic
- `rowparser.R`: Handles parsing of individual karyotype rows
- `colparser.R`: Processes individual karyotype components and chromosomal aberrations
- `merge.R`: Provides functions for merging and handling genomic intervals
- `utils.R`: Utility functions for position sorting and overlap detection
- `cytobands.R`: Manages cytoband data for translocations and insertions

## Visualization

Built-in graphing functions are available:

```r
# Source plotting functions
source("plot_cyto_graph.R")
source("cyto_graph.R")

# Create visualization
plot_cyto_graph(
  cyto_list = result$Results,
  ref_list = "GRCh38",
  ylabel = TRUE
)
```

**Plotting Parameters:**
- `cyto_list`: Results table from CytoConverter
- `list_from_cyto`: Alternative input from cyto_graph function (optional)
- `ref_list`: Reference genome for plotting coordinates (default: "GRCh38")
- `ylabel`: Enable/disable sample names on graph (default: TRUE)

## Supported Genome Builds

CytoConverter supports multiple genome builds at 850 resolution:
- **GRCh38** (default)
- **hg19**
- **hg18** 
- **hg17**

Custom cytoband lists can be provided if needed.

## Examples

### Basic Conversion
```r
# Simple trisomy 21
result <- CytoConverter("47,XY,+21")
print(result$Results)
```

### Batch Processing
```r
# Multiple samples
samples <- data.frame(
  Sample = c("Patient1", "Patient2", "Patient3"),
  Karyotype = c(
    "46,XY,del(5q13q33)",
    "47,XX,+21", 
    "46,XY,t(9;22)(q34;q11)"
  )
)
result <- CytoConverter(samples, build = "GRCh38")
```

### Error Handling
```r
result <- CytoConverter("invalid_karyotype")
if (nrow(result$Error_log) > 0) {
  print("Errors encountered:")
  print(result$Error_log)
}
```

## Contributing

When contributing to this project:
1. Follow the existing code style and structure
2. Add appropriate documentation for new functions
3. Test changes with various karyotype formats
4. Update this README if adding new features

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use CytoConverter in your research, please cite:

> CytoConverter: a web-based tool to convert karyotypes to genomic coordinates. BMC Bioinformatics 20, 467 (2019). https://doi.org/10.1186/s12859-019-3062-4