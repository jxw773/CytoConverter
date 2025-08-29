# CytoConverter

ISB-CGC-CytoConverter code is modified from a fork of the CytoConverter project: [https://github.com/jxw773/CytoConverter](https://github.com/jxw773/CytoConverter)

[CytoConverter: a web-based tool to convert karyotypes to genomic coordinates](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3062-4)

Cytogenetic nomenclature is used to describe chromosomal aberrations (or lack thereof) in a collection of cells, referred to as the cells’ karyotype. The nomenclature identifies locations on chromosomes using a system of cytogenetic bands, each with a unique name and region on a chromosome. Each band is microscopically visible after staining, and encompasses a large portion of the chromosome. More modern analyses employ genomic coordinates, which precisely specify a chromosomal location according to its distance from the end of the chromosome. Currently, there is no tool to convert cytogenetic nomenclature into genomic coordinates. Since locations of genes and other genomic features are usually specified by genomic coordinates, a conversion tool will facilitate the identification of the features that are harbored in the regions of chromosomal gain and loss that are implied by a karyotype.

## Requirements

CytoConverter requires R 4.0+. Before running the main script, make sure that required R packages
are installed by changing to the CytoConverter directory and running:

```
./init.R
```

## Running CytoConverter

### Quick Start

1. **Install dependencies:**
   ```bash
   ./init.R
   ```

2. **Run with command-line interface:**
   ```bash
   ./cytoconverter \
     --input input-file.txt \
     --output output-file.txt \
     --log log-file.txt
   ```

3. **Use directly in R:**
   ```r
   # Load the module
   library(modules)
   mod_cytoconverter <- use('modules/cytoconverter.R')
   
   # Convert a single karyotype
   result <- mod_cytoconverter$CytoConverter("47,XY,+8")
   
   # Convert a table of samples
   data <- matrix(c("Sample1", "47,XY,+8", "Sample2", "46,XX"), ncol=2, byrow=T)
   result <- mod_cytoconverter$CytoConverter(data)
   
   # Access results
   gains_losses <- result$Result
   errors <- result$Error_log
   ```

### Command Line Parameters

- **input**: Input file of sample names and karyotypes (tab-delimited, one per line)
- **output**: Output file containing genomic coordinates for gains/losses  
- **log**: Log file for warnings and errors encountered during processing

### Input Format

Input files should be tab-delimited with two columns:
```
Sample_ID    Karyotype
ABC          45,XY,der(1;19)(q10;p10)
DEF          47,XX,+der(10)t(10;21)(p13;q21)
P10          47,X,+X[30]/48,XX,+7,+9[50]
```


## Code Structure

The code has been split into multiple "modules" organized as follows:

### Core Modules
- **modules/cytoconverter.R**: Main entry point containing the `CytoConverter()` function
- **modules/rowparser.R**: Handles parsing of individual rows from karyotype tables  
- **modules/colparser.R**: Processes individual karyotype components and cell lines
- **modules/merge.R**: Utilities for handling and merging genomic intervals
- **modules/utils.R**: General utility functions for data processing
- **modules/cytobands.R**: Functions for retrieving cytoband information
- **modules/gainlossfusion.R**: Processes gain/loss events and fusion detection
- **modules/mergefusions.R**: Specialized fusion event handling

### Supporting Files
- **main.R**: Command-line wrapper script
- **init.R**: Dependency installation script
- **plot_cyto_graph.R**: Visualization functions for plotting results
- **cyto_graph.R**: Core graphing utilities

### Module Dependencies
```
   ┌─────────────────────┐      ┌────────────────────────────┐
   │ modules/            │      │ modules/merge.R            │
   │ cytoconverter.R     │      │                            │
   │                     │      │    insertSection()         │
   │  CytoConverter()    │      │    deleteIntersections()   │
   └───┬─────────────────┘      │    getContiguousSection()  │
       │                        │    mergeAdjacentSections() │
       │                        │    mergeTable()            │
   ┌───▼─────────────────┐      │    mergeDel()              │
   │ modules/rowparser.R ├──────►    mergeDelmat()           │
   │                     │      │    bigDelMerge()           │
   │      rowparse()     │      │    mergeDeletions()        │
   └───┬─────────────────┘      └────────────────────────────┘
       │        │              
       │        │               ┌─────────────────────────────┐
       │        └───────────────► modules/utils.R             │
       │                        │                             │
       │                        │    positionSorter()         │
   ┌───▼─────────────────┐      │    mergeIntOverlap()        │
   │ modules/colparser.R ├──────►    detectAdd()              │
   │                     │      └─────────────────────────────┘
   │     colparse()      │
   │     miniverter()    │      ┌─────────────────────────────┐
   └───┬─────────────────┘      │ modules/cytobands.R         │
       │                        │                             │
       └────────────────────────►    getCytoBands()           │
                                └─────────────────────────────┘
```

## CytoConverter Function Parameters

The `CytoConverter()` function accepts the following parameters:

### Required Parameters
- **in_data**: Input karyotype string or matrix with sample names (col 1) and karyotypes (col 2)

### Optional Parameters  
- **build**: Reference genome build - `"GRCh38"` (default), `"hg19"`, `"hg18"`, or `"hg17"`
- **constitutional**: Boolean - whether to include constitutional changes (`TRUE`)
- **guess**: Boolean - attempt to interpret ambiguous karyotypes (`FALSE`)
- **guess_q**: Boolean - process karyotypes with question marks (`FALSE`) 
- **forMtn**: Boolean - optimize for Mitelman database format (`TRUE`)
- **orOption**: Boolean - when "or" appears, take first option (`TRUE`)
- **sexstimate**: Boolean - estimate sex chromosome composition (`FALSE`)
- **allow_Shorthand**: Boolean - allow shorthand notation outside clones (`FALSE`)
- **count_fusions**: Boolean - detect and count fusion events (`FALSE`)
- **include_normals_graph**: Boolean - include normal samples in output (`FALSE`)

### Example Usage
```r
# Basic usage with defaults
result <- CytoConverter("46,XY,t(9;22)(q34;q11)")

# Advanced usage with custom parameters  
result <- CytoConverter(
  in_data = karyotype_matrix,
  build = "hg19", 
  constitutional = FALSE,
  guess = TRUE,
  count_fusions = TRUE
)
```

### Return Value
The function returns a list containing:
- **Result**: Data frame with genomic coordinates of gains/losses
- **Error_log**: Data frame with warnings and errors  
- **Fusion_table**: Data frame with fusion events (if `count_fusions=TRUE`)
- **List_of_samples**: Sample identifiers (if `include_normals_graph=TRUE`)

## Output Format

### Results Table Structure
The results table contains the following columns:
- **Sample ID**: Original sample identifier
- **Chr**: Chromosome (e.g., "chr1", "chrX")
- **Start**: Start genomic coordinate 
- **End**: End genomic coordinate
- **Type**: Event type ("Gain" or "Loss")
- **Percent Present**: Proportion of cells affected (e.g., "30 of 80")

### Example Output
```
Sample ID    Chr     Start      End        Type    Percent Present
ABC_1        chr1    0          249250621  Loss    unknown
DEF_1        chr10   39254935   135534747  Gain    unknown
P10_1        chrX    0          155270560  Gain    30 of 80
P10_2        chr7    0          159138663  Gain    50 of 80
```

## Visualization

Built-in functions for creating graphs displaying gains and losses:

```r
# Source the graphing functions
source("plot_cyto_graph.R")
source("cyto_graph.R")

# Create visualization
plot_cyto_graph(
  cyto_list = result$Result,      # Results from CytoConverter
  ref_list = "GRCh38",           # Reference genome build  
  ylabel = TRUE                   # Show sample names on graph
)
```

### Visualization Parameters
- **cyto_list**: Results table from CytoConverter output
- **list_from_cyto**: Pre-processed output from cyto_graph (optional)
- **ref_list**: Reference genome build for coordinate mapping
- **ylabel**: Whether to display sample names on the y-axis

## Troubleshooting

### Common Issues

**Missing dependencies**: Run `./init.R` to install required R packages

**Input format errors**: Ensure input is tab-delimited with sample names in column 1 and karyotypes in column 2

**Memory issues**: For large datasets, consider processing in smaller batches

**Complex karyotypes**: Use `guess=TRUE` to attempt parsing of ambiguous notation

### Error Messages
- **"build incorrectly specified"**: Use one of: GRCh38, hg19, hg18, hg17
- **"Warning in karyotype number not specified"**: Karyotype missing chromosome count
- **"Error in markers and other ambiguous objects"**: Contains unrecognizable elements

For additional support, check the Error_log output for detailed diagnostic information.

