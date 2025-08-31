# CytoConverter API Reference

This document provides a comprehensive reference for CytoConverter functions and their usage.

## Core Functions

### CytoConverter()

**Main function for cytogenetic analysis**

```r
CytoConverter(in_data, build = "GRCh38", constitutional = TRUE, guess = FALSE, 
              guess_q = FALSE, guess_by_first_val = FALSE, forMtn = TRUE, 
              orOption = TRUE, sexstimate = FALSE, allow_Shorthand = FALSE, 
              count_fusions = FALSE, include_normals_graph = FALSE)
```

**Parameters:**
- `in_data` - Input karyotype data (string or matrix/data.frame)
- `build` - Reference genome build ("GRCh38", "hg19", "hg18", "hg17")
- `constitutional` - Boolean for constitutional analysis mode
- `guess` - Boolean to enable guessing of ambiguous regions
- `guess_q` - Boolean for q-arm specific guessing
- `guess_by_first_val` - Boolean to guess based on first values
- `forMtn` - Boolean for Montreal nomenclature compatibility
- `orOption` - Boolean to enable OR logic in parsing
- `sexstimate` - Boolean for sex chromosome estimation
- `allow_Shorthand` - Boolean to allow shorthand notation
- `count_fusions` - Boolean to enable fusion detection and analysis
- `include_normals_graph` - Boolean for including normal samples in fusion graphs

**Returns:**
List with `Results` (data.frame) and `Error_log` (data.frame)

**Example:**
```r
# Basic usage
result <- CytoConverter("46,XY,del(7)(q22q32)")

# With fusion analysis
result <- CytoConverter(karyotype_data, count_fusions = TRUE)
```

## Visualization Functions

### plot_cyto_graph()

**Standard plotting for gains and losses**

```r
plot_cyto_graph(cyto_list = NULL, list_from_cyto = NULL, ref_list = "GRCh38", 
                ylabel = NULL, include_normals_graph = FALSE, list_of_samples = NULL)
```

**Parameters:**
- `cyto_list` - CytoConverter results data frame
- `list_from_cyto` - Pre-computed cyto_graph output (optional)
- `ref_list` - Reference genome build
- `ylabel` - Boolean to show sample labels (auto-determined if NULL)
- `include_normals_graph` - Boolean to include normal samples
- `list_of_samples` - Vector of normal sample names

**Example:**
```r
result <- CytoConverter(data)
plot_cyto_graph(result$Results)
```

### cyto_graph_fusion()

**Specialized plotting for fusion data**

```r
cyto_graph_fusion(cyto_list, ref_list = "GRCh38", include_normals_graph = FALSE, 
                  list_of_samples = NULL)
```

**Use for:** Visualizing structural rearrangements and fusion events detected with `count_fusions = TRUE`

### cyto_graph()

**Data preparation for standard plotting**

```r
cyto_graph(cyto_list, ref_list = "GRCh38", include_normals_graph = FALSE, 
           list_of_samples = NULL)
```

**Note:** Typically called internally by `plot_cyto_graph()`

## Data Processing Functions

### rowparse()

**Internal function for parsing karyotype rows**

Converts individual karyotype strings into structured genomic coordinate data.

### colparse()

**Internal function for parsing cytogenetic components**

Handles detailed interpretation of complex cytogenetic aberrations.

### miniverter()

**Internal triage function for simple vs complex aberrations**

Determines whether components require simple or complex parsing logic.

## Utility Functions

### positionSorter()

**Sorts cytogenetic band positions**

```r
positionSorter(positions)
```

### mergeIntOverlap()

**Merges overlapping genomic intervals**

```r
mergeIntOverlap(v1, v2)
```

### detectAdd()

**Detects addition/deletion patterns**

```r
detectAdd(temp_table_processed, ex_table_processed)
```

## Installation Functions

### install_libraries()

**Installs required R packages**

```r
install_libraries()
```

Called automatically by `./init.R` script.

**Required packages:**
- modules
- stringr
- stringi  
- DescTools
- dplyr
- hash
- optparse

## Input/Output Formats

### Input Format

**Single karyotype string:**
```r
CytoConverter("46,XY,del(7)(q22q32)")
```

**Table format (tab-delimited):**
```
Sample1    46,XY,del(7)(q22q32)
Sample2    47,XY,+21
Sample3    46,XX,t(9;22)(q34;q11.2)
```

### Output Format

**Standard output columns:**
- Sample ID
- Clone number  
- Chromosome
- Start coordinate
- End coordinate
- Type (Gain/Loss or fusion tags)
- Percent present

**Example output:**
```
Sample1  1  chr7  116000000  130000000  Loss  1/10
Sample2  1  chr21  15000000  48100000   Gain  1/10
```

**Fusion output includes tags:**
```
Sample3  1  chr9   130854089  130920000  #translocation_balanced|fus_1  1/10
Sample3  1  chr22  23500000   23632600   #translocation_balanced|fus_2  1/10
```

## Fusion Analysis

### Supported Fusion Types

- **Translocations:** `t(9;22)(q34;q11.2)`, `der(10)t(10;21)(p13;q21)`
- **Inversions:** `inv(16)(p13.1q22)`
- **Ring chromosomes:** `r(7)(p22q36)`
- **Dicentric chromosomes:** `dic(X;Y)(p22.3;p11.3)`
- **Robertsonian translocations:** `rob(13;14)(q10;q10)`
- **And many others** - see Fusion_tags_WIP.txt for complete list

### Fusion Tags

Fusion events are marked with descriptive tags:
- `#translocation_balanced` - Balanced translocations
- `#derivative_chrom` - Derivative chromosomes
- `#ring_chrom` - Ring chromosomes
- `#dicentric` - Dicentric chromosomes

## Reference Genome Builds

Supported builds with 850-resolution cytobands:
- **GRCh38** (default) - Latest human genome reference
- **hg19** - Previous standard reference
- **hg18** - Older reference build
- **hg17** - Legacy reference build

Build files are located in `Builds/` directory:
- `cytoBand_GRCh38.txt`
- `cytoBand_hg19.txt`  
- `cytoBand_hg18.txt`
- `cytoBand_hg17.txt`

## Error Handling

CytoConverter provides comprehensive error logging:

```r
result <- CytoConverter(data)
errors <- result$Error_log  # Check for parsing issues
results <- result$Results   # Main results table
```

Common error types:
- Malformed karyotype strings
- Unrecognized cytogenetic notation
- Invalid chromosome names
- Conflicting aberration descriptions

## Performance Considerations

- **Large datasets:** Use command-line interface for better memory management
- **Fusion analysis:** Significantly more computationally intensive than standard analysis
- **Threading:** Command-line interface supports parallel processing
- **Memory:** Fusion analysis requires more memory for complex karyotypes

## See Also

- [README.md](README.md) - General usage and examples
- [Fusion_tags_WIP.txt](Fusion_tags_WIP.txt) - Detailed fusion tag documentation
- Module files in `modules/` directory for implementation details