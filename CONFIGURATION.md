# CytoConverter Configuration Guide

This guide explains the various configuration options and parameters available in CytoConverter for different analysis scenarios.

## Analysis Mode Configuration

### Constitutional vs Somatic Analysis

**Constitutional Analysis** (constitutional = TRUE)
- Use for germline karyotypes and constitutional chromosomal disorders
- Handles constitutional cytogenetic notation (e.g., markers ending with 'c')
- More conservative interpretation of ambiguous regions
- Default for interactive R usage

**Somatic Analysis** (constitutional = FALSE)  
- Use for cancer karyotypes and acquired chromosomal aberrations
- Optimized for complex tumor karyotypes with multiple aberrations
- More permissive interpretation allowing for tumor evolution patterns
- Default for command-line interface

```r
# Constitutional analysis
result <- CytoConverter(data, constitutional = TRUE)

# Somatic analysis  
result <- CytoConverter(data, constitutional = FALSE)
```

## Ambiguous Region Handling

### Guessing Parameters

**Basic Guessing** (guess = FALSE/TRUE)
- Controls overall guessing behavior for ambiguous cytogenetic regions
- When TRUE: Attempts to resolve unclear breakpoints using heuristics
- When FALSE: Reports uncertain regions as errors

**Q-arm Guessing** (guess_q = FALSE/TRUE)  
- Specifically handles ambiguous q-arm regions marked with '?'
- When TRUE: Strips '?' symbols and attempts coordinate assignment
- Use with caution as it may introduce false precision

**First Value Guessing** (guess_by_first_val = FALSE/TRUE)
- Uses first coordinate when ranges are ambiguous
- Provides consistent but potentially inaccurate coordinate assignment
- Useful for standardizing output when exact precision isn't critical

```r
# Conservative approach - report ambiguities as errors
result <- CytoConverter(data, guess = FALSE, guess_q = FALSE)

# Permissive approach - attempt to resolve ambiguities  
result <- CytoConverter(data, guess = TRUE, guess_q = TRUE, guess_by_first_val = TRUE)
```

## Nomenclature Compatibility

### Montreal vs ISCN Standards

**Montreal Nomenclature** (forMtn = TRUE)
- Compatible with Montreal chromosomal nomenclature conventions
- Handles specific Montreal notation variations
- Default setting for broad compatibility

**OR Logic Handling** (orOption = TRUE/FALSE)
- Controls handling of "or" statements in karyotypes (e.g., "del(5q) or del(7q)")
- When TRUE: Takes first option from OR statements
- When FALSE: May report OR statements as parsing errors

**Shorthand Notation** (allow_Shorthand = FALSE/TRUE)
- Enables parsing of abbreviated cytogenetic notation
- Use when input contains non-standard shorthand forms
- May increase false positive rate

```r
# Standard ISCN-compatible parsing
result <- CytoConverter(data, forMtn = TRUE, orOption = TRUE, allow_Shorthand = FALSE)

# Permissive parsing for varied notation styles
result <- CytoConverter(data, forMtn = TRUE, orOption = TRUE, allow_Shorthand = TRUE)
```

## Sex Chromosome Analysis

### Sex Chromosome Estimation (sexstimate = FALSE/TRUE)

When enabled, provides enhanced analysis of sex chromosome aberrations:
- Estimates normal sex chromosome complements
- Adjusts gain/loss calling for X and Y chromosomes  
- Accounts for sex chromosome dosage compensation
- Particularly useful for constitutional analysis

```r
# Enhanced sex chromosome analysis
result <- CytoConverter(data, sexstimate = TRUE, constitutional = TRUE)
```

## Fusion Analysis Configuration

### Basic Fusion Detection (count_fusions = FALSE/TRUE)

**Standard Analysis** (count_fusions = FALSE)
- Focuses on gains and losses
- Faster processing for large datasets
- Suitable for copy number analysis

**Fusion Analysis** (count_fusions = TRUE)
- Detects structural rearrangements and chromosomal fusions
- Provides detailed breakpoint information
- Applies fusion-specific tags to results
- Significantly more computationally intensive

```r
# Standard gain/loss analysis
result <- CytoConverter(data, count_fusions = FALSE)

# Comprehensive fusion analysis
result <- CytoConverter(data, count_fusions = TRUE)
```

### Fusion Visualization (include_normals_graph = FALSE/TRUE)

Controls inclusion of normal samples in fusion-specific plots:
- Useful for comparative analysis with control samples
- Requires list_of_samples parameter when enabled
- Only applicable when using cyto_graph_fusion()

```r
# Fusion analysis with normal sample comparison
result <- CytoConverter(data, count_fusions = TRUE, include_normals_graph = TRUE)
plot_data <- cyto_graph_fusion(result$Results, include_normals_graph = TRUE, 
                              list_of_samples = c("Normal1", "Normal2"))
```

## Reference Genome Configuration

### Build Selection

**GRCh38** (default)
- Latest human genome reference
- Recommended for new analyses
- Most comprehensive cytoband annotation

**Legacy Builds** (hg19, hg18, hg17)
- Use for compatibility with older datasets
- Required when comparing with historical analyses
- May have less detailed cytoband information

```r
# Latest reference
result <- CytoConverter(data, build = "GRCh38")

# Legacy compatibility  
result <- CytoConverter(data, build = "hg19")
```

### Custom Cytoband Data

For specialized applications, custom cytoband reference data can be provided:
- Must match the expected format of built-in references
- Useful for non-human species or specialized coordinate systems
- Place custom files in Builds/ directory

## Recommended Configurations

### Clinical Constitutional Analysis
```r
result <- CytoConverter(data, 
                       constitutional = TRUE,
                       guess = FALSE,
                       sexstimate = TRUE,
                       build = "GRCh38")
```

### Cancer Research Analysis  
```r
result <- CytoConverter(data,
                       constitutional = FALSE, 
                       guess = TRUE,
                       count_fusions = TRUE,
                       build = "GRCh38")
```

### High-Throughput Screening
```r
result <- CytoConverter(data,
                       constitutional = FALSE,
                       guess = TRUE,
                       count_fusions = FALSE,  # Faster processing
                       build = "GRCh38")
```

### Permissive Research Analysis
```r
result <- CytoConverter(data,
                       constitutional = FALSE,
                       guess = TRUE,
                       guess_q = TRUE,
                       guess_by_first_val = TRUE,
                       allow_Shorthand = TRUE,
                       count_fusions = TRUE,
                       build = "GRCh38")
```

### Conservative Clinical Analysis
```r
result <- CytoConverter(data,
                       constitutional = TRUE,
                       guess = FALSE,
                       guess_q = FALSE,
                       sexstimate = TRUE,
                       count_fusions = FALSE,
                       build = "GRCh38")
```

## Performance Considerations

- **Fusion analysis** significantly increases processing time
- **Guessing parameters** can slow processing for complex karyotypes
- **Large datasets** benefit from command-line interface with threading
- **Memory usage** scales with input size and fusion complexity

## Error Handling

Always check the Error_log component of results:
```r
result <- CytoConverter(data, ...)
if (nrow(result$Error_log) > 0) {
    print("Warnings or errors encountered:")
    print(result$Error_log)
}
```

Common issues and solutions:
- **Parse errors**: Try enabling guess parameters
- **Unknown notation**: Check allow_Shorthand setting
- **Performance issues**: Disable fusion analysis for large datasets
- **Coordinate mismatches**: Verify correct build parameter