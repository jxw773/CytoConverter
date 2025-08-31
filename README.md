# CytoConverter

[CytoConverter: a web-based tool to convert karyotypes to genomic coordinates](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3062-4)

Cytogenetic nomenclature is used to describe chromosomal aberrations (or lack thereof) in a collection of cells, referred to as the cells’ karyotype. The nomenclature identifies locations on chromosomes using a system of cytogenetic bands, each with a unique name and region on a chromosome. Each band is microscopically visible after staining, and encompasses a large portion of the chromosome. More modern analyses employ genomic coordinates, which precisely specify a chromosomal location according to its distance from the end of the chromosome. Currently, there is no tool to convert cytogenetic nomenclature into genomic coordinates. Since locations of genes and other genomic features are usually specified by genomic coordinates, a conversion tool will facilitate the identification of the features that are harbored in the regions of chromosomal gain and loss that are implied by a karyotype.

## Requirements

CytoConverter requires R 4.0+. Before running the main script, make sure that required R packages
are installed by changing to the CytoConverter directory and running:

```
./init.R
```

## Running CytoConverter

Run CytoConverter with the wrapper script using the following command:

```
./cytoconverter \
  --input input-file.txt \
  --threads 4 \
  --output output-file.txt \
  --log log-file.txt
```

Adjust parameters for your specific run:

- input: Input file of sample names and associated karyotypes, one per line, tab delimited.
- threads: Number of parallel threads to run. The input file will be split into pieces accordingly.
- output: Output file containing genomic coordinates and indications of gain or loss for all samples.
- log: Log file containing any warnings or errors encountered during processing.


## Code Structure

The code has been split into multiple "modules" and structured as follows.

- **cytoconverter.R**: includes the main entrypoint for CytoConverter. 
- **rowparser.R**: includes control flow for parsing rows of a karyotype table.
- **colparser.R**: includes control flow for parsing each cell line or component of a karyotype. 
- **merge.R**: includes helper functions for handling and merging intervals.
- **mergefusions.R**: includes specialized merge functions for fusion data and structural rearrangements.
- **gainlossfusion.R**: includes functions for detecting and classifying gains, losses, and fusion events.
- **utils.R**: includes utility functions.  
- **cytobands.R**: includes a function for getting cytobands for translocations and insertions.  

### Additional Files
- **cyto_graph_fusion.R**: specialized plotting functions for fusion visualization
- **plot_cyto_graph.R**: standard plotting functions for gains and losses
- **Fusion_tags_WIP.txt**: documentation of fusion tag system and supported fusion types  


```
   ┌────────────────────┐      ┌───────────────────────────┐
   │ cytoconverter.R    │      │ merge.R                   │
   │                    │      │                           │
   │    CytoConverter() │      │    insertSection()        │
   └───┬────────────────┘      │    deleteIntersections()  │
       │                       │    getContiguousSection() │
       │                       │    mergeAdjacentSections()│
       │                       │    mergeTable()           │
       │                       │    mergeDel()             │
   ┌───▼───────────┐           │    mergeDelmat()          │
   │ rowparser.R   ├───────────►    bigDelMerge()          │
   │               │           │    mergeDeletions()       │
   │    rowparse() │           └───────────────────────────┘
   └───┬────────┬──┘
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

   ┌─────────────────────────────────────────────────────────┐
   │ Fusion Detection & Analysis Modules                     │
   ├─────────────────────────────────────────────────────────┤
   │ gainlossfusion.R        │ mergefusions.R                │
   │                         │                               │
   │   gainloss()            │   insertSection_Fus()         │
   │   fusion()              │   deleteIntersections_Fus()   │
   │                         │   mergeAdjacentSections_Fus() │
   │                         │   mergeTable_Fus()            │
   └─────────────────────────────────────────────────────────┘

   ┌─────────────────────────────────────────────────────────┐
   │ Visualization Modules                                   │
   ├─────────────────────────────────────────────────────────┤
   │ cyto_graph_fusion.R     │ plot_cyto_graph.R             │
   │                         │                               │
   │   cyto_graph_fusion()   │   plot_cyto_graph()           │
   │                         │                               │
   └─────────────────────────────────────────────────────────┘
```

## Fusion Capabilities

CytoConverter has comprehensive support for detecting and analyzing chromosomal fusions and structural aberrations. The fusion detection module can identify and classify various types of chromosomal rearrangements, providing detailed information about their genomic coordinates and fusion characteristics.

### Supported Fusion Types

CytoConverter can detect and analyze the following types of chromosomal fusions and structural aberrations:

#### Basic Fusion Types
- **Translocations** - Both balanced and unbalanced translocations
  - `t(chromosome1;chromosome2)(breakpoint1;breakpoint2)` - Balanced translocations
  - `der(chromosome)t(chromosome1;chromosome2)(breakpoint1;breakpoint2)` - Derivative chromosomes from translocations
  
- **Inversions** - Chromosomal inversions
  - `inv(chromosome)(breakpoint1breakpoint2)` - Paracentric or pericentric inversions
  
- **Insertions** - Material inserted from one chromosome to another
  - `ins(chromosome1;chromosome2)(insertion_point;breakpoint1breakpoint2)` - Insertions

#### Complex Fusion Types
- **Derivative Chromosomes** - Chromosomes derived from structural rearrangements
  - `der(chromosome)` - General derivative chromosome notation
  - `rec(chromosome)` - Recombinant chromosomes
  
- **Ring Chromosomes** - Circular chromosomes formed by terminal deletions and fusion
  - `r(chromosome)(breakpoint1breakpoint2)` - Ring chromosomes
  
- **Multicentric Chromosomes** - Chromosomes with multiple centromeres
  - **Dicentric Chromosomes** - `dic(chromosome1;chromosome2)` - Two centromeres
  - **Tricentric Chromosomes** - `trc(chromosome1;chromosome2;chromosome3)` - Three centromeres
  
- **Isochromosomes** - Chromosomes with identical arms
  - `i(chromosomearm)` - Isochromosomes
  - `ider(chromosome)` - Isochromosome for derivative chromosome
  - `idic(chromosome)` - Isodicentric chromosomes
  
- **Robertsonian Translocations** - Fusion of acrocentric chromosomes
  - `rob(chromosome1;chromosome2)` - Robertsonian translocations

#### Additional Structural Aberrations
- **Duplications** - `dup(chromosome)(breakpoint1breakpoint2)` 
- **Triplications** - `trp(chromosome)(breakpoint1breakpoint2)`
- **Quadruplications** - `qdp(chromosome)(breakpoint1breakpoint2)`
- **Fragile Sites** - `fra(chromosome)(breakpoint)`
- **Centromere Fission** - `fis(chromosome)(breakpoint)`

### Fusion Input Format

Fusion karyotypes should follow standard cytogenetic nomenclature. Examples:

```
Sample1    47,XY,+der(10)t(10;21)(p13;q21)
Sample2    46,XX,der(1;19)(q10;p10)
Sample3    46,XY,t(9;22)(q34;q11.2)
Sample4    45,X,dic(X;Y)(p22.3;p11.3)
Sample5    46,XX,inv(16)(p13.1q22)
Sample6    47,XY,+r(7)(p22q36)
Sample7    46,XY,ins(2;5)(p13;q14q33)
```

### Fusion Output Format

When fusion analysis is enabled, CytoConverter outputs additional information about detected fusions:

- **Fusion Type Classification** - Each fusion is tagged with its specific type (e.g., `#translocation_balanced`, `#derivative_chrom`, `#ring_chrom`)
- **Chromosomal Breakpoints** - Precise genomic coordinates of fusion breakpoints
- **Fusion Partners** - Identification of chromosomes involved in each fusion
- **Fusion Orientation** - Direction and orientation of fused segments

### Using Fusion Analysis

To enable fusion detection and analysis:

```r
# Enable fusion counting and analysis
result <- CytoConverter(input_data, count_fusions = TRUE)

# Access fusion-specific results
fusion_table <- result$Results[grepl("#", result$Results$Type), ]

# Plot fusion data using fusion-specific plotting
plot_result <- cyto_graph_fusion(fusion_table, ref_list = "GRCh38")
```

### Fusion Visualization

CytoConverter provides specialized visualization capabilities for fusion data:

- **cyto_graph_fusion()** - Creates fusion-specific plots highlighting structural rearrangements
- **Fusion-specific color coding** - Different colors for different types of fusions
- **Breakpoint visualization** - Precise marking of fusion breakpoints on chromosome plots

## Additional Information

Builds are at 850 resolution and provided for human genome builds GRCh38, hg19, hg18, and hg17
if wanted, the user can supply thier own list of cytobands to process as CytoConverter uses the 
bands at 850 resolution for build GRCh38 as default.

The function CytoConverter will output a list with the first element being the results table and 
the second element being the error and warning table. 

CytoConverter has multiple parameters:

- **in_data** - input karyotype or karyotype table
- **build** - a string of the build used for reference containing chromosome, chromosome position,
corresponding cytoband, and staining pattern (GRCh38, hg19, hg18, hg17). This parameter is set to
GRCh38 by default.
- **count_fusions** - boolean flag to enable fusion detection and analysis (default: FALSE)
- **constitutional** - boolean flag for constitutional analysis (default: TRUE)
- **guess** - boolean flag to enable guessing of ambiguous regions (default: FALSE)

To access each of the elements, place the result into an R variable like so:

```
Variable_name <- CytoConverter(in_data);
```

To get the results table use ```Variable_name$Results```
To get the error log use ```Variable_name$Error_log```

The results table consists of the sample name followed by the clone line number, the start genomic
coordinate of a gain or loss, the end coordinate of a gain or loss, an indicator if the sample is a
gain, loss, or fusion type, and the number of cells in a clone out of the total cells in a sample.

For fusion analysis, additional columns may include fusion type classifications (marked with # symbols)
and breakpoint information for structural rearrangements.

### Visualization Functions

Built-in functions are provided to create graphs displaying samples with gains, losses, and fusions:

#### Standard Plotting
- **plot_cyto_graph()** - Standard plotting for gains and losses
- **cyto_graph()** - Data preparation for standard plotting

#### Fusion-Specific Plotting  
- **plot_cyto_graph_fusion()** - Specialized plotting for fusion data with enhanced visualization
- **cyto_graph_fusion()** - Data preparation for fusion plotting with additional fusion metadata

#### Parameters for Plotting Functions
- **cyto_list** - table output from CytoConverter
- **list_from_cyto** - output from cyto_graph (unnecessary if cyto_list is used)
- **ref_list** - sets reference to use for plotting coordinates, default is GRCh38
- **ylabel** - option to enable or disable printing sample names on the graph
- **include_normals_graph** - option to include normal samples in fusion graphs (fusion plotting only)

