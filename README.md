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

## Additional Information

Builds are at 850 resolution and provided for human genome builds GRCh38, hg19, hg18, and hg17
if wanted, the user can supply thier own list of cytobands to process as CytoConverter uses the 
bands at 850 resolution for build GRCh38 as default.

The function CytoConverter will output a list with the first element being the results table and 
the second element being the error and warning table. 

CytoConverter has two required paramenters:

- in_data - input karyotype or karyotype table
- build - a string of the build used for reference containing chromosome, chromosome position,
corresponding cytoband, and staining pattern (GRCh38, hg19, hg18, hg17). This parameter is set to
GRCh38 by default.

To access each of the elements, place the result into an R variable like so:

```
Variable_name <- CytoConverter(in_data);
```

To get the results table use ```Variable_name$Results```
To get the error log use ```Variable_name$Error_log```

The results table consists of the sample name followed by the clone line number, the start genomic
coordinate of a gain or loss, the end coordinate of a gain or loss, an indicator if the sample is a
gain or a loss, and the number of cells in a clone out of the total cells in a sample.

A built in function to create a graph displaying samples with gains and losses is provided:

- Source functions plot_cyto_graph and cyto_graph. 
- plot_cyto_graph contains 4 parameters:
  - cyto_list - table output from CytoConverter
  - list_from_cyto - output from cyto_graph (unnessesary if cyto_list is used)
  - ref_list - sets reference to use for plotting coordinates, default is GRCh38
  - ylabel - option to enable or disable printing sample names on the graph

