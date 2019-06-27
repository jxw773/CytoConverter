# CytoConverter
Change directory to file directory script is in before running script
Packages stringr and stringi are required for package. The script will attempt to install and load packages by default

Run script once to output function CytoConverter which accepts a string consisting of a karyotype or a table with one column denoting the sample name and the second column denoting the sample karyotype. Sample names must be different.



Builds are at 850 resolution and provided for human genome builds GRCh38, hg19, hg18, and hg17
if wanted, the user can supply thier own list of cytobands to process as CytoConverter uses the bands at 850 resolution for build GRCh38 as default.


Function CytoConverter will output a list with the first element being the results table and the second element being the error and warning table 
CytoConverter has two paramenters
in_data - input karyotype or karyotype table
build - a string of the build used for reference containing chromosome, chromosome position, corresponding cytoband, and staining pattern (GRCh38, hg19, hg18, hg17). Parameter is set to GRCh38 by default.
To access each of the elements, place the result into an R variable like so
##Variable_name<-CytoConverter(in_data);
To get results table use Variable_name$Results
To call the error log use Variable_name$Error_log

Results table consists of the sample name followed by the clone line number , the start genomic coordinate of a gain or loss, the end coordinate of a gain or loss, an indicator if the sample is a gain or a loss, and the number of cells in a clone out of total cells in a sample. 

A built in function to create a graph displaying samples with gains and losses is provided. 
1) Source functions plot_cyto_graph and cyto_graph. 
2) plot_cyto_graph contains 4 parameters 
cyto_list - table output from CytoConverter
list_from_cyto - output from cyto_graph (unnessesary if cyto_list is used)
ref_list - sets reference to use for plotting coordinates, default is GRCh38
ylabel - option to enable or disable printing sample names on the graph
