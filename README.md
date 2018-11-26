# CytoConverter
Change directory to file directory script is in before running script
Packages stringr and stringi are required for package. The script will attempt to install and load packages by default

Run script once to output function CytoConverter which accepts a string consisting of a karyotype or a table with one column denoting the sample name and the second column denoting the sample karyotype. Sample names must be different.

Supplementary file is cytoband.txt
if wanted, the user can supply thier own list of cytobands to process as CytoConverter uses the bands at 800 resolution for build GCHr38.


Function CytoConverter will output a list with the first element being the results table and the second element being the error and warning table 
To acess each of the elements, place the result into an R variable like so
##Variable_name<-CytoConverter(in_data);
To get results table use Variable_name$Results
To call the error log use Variable_name$Error_log

Results table consists of the sample name followed by the clone line number , the start genomic coordinate of a gain or loss, the end coordinate of a gain or loss, an indicator if the sample is a gain or a loss, and the number of cells in a clone out of total cells in a sample. 

