#!/usr/bin/env Rscript
# run_cytoconverter.R
# runs CytoConverter script with CYTO_SAMPLE.TXT as an input
# result files (result.txt and error.txt) are stored in tsv format
# result files need to be processed using process_cyto_data.sh to format the data before loading into BQ
#

options(scipen = 999)

library(modules)
mod <- modules::use('modules')

# load other libraries/requirements
mod$init$loadLibraries()

source('cytoscript_vinput.R')
input_m<-as.matrix(as.data.frame(read.delim(file="input/CytoConverter_Example_File.txt", header=F, sep="\t")))
result<-CytoConverter(input_m,constitutional=F,guess=T)
write.table(result[[1]],file='output/CYTO_CONVERTED.TXT', quote=FALSE, sep='\t', col.names=F, row.names=F)
write.table(result[[2]], file='output/CYTO_CONVERTED_LOG.TXT', quote=FALSE, sep='\t', col.names=F, row.names=F)

