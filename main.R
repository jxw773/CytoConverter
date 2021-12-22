#!/usr/bin/env Rscript
# run_cytoconverter.R
# runs CytoConverter script with CYTO_SAMPLE.TXT as an input
# result files (result.txt and error.txt) are stored in tsv format
# result files need to be processed using process_cyto_data.sh to format the data before loading into BQ
#

options(scipen = 999)

library(modules)
mod_init <- modules::use('modules/init.R')
# load other libraries/requirements
mod_init$loadLibraries()

mod_cytobands <- modules::use('modules/cytobands.R')
mod_utils <- modules::use('modules/utils.R')
mod_merge <- modules::use('modules/merge.R')
mod_rowparser <- modules::use('modules/rowparser.R')
mod_parser <- modules::use('modules/parser.R')

source('cytoscript.R')
input_m<-as.matrix(as.data.frame(read.delim(file="input/CytoConverter_Example_File.txt", header=F, sep="\t")))
result<-CytoConverter(input_m,constitutional=F,guess=T)
write.table(result[[1]],file='output/CYTO_CONVERTED.TXT', quote=FALSE, sep='\t', col.names=F, row.names=F)
write.table(result[[2]], file='output/CYTO_CONVERTED_LOG.TXT', quote=FALSE, sep='\t', col.names=F, row.names=F)

