#!/usr/bin/env Rscript
# run_cytoconverter.R
# runs CytoConverter script with CYTO_SAMPLE.TXT as an input
# result files (result.txt and error.txt) are stored in tsv format
# result files need to be processed using process_cyto_data.sh to format the data before loading into BQ
#
# s ./input/cyto_input_* | awk -F_input_ '{print $2};' | xargs -n 1 -t -I {} -P 5 bash -c './main_sub.R -i ./input/cyto_input_{} -o ./output/output_{}.txt -l ./output_log_{}.txt > ./output/log_{}'


options(scipen = 999)


# Load modules

mod_cytoscript <- modules::use('cytoscript.R')


# Parse arguments

option_list = list(
    optparse::make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        help="input file name",
        metavar="character"
    ),
    optparse::make_option(
        c("-o", "--output"),
        type="character",
        default=NULL,
        help="output file name",
        metavar="character"
    ),
    optparse::make_option(
        c("-l", "--log"),
        type="character",
        default=NULL,
        help="log file name",
        metavar="character"
    )
)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$output)) {
    optparse::print_help(opt_parser)
    stop("Missing required argument: --output -> output file name")
}



# Call CytoConverter

#source('cytoscript.R')
input_m <- if (is.null(opt$input)) {
    # Read from stdin if no input specified
    as.matrix(as.data.frame(read.delim(file=file("stdin"), header=F, sep="\t")))
} else {
    as.matrix(as.data.frame(read.delim(file=opt$input, header=F, sep="\t")))
}

result <- mod_cytoscript$CytoConverter(input_m, constitutional=F, guess=T)

write.table(result[[1]], file=opt$output, quote=FALSE, sep='\t', col.names=F, row.names=F)
if (!is.null(opt$log)) {
    write.table(result[[2]], file=opt$log, quote=FALSE, sep='\t', col.names=F, row.names=F)
}

