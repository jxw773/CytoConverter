#!/usr/bin/env Rscript

# install dependencies for CytoConverter

if (!require("modules")) {
    utils::install.packages("modules")
    if (!require("modules")) {
        stop("could not install modules")
    }
}

mod_install <- modules::use('modules/install.R')
mod_install$install_libraries()

