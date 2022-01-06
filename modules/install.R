# install/load libraries as needed

install_libraries <- function() {
    if (!require("stringr")) {
        utils::install.packages("stringr")
        if (!require("stringr")) {
            stop("could not install stringr")
        }
    }
  
    if (!require("DescTools")) {
        utils::install.packages("DescTools")
        if (!require("DescTools")) {
            stop("could not install DescTools")
        }
    }
  
    if (!require("stringi")) {
        utils::install.packages("stringi")
        if (!require("stringi")) {
            stop("could not install stringi")
        }
    }
  
    if (!require("dplyr")) {
        utils::install.packages("dplyr")
        if (!require("dplyr")) {
            stop("could not install dplyr")
        }
    }
  
    if (!require("hash")) {
        utils::install.packages("hash", repos = 'https://mirrors.nics.utk.edu/cran')
        if (!require("hash")) {
            stop("could not install hash")
        }
    }

    if (!require("optparse")) {
        utils::install.packages("optparse")
        if (!require("optparse")) {
            stop("coud not install optparse")
        }
    }
}

