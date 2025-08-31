#' Package Installation Module for CytoConverter
#' 
#' @description
#' This module handles the automated installation and loading of required R packages
#' for CytoConverter. It provides a robust installation system that handles missing
#' packages gracefully and includes fallback CRAN mirrors for improved reliability.
#' 
#' @details
#' The installation system handles the following required dependencies:
#' \itemize{
#'   \item **stringr** - String manipulation and pattern matching
#'   \item **stringi** - Advanced string processing (required by stringr)
#'   \item **DescTools** - Statistical tools and descriptive analysis functions
#'   \item **dplyr** - Data manipulation and transformation
#'   \item **hash** - Hash table data structures for efficient merging
#'   \item **optparse** - Command-line argument parsing
#' }
#' 
#' Features include:
#' \itemize{
#'   \item Automatic detection of missing packages
#'   \item Graceful handling of installation failures
#'   \item Alternative CRAN mirror for hash package
#'   \item Clear error messages for troubleshooting
#' }

# install/load libraries as needed

#' Install Required Libraries for CytoConverter
#' 
#' @description
#' This function automatically installs and loads all required R packages for
#' CytoConverter. It checks for package availability and attempts installation
#' if packages are missing, with appropriate error handling for installation failures.
#' 
#' @return NULL (called for side effects of package installation/loading)
#' 
#' @details
#' The function implements a robust installation strategy:
#' 
#' **For each required package**:
#' \enumerate{
#'   \item Check if package is already installed and loadable
#'   \item If not available, attempt installation from CRAN
#'   \item Verify successful installation by attempting to load
#'   \item Provide clear error message if installation fails
#' }
#' 
#' **Special cases**:
#' \itemize{
#'   \item **hash package**: Uses alternative CRAN mirror (mirrors.nics.utk.edu)
#'     for improved reliability in some network environments
#'   \item **stringi dependency**: Required by stringr, installed automatically
#' }
#' 
#' **Error handling**: If any package fails to install, the function stops
#' with a descriptive error message to help users troubleshoot installation issues.
#' 
#' @examples
#' \dontrun{
#' # Install all required packages
#' install_libraries()
#' }
#' 
#' @note This function requires internet connectivity and appropriate permissions
#' to install packages in the R library directory. Use sudo if permission errors occur.
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

