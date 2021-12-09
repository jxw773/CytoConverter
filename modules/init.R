# install/load libraries as needed

loadLibraries <- function() {
  if (require("stringr")) {
    print("stringr is loaded correctly")
  } else {
    print("trying to install stringr")
    install.packages("stringr")
    if (require("stringr")) {
      print("stringr installed and loaded")
    } else {
      stop("could not install stringr")
    }
  }
  
  if (require("DescTools")) {
    print("DescTools is loaded correctly")
  } else {
    print("trying to install DescTools")
    install.packages("DescTools")
    if (require("DescTools")) {
      print("DescTools installed and loaded")
    } else {
      stop("could not install DescTools")
    }
  }
  
  if (require("stringi")) {
    print("stringi is loaded correctly")
  } else {
    print("trying to install stringr")
    install.packages("stringi")
    if (require("stringi")) {
      print("stringi installed and loaded")
    } else {
      stop("could not install stringi")
    }
  }
  
  if (require("dplyr")) {
    print("dplyr is loaded correctly")
  } else {
    print("trying to install dplyr")
    install.packages("dplyr")
    if (require("dplyr")) {
      print("dplyr installed and loaded")
    } else {
      stop("could not install dplyr")
    }
  }
  
  if (require("hash")) {
    print("hash is loaded correctly")
  } else {
    print("trying to install hash")
    install.packages("hash", repos = 'https://mirrors.nics.utk.edu/cran')
    if (require("hash")) {
      print("hash installed and loaded")
    } else {
      stop("could not install hash")
    }
  }
}