#Add packages to the following arrays, then run the code below
cran.pkgs <-
  c("DBI", 
    "here", 
    "scales", 
    "caret", 
    "janitor", 
    "fitdistrplus", 
    "ggthemes", 
    "WRS2", 
    "tidyverse" 
  )

bioc.pkgs <-
  c(
  )

install.packages(cran.pkgs,dependencies=TRUE, repos="http://cran.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite(bioc.pkgs)
