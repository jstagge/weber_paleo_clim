# *------------------------------------------------------------------
# | PROGRAM NAME: 00-prepare_file_system
# | FILE NAME: 00-prepare_file_system.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This code installs all required packages to run the Paleoflow app.
# | 
# | 
# *------------------------------------------------------------------



###########################################################################
###  Check for necessary packages and install if needed
###########################################################################
### Set a list of packages
list_of_packages <- c("lubridate", "tidyverse", "ggplot2", "dataRetrieval", "devtools", "readxl", "scales", "svglite", "here", "data.table", "ggrepel", "cluster", "ggdendro", "ade4", "dendextend", "ggsci", "reshape2")

	
### Determine which packages are missing
package_list <- installed.packages()[,"Package"]
installed_test <- (list_of_packages %in% package_list)
packages_needed <- list_of_packages[!installed_test]

### If packages are missing, install them
if(length(packages_needed)) install.packages(packages_needed)

### A few packages have issues in the CRAN repository, must be installed from github
require(devtools)

devtools::install_github("jstagge/staggefuncs")
devtools::install_github("jstagge/paleoAPR")

# install_github('jstagge/paleoAPR', lib = "/home/jhstagge/R/x86_64-pc-linux-gnu-library/3.6")


