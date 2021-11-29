# BASIC INFO ---------------------------
##
## Script name: InstallPackage.R
##
## Purpose of script: This is package for package install
##
## Author: JijunYu
##
## Date Created: 2021-11-28
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse",
              "stringr",
              "devtools",
              "mlr3verse",
              "praznik",
              "kknn",
              "xgboost",
              "forestmodel",
              "rpart",
              "stats",
              "FSelectorRcpp",
              "care",
              "iml",
              "DALEX",
              "apcluster",
              "mlr3tuningspaces",
              "precrec") # list of packages to load
n_packages <- length(packages) # count how many packages are required
new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}
cat("ALL PACKAGES HAVE BEEN INSTALLED... \n\n", sep = "")
#--------------Main -----------------------
