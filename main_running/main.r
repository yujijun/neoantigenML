# BASIC INFO ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: JijunYu
##
## Date Created: 2021-11-27
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##

# SET WORKING DIRECTORY -----------------------------
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "/cloud/project"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
output_path <- "./result/yjj/"
## SET OPTIONS ---------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse", "stringr", "readxl") # list of packages to load
n_packages <- length(packages) # count how many packages are required
new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}
# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}

#--------------Main -----------------------
load_all()
#### Original dataset input ####
#### Property calculation and selected manually ####
ManualPropertySelection <- c("aaComp","boman","crucianiProperties",
                             "fasgaiVectors","instaIndex","kideraFactors",
                             "protFP","stScales","tScales")
output.property <- PropertyofPepSingle(peptides = peptides,PropAll = ManualPropertySelection)
#### Create machine learning dataset ####
MLtestData <- output.property %>%
  mutate(judge = as.factor(Neodataset$judge)) ## Why some is NA
MLtestData <- MLtestData[which(rowSums(is.na(MLtestData)) == 0),]
#### Feature Selection ####
FeatureSeByFilter.result <- FeatureSeByFilter(dataset = MLtestData)
FeatureSebyVariImpFil.result <- FeatureSebyVariImpFil(dataset = MLtestData)
#### Training ####
Training.rpart <- neoML.rpart(MLtestData)
#### Model Tuning ####
#### Performance Evaluation and Comparison by benmark ####
#### Visualization ####
#### output ####
write.csv(output.property,file = paste0(output_path,
                                          "output.property.csv"),
          quote = F)


