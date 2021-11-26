# This is main function for whole running.
# Author: jijunyu
# email: jijunyuedu@outlook.com
# Update Time: Fri Nov 26 12:23:06 2021

#### Library and hyperparameter ####
library(devtools)
library(tidyverse)
library(mlr3verse)
load_all()
output_path <- "./result/yjj/"
#### Original dataset input ####
#### Property calculation ####
output.property <- PropertyofPepSingle(peptides = peptides)
#### Feature Selection ####
#### Training ####
#### Performance Evaluation and Comparison ####
#### output ####
write.csv(output.property,file = paste0(output_path,
                                          "output.property.csv"),
          quote = F)
