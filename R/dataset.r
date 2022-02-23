# This script is used for preprocessing datasets from data-raw folder.
# Start Time: Fri Nov 26 12:46:44 2021
# Update Time:
# Author:JijunYu

#### Preprocess of simulate datasets of neopeptides #####
library(tidyverse)
positive <- read.delim2("./data-raw/final.pos.neoantigens.filtered.tsv")
negative <-  read.delim2("./data-raw/final.neg.neoantigens.filtered.tsv")
positive.pep  <- positive[,1:3] %>%
  mutate(judge = 1)
negative.pep <- negative[,1:3] %>%
  mutate(judge = 0)
Neodataset <- rbind(positive.pep,negative.pep) %>%
  mutate(Length = str_length(wild_Peptide)) %>%
  filter(Length == 9) %>%
  distinct(wild_Peptide,.keep_all = TRUE)
usethis::use_data(Neodataset,overwrite = T)

#### example peptides ####
peptides <- Neodataset$wild_Peptide
usethis::use_data(peptides,overwrite = T)

#### ML example datasets ####
MLtestData <- read.csv("./data-raw/output.property.csv")
MLtestData <- MLtestData %>%
  mutate(judge = as.factor(Neodataset$judge)) %>%
  dplyr::select(!X)
MLtestData <- MLtestData[which(rowSums(is.na(MLtestData)) == 0),]
usethis::use_data(MLtestData,overwrite = T)

# #### IEDB human host and mus host datasets ####
# deredundancy.allhuman <- read.delim("./result/yjj/Tcell_deredundancy_all.human.clean.tsv")
# peptides.allhuman <- deredundancy.allhuman %>%
#   pull(description)
# deredundancy.allmouse <- read.delim("./result/yjj/Tcell_deredundancy_all.mouse.clean.tsv")
# peptides.allmouse <- deredundancy.allmouse %>%
#   pull(description)
# deredundancy.all <- read.delim("./result/yjj/Tcell_deredundancy_all.clean.tsv")
# peptides.allmouse <- deredundancy.allmouse %>%
#   pull(description)
# usethis::use_data(deredundancy.allhuman,deredundancy.allmouse,deredundancy.all,
#                   peptides.allhuman,peptides.allmouse,peptides.all,overwrite = T)
#

# #### Feature Matrix ####
# #use_data of calculated by CreatFeatureMatrix.R in ./main_running folder
# usethis::use_data(mousetestML,overwrite = TRUE)
