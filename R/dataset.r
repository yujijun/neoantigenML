# This script is used for preprocessing datasets from data-raw folder.
# Start Time: Fri Nov 26 12:46:44 2021
# Update Time:
# Author:JijunYu

#### Preprocess of simulate datasets of neopeptides #####
positive <- read.delim2("./data-raw/final.pos.neoantigens.filtered.tsv")
negative <-  read.delim2("./data-raw/final.neg.neoantigens.filtered.tsv")
positive.pep  <- positive %>%
  select(c(1,2,3)) %>%
  mutate(judge = 1)
negative.pep <- negative %>%
  select(c(1,2,3)) %>%
  mutate(judge = 0)
Neodataset <- rbind(positive.pep,negative.pep)
usethis::use_data(Neodataset,overwrite = T)
