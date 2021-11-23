#### create  posivive and negative test dataset
library(tidyverse)
positive <- read.delim2("./data-raw/final.pos.neoantigens.filtered.tsv")
negative <-  read.delim2("./data-raw/final.neg.neoantigens.filtered.tsv")
positive.pep  <- positive %>% 
  select(c(1,2,3)) %>% 
  mutate(judge = 1)
negative.pep <- negative %>% 
  select(c(1,2,3)) %>% 
  mutate(judge = 0)
Neodataset <- rbind(positive.pep,negative.pep)
usethis::use_data(negative,positive,Neodataset,overwrite = TRUE)
