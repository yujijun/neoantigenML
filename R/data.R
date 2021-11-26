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

#' positive and negative neopeptides from  netMHCpan
#' 
#' A dataset of HLA-A*02:01 related simulated neopeptides from netMHCpan result
#' @format A data frame with 755 rows and 7 variables in positive; A negative data frame with 6032 rows and  columns
#' 
"positive"
"negative"

#' formated neoantigens dataset for machine learning running
#' 
#' A dataset of HLA-A*02:01 related simulated neopeptides.
#' @format A data frame with 6787 rows and 4 rows.
#' 
"Neodataset"

