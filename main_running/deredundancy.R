# BASIC INFO ---------------------------
##
## Script name: deredundance.R
##
## Purpose of script: This is script for dataset deredundancy
##
## Author: JijunYu
##
## Date Created: 2021-12-04
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main -----------------------

#### data input and clean colnames ####
library(tidyverse)
library(Biostrings)
library(motifStack)
library(ggpubr)
input_folder <- "/Users/yujijun/Documents/work/4_AntigenML/NeoantigenML/input/"
output_folder <- "/Users/yujijun/Documents/work/4_AntigenML/NeoantigenML/result/yjj/"
Tcell_v3 <- read.delim(stringr::str_c(input_folder,"tcell_full_v3.csv"),sep = ",",header = T)
Tcell_v3_clean <- Tcell_v3 %>%
  select(matches("Epitope")|matches("Host")|matches("Assay")|matches("MHC"))
colnames(Tcell_v3_clean) <- Tcell_v3_clean[1,]
Tcell_v3_clean <- Tcell_v3_clean[-1,]
Tcell_v3_clean <- janitor::clean_names(Tcell_v3_clean)

####delete non-fit peptides ####
#1. delete NA rows
Tcell_v3_clean <- Tcell_v3_clean %>%
  filter(!if_any(.cols =  c(description,qualitative_measure,allele_name),.fns = `==` , "")) %>% mutate(pep_length = stringr::str_length(description))
#312710

#2. delete non Linear peptide#
Tcell_v3_clean <- Tcell_v3_clean %>%
  arrange(description) %>%
  filter(object_type == "Linear peptide") %>%
  arrange(desc(pep_length))
#3.delete peptide with [+] in them
Tcell_v3_clean <- Tcell_v3_clean[-c(grep("[+]",Tcell_v3_clean$description)),]
#302480

#### filter peptides ####
#1. filter replication columns
Tcell_v3_clean <- Tcell_v3_clean %>%
  mutate(pep_length == str_length(description))
  filter(pep_length == 9) %>%
  distinct(description,allele_name,
               qualitative_measure,
               organism_name,name,pep_length)
test <- Tcell_v3_clean
peptides.single <- names(table(test$description)[table(test$description) == 1])
peptides.rep <- names(table(test$description)[table(test$description) > 1])
Tcell_deredundancy_single <- Tcell_v3_clean %>%
  filter(description %in% peptides.single)
Tcell_deredundancy_rep <- Tcell_v3_clean %>%
  filter(description %in% peptides.rep)
# 肽段去冗余过程如下 #
#如果有一个为阳性，则将该阳性结果留下；
#如果超过一个阳性结果，则判断mhc 是否相同，取所有MHC亚型留下
#如果没有阳性，则将所有MHC亚型留下。
df.empty <- Tcell_deredundancy_single[1,]
for(peptide in peptides.rep){
  df.tmp <- Tcell_deredundancy_rep %>%
    filter(description == peptide)
  if(is_empty(grep("Positive",df.tmp$qualitative_measure)) == FALSE){
    df.tmp <- df.tmp[grep("Positive",df.tmp$qualitative_measure),]
    df.tmp <- df.tmp %>%
      distinct(allele_name,.keep_all =T) %>%
      filter(!(allele_name == ""))
  }else{
    df.tmp <- df.tmp %>%
      distinct(allele_name,.keep_all=T) %>%
      filter(!(allele_name == ""))
  }
  df.empty <- rbind(df.empty,df.tmp)
}
df.empty <- df.empty[-1,]
Tcell_deredundancy_all <- rbind(Tcell_deredundancy_single,df.empty)
Tcell_deredundancy_all$name[str_detect(Tcell_deredundancy_all$name,pattern = "Homo sapiens")]<- "Homo sapiens"
# Homo sapiens
# 30
# Homo sapiens (human)
# 24398
# Homo sapiens Black
# 29
# Homo sapiens Caucasian
Tcell_deredundancy_all.judge <- Tcell_deredundancy_all %>%
  mutate(judge = if_else(qualitative_measure == "Negative",0,1)) %>%
  group_by(description) %>%
  mutate(judge.sum = sum(judge)) %>%
  ungroup() %>%
  mutate(judge.final = if_else(judge.sum == 0,0,1))

#### final #####
# host is human
Tcell_deredundancy_all.hosthuman = Tcell_deredundancy_all.judge %>%
  filter(!is.na(str_match(pattern = "Homo sapiens",string=name))) %>%
  distinct(description,.keep_all = T)

Tcell_deredundancy_all.hostMus = Tcell_deredundancy_all.judge %>%
  filter(!is.na(str_match(pattern = "Mus musculus",string=name))) %>%
  distinct(description,.keep_all = T)

Tcell_deredundancy_all.distinct = Tcell_deredundancy_all.judge %>%
  distinct(description,.keep_all = T)

save(Tcell_deredundancy_all.distinct, file = "./result/yjj/Tcell_deredundancy_all.distinct.RData")
save(Tcell_deredundancy_all.hosthuman, file = "./result/yjj/Tcell_deredundancy_hosthuman.RData")
save(Tcell_deredundancy_all.hostMus, file = "./result/yjj/Tcell_deredundancy_hostMus.RData")



host_table <- as.data.frame(table(Tcell_deredundancy_all.distinct$name))



