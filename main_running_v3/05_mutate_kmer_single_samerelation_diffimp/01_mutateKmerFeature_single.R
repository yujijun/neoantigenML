#######--------------library and load dataset------------------------------- ####
library(devtools)
library(tidyverse)
library(data.table)
library(Repitope)
library(parallel)
library(foreach)
library(doSNOW)
options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
load_all()
setwd("/public/home/yujijun01/01_NeoantigenML/")

#### ---- 5.5 包含突变前后的所有特征集 --------####
setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
output_path <- "./result_v3/05_mutate_kmer_single_samerelation_diffimp_result/02_特征值计算/v_3/"
final_result <- readr::read_delim("../02_NeoSimData/output/StimNeoantigen/neoantigen_all_0115.tsv") %>%
  distinct()
set.seed(12)
final_result$judge <- plyr::mapvalues(x = final_result$type,
                                      from = c("neg","pos"),
                                      to = c(0,1))
final_result$judge <- as.factor(final_result$judge)
final_result <- final_result %>%
  rownames_to_column(var = "term_id") %>%
  relocate(term_id,.before = 1)
final_result_v2 <- final_result %>%
  add_count(mutate_peptide) %>%
  group_by(mutate_peptide) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  as.data.table()
final_result_v2 <- final_result_v2  %>%
  select(-c("term_id","n"))

peptide <- final_result_v2$mutate_peptide
fragLen <- 3:8
PepDesc <- Features_PeptDesc(peptideSet = peptide,
                             fragLenSet = fragLen)
PepDesc <- PepDesc[,-c(2)]
#1.2 属性筛选 - 使用作者的方法进行相关性 和 重要性筛选 #####
PepDesc <- PepDesc %>%
  na.omit()
print(dim(PepDesc))
featureDF_MHCI <- PepDesc
peptideinfo <- final_result_v2[,c("mutate_peptide","judge")]
colnames(peptideinfo) <- c("Peptide","Immunogenicity")

minFeatureSet <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCI),
  metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1,
  corThreshold=0.8,
  featureNSet = 5000,
  #featureNSet=feature_tmp,
  criteria="intersect",
  returnImpDFList=T
)
PepDesc_withjudge <- PepDesc %>%
  left_join(peptideinfo,by = "Peptide") %>%
  rename(judge = Immunogenicity)
Import_matrix <- minFeatureSet$ImpDFList_Asc$`.1`$imp[,c(1,2)]

for(feature_tmp in c(seq(10,50,2))){
  PepDesc_select <- Import_matrix %>%
    top_n(n = feature_tmp,wt = Importance) %>%
    pull(FeatureID)
  PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
  feature_num <- ncol(PepDesc_filtered) -1
  save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered_80","_",feature_num,".RData"))
}


