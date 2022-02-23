########--------------library and load dataset------------------------------- ####
library(devtools)
library(tidyverse)
library(data.table)
library(Repitope)
library(parallel)
library(foreach)
library(doSNOW)
### ---- 5.5 包含突变前后的所有特征集 --------####
setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
load_all()
output_path <- "./result_v3/04_mutate_kmer_single_differelation_result/02_特征值计算/"
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

Filter_fun <- function(cutoff_tmp){
  minFeatureSet <- Features_MinimumFeatures(
    featureDFList=list(featureDF_MHCI),
    metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
    seedSet=1,
    corThreshold=cutoff_tmp,
    featureNSet=5000,
    criteria="intersect",
    returnImpDFList=T
  )
  PepDesc_withjudge <- PepDesc %>%
    left_join(peptideinfo,by = "Peptide") %>%
    rename(judge = Immunogenicity)
  PepDesc_select <- minFeatureSet$MinimumFeatureSet
  PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
  cutoff_num <- cutoff_tmp * 100
  feature_num <- ncol(PepDesc_filtered) - 1
  save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered","_",cutoff_num,"_",feature_num,".RData"))
}

for(cutoff_tmp in seq(0,1,0.1)){
  Filter_fun(cutoff_tmp = cutoff_tmp)
}
