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


# #### ---- 5.0 思考要做测试的数据集 ---- ####
load("./data/deredundancy.allclean.rda")
load("./data/deredundancy.allhuman.rda")
load("./data/deredundancy.allmouse.rda")
load("./data/peptides.allclean.rda")
load("./data/peptides.allhuman.rda")
load("./data/peptides.allmouse.rda")


#### ---- 5.5 包含突变前后的所有特征集 --------####
setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
output_path <- "./result_v3/06_mutate_kmer_single_combind_perfect_result/02_特征值计算/"
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
save(PepDesc,file = paste0(output_path,"mutate.kmer.withpeptide.RData"))
save(peptideinfo,file = paste0(output_path,"peptideinfo.RData"))


colnames(peptideinfo) <- c("Peptide","Immunogenicity")
minFeatureSet <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCI),
  metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1,
  corThreshold=0.75,
  featureNSet=500,
  criteria="intersect",
  returnImpDFList=T
)
PepDesc_withjudge <- PepDesc %>%
  left_join(peptideinfo,by = "Peptide") %>%
  rename(judge = Immunogenicity)
PepDesc_select <- minFeatureSet$MinimumFeatureSet
PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered75_90.RData"))


minFeatureSet <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCI),
  metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1,
  corThreshold=0.85,
  featureNSet=500,
  criteria="intersect",
  returnImpDFList=T
)
PepDesc_withjudge <- PepDesc %>%
  left_join(peptideinfo,by = "Peptide") %>%
  rename(judge = Immunogenicity)
PepDesc_select <- minFeatureSet$MinimumFeatureSet
PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered85_146.RData"))

minFeatureSet <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCI),
  metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1,
  corThreshold=0.95,
  featureNSet=500,
  criteria="intersect",
  returnImpDFList=T
)
PepDesc_withjudge <- PepDesc %>%
  left_join(peptideinfo,by = "Peptide") %>%
  rename(judge = Immunogenicity)
PepDesc_select <- minFeatureSet$MinimumFeatureSet
PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered95_459.RData"))
