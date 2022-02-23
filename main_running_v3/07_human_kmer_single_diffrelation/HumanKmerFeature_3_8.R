
########--------------library and load dataset------------------------------- ####
library(devtools)
library(tidyverse)
library(data.table)
library(Repitope)
library(parallel)
library(foreach)
library(doSNOW)
#options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
#load_all()
setwd("/public/home/yujijun01/01_NeoantigenML/")

# ------------------------------ action ------------------------------------------
# #############################################################################
# ########---- Prepare and generate all datasets for ML method ------------####
# #############################################################################
# #### ---- 5.0 思考要做测试的数据集 ---- ####
load("./data/deredundancy.allclean.rda")
load("./data/deredundancy.allhuman.rda")
load("./data/deredundancy.allmouse.rda")
load("./data/peptides.allclean.rda")
load("./data/peptides.allhuman.rda")
load("./data/peptides.allmouse.rda")

#### ---- 1 计算生成的所有的不包含CPP特征的matrix----####
# generate parameter #
peptide = peptides.allhuman
fragLen = 3:8
peptideinfo = deredundancy.allhuman
output_path = "result_v3/07_human_kmer_single_diffrelation_result/02_特征值计算/"
#output_prefix = "Human"


peptideinfo$Immunogenicity <- plyr::mapvalues(x = peptideinfo[,3],
                                              from = c(1,0),
                                              to = c("Positive","Negative"))
colnames(peptideinfo) <- c("Peptide","allele_name","judge","Immunogenicity")
peptideinfo$Immunogenicity <- factor(peptideinfo$Immunogenicity,levels = c("Positive","Negative"))
peptideinfo <- as.data.table(peptideinfo)

# 1.1  属性计算 -----------------####
PepDesc <- Features_PeptDesc(peptideSet = peptide,
                             fragLenSet = fragLen)
PepDesc <- PepDesc[,-c(2)]
#1.2 属性筛选 - 使用作者的方法进行相关性 和 重要性筛选 #####
PepDesc <- PepDesc %>%
  na.omit()
print(dim(PepDesc))
featureDF_MHCI <- PepDesc

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
    left_join(peptideinfo[,c(1,3)],by = "Peptide")
  PepDesc_select <- minFeatureSet$MinimumFeatureSet
  PepDesc_filtered <- PepDesc_withjudge[,c(PepDesc_select,"judge"),with = F]
  cutoff_num <- cutoff_tmp * 100
  feature_num <- ncol(PepDesc_filtered) - 1
  save(PepDesc_filtered,file =  paste0(output_path,"PepDesc_filtered","_",cutoff_num,"_",feature_num,".RData"))
}

for(cutoff_tmp in seq(0,1,0.1)){
  Filter_fun(cutoff_tmp = cutoff_tmp)
}


