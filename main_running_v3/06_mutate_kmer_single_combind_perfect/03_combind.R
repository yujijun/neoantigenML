library(tidyverse)
library(devtools)
library(data.table)
load_all()
setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
output_path <- "./result_v3/06_mutate_kmer_single_combind_perfect_result/02_特征值计算/"
load(paste0(output_path,"mutate.single.property.withpeptide.RData"))
load(paste0(output_path,"mutate.kmer.withpeptide.RData"))
load(paste0(output_path,"peptideinfo.RData"))

combine_feature <- PepDesc %>%
  left_join(mutate.single.property.withpeptide,by = c("Peptide" = "peptide"))

colnames(peptideinfo) <- c("Peptide","Immunogenicity")

for(cutoff_tmp in seq(0,1,0.1)){
  featureDF_MHCI <- combine_feature
  minFeatureSet <- Features_MinimumFeatures(
    featureDFList=list(featureDF_MHCI),
    metadataDF=peptideinfo[,.(Peptide,Immunogenicity)][,Cluster:=.I],
    seedSet=1,
    corThreshold=cutoff_tmp,
    featureNSet=5000,
    criteria="intersect",
    returnImpDFList=T
  )
  combind_withjudge <- combine_feature %>%
    left_join(peptideinfo,by = "Peptide") %>%
    rename(judge = Immunogenicity)
  combind_select <- minFeatureSet$MinimumFeatureSet
  combind_filtered <- combind_withjudge[,c(combind_select,"judge"),with = F]
  cutoff_num <- cutoff_tmp * 100
  feature_num <- ncol(combind_filtered) - 1
  save(combind_filtered,file =  paste0(output_path,"combind_filtered_",cutoff_num,"_",feature_num,".RData"))
}


