# BASIC INFO ---------------------------
##
## Script name:CreatFeatureMatrix.R
##
## Purpose of script: This is script for ready of feature matrix
##
## Author: JijunYu
##
## Date Created: 2021-12-11
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main -----------------------
##---------------library and load_all and parameter-----------------#####
library(tidyverse)
library(devtools)
library(data.table)
load_all()

setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
final_result <- readr::read_delim("../02_NeoSimData/output/StimNeoantigen/neoantigen_all_0115.tsv") %>%
  distinct()
set.seed(12)
final_result$judge <- plyr::mapvalues(x = final_result$type,
                                      from = c("neg","pos"),
                                      to = c(0,1))

final_result$judge <- as.factor(final_result$judge)
final_result <- final_result %>%
  rownames_to_column(var = "term_id")
final_result_v2 <- final_result %>%
     add_count(mutate_peptide) %>%
     group_by(mutate_peptide) %>%
     sample_n(size = 1) %>%
     ungroup()
#### 01 肽段特征计算 #####
# ManualPropertySelection <- c("aaComp","boman","crucianiProperties",
#                              "fasgaiVectors","instaIndex","kideraFactors",
#                              "protFP","stScales","tScales")
ManualPropertySelection <- c("aIndex","blosumIndices","boman","charge","crucianiProperties",
                      "fasgaiVectors","hmoment","hydrophobicity","instaIndex",
                      "kideraFactors","mswhimScores","pl","protFP","vhseScales","zScales")
cat("Calculating property of peptides...\n\n", sep = "")
mutate.property <- PropertyofPepSingle(peptides = final_result_v2$mutate_peptide,PropAll = ManualPropertySelection)

rownames(mutate.property) <- final_result_v2$mutate_peptide
mutate.property <- mutate.property %>%
  janitor::clean_names() %>%
  na.omit()

mutateML <- mutate.property %>%
  mutate(judge = as.factor(final_result_v2$judge))


#### 02 基于相关性的特征筛选 ####
# 2.1 未经过特征筛选的相关性展示 ###
output_path <- "./result_v3/02_mutate_wholepeptide_single_samerelation_diffimp_result/01_相关性计算/"
CorPlot(mutate.property,
        outputname =paste0(output_path,"AllCorrelation.pdf"),
        width = 24,
        height = 24,
        #order = "hclust",
        corvismethod = "number")
CorPlot(mutate.property,
        outputname =paste0(output_path,"AllCorrelation_hclust_circle.pdf"),
        width = 48,
        height = 48,
        order = "hclust",
        corvismethod = "circle")

# 2.2 基于特征之间的相关性删去特征
# 所有这些方法可以用两个函数实现。
parCor <- function(mat, nblocks=10, coreN=parallel:::detectCores(logical=F)){
  ## Add dummy variables if ncol(x) %% nblocks does not give a remainder of 0
  if (ncol(mat) %% nblocks != 0){
    DUMMY <- data.table::as.data.table(matrix(data=0, nrow=nrow(mat), ncol=nblocks-(ncol(mat) %% nblocks)))
    colnames(DUMMY) <- paste0("DUMMY_", 1:ncol(DUMMY))
    x <- cbind(data.table::as.data.table(mat), DUMMY)
  }else{
    x <- data.table::as.data.table(mat)
  }

  ## Split blocks
  NCOL <- ncol(x)
  SPLIT <- split(1:NCOL, rep(1:nblocks, each=NCOL/nblocks))
  COMBS <- data.table::CJ(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  ## Parallel computation
  cat("Starting parallelized computation of sub-matrices.\n")
  cl <- parallel::makeCluster(coreN, type="SOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=nrow(COMBS), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  results <- foreach::foreach(i=1:nrow(COMBS), .packages="data.table", .inorder=T, .options.snow=opts)%dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    cor(x[, ..G1], x[, ..G2])
  }
  cat("\n")
  close(pb)
  parallel::stopCluster(cl)
  gc();gc()

  ## Integrate sub-matrices into one correlation matrix
  corMAT <- matrix(0, ncol=NCOL, nrow=NCOL)
  for(i in 1:nrow(COMBS)){
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    COR <- results[[i]]
    corMAT[G1, G2] <- COR
  }
  close(pb)
  corMAT <- Matrix::forceSymmetric(corMAT, uplo="U")
  corMAT <- as.matrix(corMAT)
  colnames(corMAT) <- colnames(x)
  rownames(corMAT) <- colnames(x)
  P <- setdiff(1:ncol(x), grep("^DUMMY_", colnames(x)))
  corMAT <- corMAT[P,P]

  ## Replace diagonal NA to one
  for(i in 1:nrow(corMAT)){
    corMAT[i,i] <- 1
  }

  ## Finish
  rm(list=setdiff(ls(), c("corMAT")))
  gc();gc()
  return(corMAT)
}
mutate.property_Cor <- parCor(mutate.property)
mutate.property_Cor[which(is.na(mutate.property_Cor))] <- 0
mutate.property_Cor_filtered <- caret::findCorrelation(mutate.property_Cor,cutoff = 0.80,verbose = F,names = T)
mutate.property_filtered <- as.data.table(mutate.property)[,-mutate.property_Cor_filtered,with=F]


#### 03 基于机器学习的特征重要性计算 ####
cat("Feature selection by Variant important Fileter method...\n\n", sep = "")
mutate.property_filtered$judge <- final_result_v2$judge
FeatureSebyVariImpFil.result <- FeatureSebyVariImpFil(dataset = mutate.property_filtered)

#### 04 基于机器学习计算的重要性值，导出特征矩阵 ####
output_path <- "./result_v3/02_mutate_wholepeptide_single_samerelation_diffimp_result/02_特征值计算/"
ImportScore <- FeatureSebyVariImpFil.result[,1:2]
for(top_n_tmp in c(seq(5,41,5),41)){
  featureselect <- ImportScore %>%
    top_n(wt=classif.ranger.score,n = top_n_tmp) %>%
    pull(classif.ranger.feature)
  mutateML <- mutate.property_filtered %>%
    as.data.table() %>%
    .[,featureselect,with=F] %>%
    mutate(judge = final_result_v2$judge)
  save(mutateML,file = paste0(output_path,"mutateML_80_1160_",top_n_tmp,".RData"))
}


