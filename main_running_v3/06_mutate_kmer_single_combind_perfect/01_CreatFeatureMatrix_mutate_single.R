
##---------------library and load_all and parameter-----------------#####
library(tidyverse)
library(devtools)
library(data.table)
load_all()
set.seed(12)
setwd("E:/Projects/01_Neoantigen/01_NeoantigenML_v2/")
output_path <- "./result_v3/06_mutate_kmer_single_combind_perfect_result/02_特征值计算/"
final_result <- readr::read_delim("../02_NeoSimData/output/StimNeoantigen/neoantigen_all_0115.tsv") %>%
  distinct()
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
#save(final_result_v2,file = paste0(output_path,"final_result_v2_mutate_quchong.RData"))

#### 01 肽段特征计算 #####
ManualPropertySelection <- c("aIndex","blosumIndices","boman","charge","crucianiProperties",
                      "fasgaiVectors","hmoment","hydrophobicity","instaIndex",
                      "kideraFactors","mswhimScores","pl","protFP","vhseScales","zScales")
cat("Calculating property of peptides...\n\n", sep = "")
mutate.single.property <- PropertyofPepSingle(peptides = final_result_v2$mutate_peptide,PropAll = ManualPropertySelection)

rownames(mutate.single.property) <- final_result_v2$mutate_peptide
mutate.single.property <- mutate.single.property %>%
  janitor::clean_names() %>%
  na.omit()

mutate.single.property.withpeptide <- mutate.single.property %>%
  mutate(peptide = final_result_v2$mutate_peptide) %>%
  relocate(peptide,.before = 1)

mutate.single.property.withjudge <- mutate.single.property.withpeptide %>%
  left_join(final_result_v2[,c("mutate_peptide","judge")],by = c("peptide"="mutate_peptide"))
save(mutate.single.property.withpeptide,file = paste0(output_path,"mutate.single.property.withpeptide.RData"))
save(mutate.single.property.withjudge,file = paste0(output_path,"mutate.single.property.withjudge.RData"))

# 2.2 基于特征之间的相关性删去特征
CorFiltered <- function(propertyMatrix = mutate.single.property.withpeptide, cutoff_tmp = 0.75){
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
  propertyMatrix_Cor <- parCor(propertyMatrix[,!(colnames(propertyMatrix) %in% "peptide")])
  propertyMatrix_Cor[which(is.na(propertyMatrix_Cor))] <- 0
  propertyMatrix_Cor_filtered <- caret::findCorrelation(propertyMatrix_Cor,cutoff = cutoff_tmp,verbose = F,names = T)
  propertyMatrix_filtered <- as.data.table(propertyMatrix)[,-propertyMatrix_Cor_filtered,with=F]
  return(propertyMatrix_filtered)
}
mutate.single.property_filtered75 <- CorFiltered(propertyMatrix = mutate.single.property.withpeptide,0.75)

