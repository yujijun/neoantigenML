## ---------------------------
##
## Script name: PropertyofPepMutation.R
##
## Purpose of script:
##
## Author: Dr. JijunYu
##
## Date Created: 2021-12-22
##
## Copyright (c) Timothy Farewell, 2021
## Email: jijunyuedu@outlook.com
##
## ---------------------------
##
## Notes: 研究局部突变和全局突变对于免疫原性的影响
##
##
## ---------------------------

## set working directory for Mac and
## ---------------------------

# options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
# memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on macs.

## ---------------------------

# Description: This script mainly used to calculate property of peptides from mutation situation.
# Properties
# 单点突变  应该考虑突变对于周围几个位点的影响。 肽段背景不同，对预测精确度有显著的影响。
# 结构生物学可以做的分析，就是肽段也不变，然后对于单点进行泛突变，然后看看突变对于结合的影响。
# 体积大小 / 正负/芳香性质 单点突变考虑的性质有哪些
#
#### ------ mutation peptide fragment -------#####
PepFrag.mutation <- function(
  peptideSet = peptides,
  coreN=parallel::detectCores(logical=F),
  fragLenSet = 3:6, #make sure the length of peptide longer than 8
  CutStart = NULL, # a numeric of cut peptides from default is NULL
  Cutend = NULL
){
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(tidyverse)
  library(data.table)
  #start calculation
  set.seed(12)
  time.start <- proc.time()

  if(is.numeric(CutStart) & is.numeric(Cutend)){
    pepCut = sapply(1:length(peptideSet),function(i){stringr::str_sub(peptideSet[i],CutStart,Cutend)})
    peptideSet = pepCut
  }
  peptide.Frag.single <- function(peptide,fragLen){
    library(data.table)
    f <- sapply (1:max(nchar(peptide) - fragLen + 1),function(i){stringr::str_sub(peptide,i,i+fragLen-1)})
    start_pos <- 1:(nchar(peptide) - fragLen + 1)
    end_pos <- fragLen:nchar(peptide)
    df <- data.table(pep = peptide, f=f,
                     fragLen = fragLen,start_pos=start_pos,end_pos=end_pos)
    return(df)
  }

  pep.list <- list()
  n = 1
  for(fragLentmp in fragLenSet){
    print(paste0("Generate fragment with length == ",fragLentmp))
    parameterDT <- data.table::CJ(peptideSet,fragLentmp) %>%
      magrittr::set_colnames(c("Peptide","FragLen"))
    cl <- parallel::makeCluster(coreN,type = "PSOCK")
    doSNOW::registerDoSNOW(cl)
    sink(tempfile())
    pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
    sink()
    opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
    dt_peptdesc <- foreach::foreach(i=1:nrow(parameterDT), .inorder=T, .options.snow=opts)%dopar%{
      peptide.Frag.single(parameterDT$"Peptide"[[i]], parameterDT$"FragLen"[[i]])
    }
    dt_df <- dt_peptdesc %>% rbindlist()
    #%>%  cbind.data.frame() #Attention, this could not be unique
    close(pb)
    parallel::stopCluster(cl)
    gc();gc()
    pep.list[[n]] <- dt_df
    n = n +1
  }
  pep.df <- pep.list %>% rbindlist()
  # Finish the timer
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")

  # Clear logs
  rm(list=setdiff(ls(), c("pep.df")))
  return(pep.df)
}

#filtering all fragment by position of mutation
mutation.pos <- function(final_result){
  # first column and second column must be wild and mutation peptides
  wild <- strsplit(final_result[[1]],split = "")
  mutation <- strsplit(final_result[[2]],split = "")
  pos <-  c()
  for(i in 1:length(wild)){
    match_tmp <- wild[[i]] == mutation[[i]]
    pos_tmp <- which(match_tmp == FALSE)
    pos <- c(pos,pos_tmp)
  }
  final_result$pos <- pos
  return(final_result)
}



#### peptide description ####
# only choose peptide including mutation information.
#' Title peptide Feature Description
#' @export
#' @example
PepFeatureDesc.mutate <- function(
  peptideSet,
  pair_peptideSet,
  featureSet = NULL,
  coreN = parallel::detectCores(logical=F),
  fragLenSet = 2:8,
  CutStart = NULL, # a numeric of cut peptides from default is NULL
  Cutend = NULL,
  peptideSet.pos
){
  require(data.table)
  require(tidyverse)
  require(parallel)
  library(foreach)
  library(doSNOW)
  # Start calculation
  set.seed(12345)
  time.start <- proc.time()
  if(is.numeric(CutStart) & is.numeric(Cutend)){
    pepCut = sapply(1:length(peptideSet),function(i){stringr::str_sub(peptideSet[i],CutStart,Cutend)})
    peptideSet = pepCut
  }
  # Working functions
  peptideDescriptor.Batch <- function(peptide){
    unlist(list(
      Peptides::aIndex(peptide),
      Peptides::blosumIndices(peptide),
      Peptides::boman(peptide),
      Peptides::charge(peptide),
      Peptides::crucianiProperties(peptide),
      Peptides::fasgaiVectors(peptide),
      Peptides::hmoment(peptide),
      Peptides::hydrophobicity(peptide),
      #Peptides::instaIndex(peptide),
      Peptides::kideraFactors(peptide),
      Peptides::mswhimScores(peptide),
      Peptides::pI(peptide),
      Peptides::protFP(peptide),
      Peptides::vhseScales(peptide),
      Peptides::zScales(peptide)
    ))
  }
  peptideDescriptor.NameSet <- c("AliphaticIndex","BLOSUM1","BLOSUM2","BLOSUM3","BLOSUM4","BLOSUM5","BLOSUM6","BLOSUM7","BLOSUM8","BLOSUM9","BLOSUM10","Boman","Charge","PP1","PP2","PP3","F1","F2","F3","F4","F5","F6","HydrophobicMoment","Hydrophobicity","KF1","KF2","KF3","KF4","KF5","KF6","KF7","KF8","KF9","KF10","MSWHIM1","MSWHIM2","MSWHIM3","pI","ProtFP1","ProtFP2","ProtFP3","ProtFP4","ProtFP5","ProtFP6","ProtFP7","ProtFP8","VHSE1","VHSE2","VHSE3","VHSE4","VHSE5","VHSE6","VHSE7","VHSE8","Z1","Z2","Z3","Z4","Z5")
  peptideDescriptor.FragStat.Single <- function(peptide,pair_peptide,fragLen,peptideSet.pos){
    library(data.table)
    library(tidyverse)
    f <- sapply (1:max(nchar(peptide) - fragLen + 1),function(i){stringr::str_sub(peptide,i,i+fragLen-1)})
    start_pos <- 1:(nchar(peptide) - fragLen + 1)
    end_pos <- fragLen:nchar(peptide)
    df <- data.table(pep = peptide, f=f,
                     fragLen = fragLen,start_pos=start_pos,end_pos=end_pos)
    df <- df %>%
      filter(start_pos <= peptideSet.pos,end_pos >=peptideSet.pos)
    f <- df$f
    d <- sapply(f[nchar(f)==fragLen], function(s){peptideDescriptor.Batch(s)})
    data.table::data.table(
      "Peptide"=stringr::str_c(peptide,"_",pair_peptide),
      "FragLen"=fragLen,
      "AADescriptor"=peptideDescriptor.NameSet,
      "Min"=matrixStats::rowMins(d),
      "Max"=matrixStats::rowMaxs(d),
      "Mean"=matrixStats::rowMeans2(d),
      "Median"=matrixStats::rowMedians(d)
    )
  }
  # Parallelized calculation of descriptive statistics
  flag1 <- is.null(featureSet)
  flag2 <- length(grep("PeptDesc_", featureSet))>=1
  flag_comp <- flag1==T | (flag1==F & flag2==T)
  if(flag_comp==T){
    parameterDT <- data.table::CJ(peptideSet = stringr::str_c(peptideSet,"_",pair_peptideSet,"_",peptideSet.pos),fragLenSet) %>%
      magrittr::set_colnames(c("peptide","fragLen")) %>%
      tidyr::separate(col = peptide,into = c("peptide","pair_peptide","pos"),sep = "_")
    cl <- parallel::makeCluster(coreN, type="PSOCK")
    doSNOW::registerDoSNOW(cl)
    sink(tempfile())
    pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
    sink()
    opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
    dt_peptdesc <- foreach::foreach(i=1:nrow(parameterDT), .inorder=T, .options.snow=opts)%dopar%{
      peptideDescriptor.FragStat.Single(parameterDT$"peptide"[[i]],
                                        parameterDT$"pair_peptide"[[i]],
                                        parameterDT$"fragLen"[[i]],
                                        parameterDT$"pos"[[i]])
    } %>% data.table::rbindlist()
    close(pb)
    parallel::stopCluster(cl)
    gc();gc()
    col_id <- c("Peptide", "AADescriptor", "FragLen")
    col_val <- setdiff(colnames(dt_peptdesc), col_id)
    dt_peptdesc <- data.table::melt.data.table(dt_peptdesc, id=col_id, measure=col_val, variable.name="Stat", value.name="Value")
    dt_peptdesc[,"Feature":=paste0("PeptDesc_", AADescriptor, "_", Stat, "_", FragLen)][,"AADescriptor":=NULL][,"FragLen":=NULL][,"Stat":=NULL]
    dt_peptdesc <- data.table::dcast.data.table(dt_peptdesc, Peptide~Feature, value.var="Value", fun=mean)

  }else{
    dt_peptdesc <- data.table::data.table("Peptide"=peptideSet)
  }

  # Minimum set of features (optional)
  if(!is.null(featureSet)){
    featureSet <- intersect(colnames(dt_peptdesc), featureSet)
    dt_peptdesc <- dt_peptdesc[, c("Peptide", featureSet), with=F]
  }

  # Finish the timer
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")

  # Clear logs
  rm(list=setdiff(ls(), c("dt_peptdesc")))
  gc();gc()

  # Output
  return(dt_peptdesc)
}



# 首先考虑，蛋白质考虑突变影响的周围的几个长度 ？？？
# 突变位点
# 搞成K-mer
# 计算属性
# 计算比例差异 vs 属性的影响


