## ---------------------------
##
## Script name: pepKmerFeature.R
##
## Purpose of script:
##
## Author: Dr. JijunYu
##
## Date Created: 2021-12-20
##
## Copyright (c) Timothy Farewell, 2021
## Email: jijunyuedu@outlook.com
## ---------------------------
##
## Notes: This is a function group for all kmer feature generation
##
##
## ---------------------------


## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
library(devtools)
library(parallel)
library(snow)
library(foreach)
#### Create Fake code of peptide sets ####
#' Fragment feature calculation
#'
#' \code{pepFakecode} Rscript Coding peptide with specific rule of code.
#' \code{pepFrag} Rscript to generate Fragment of peptides.
#' \Code{pepFragSign} filtering all significant fragment
#' \code{PepFeatureDesc} pepFeature description
#'
#' @param peptideSet A vector a group of peptides
#' @param coreN The number of cores, default is all core in your machine
#' @param peptides A vector a group of peptides
#' @param JudgeMatrix A data.table with immune and non-immne info.
#' @param prop A proportion of fragment in all peptides/
#' @param coreN The number of cores, default is all core in your machine
#' @param fragLenSet A set of sliding window
#' @param featureSet  A vector to select feature,default is NULL
#' @param Fake a logical TRUE or FALSE
#' @param FakeMatrix A dataframe of fakacode info, if Fake is ture,then this parameter need to be gaven. Colnames
#' @param CutStart # a numeric of cut peptides from, default is NULL
#' @param Cutend # a number of cut peptides to, default is NULL
#' @export
#' @name  Feature_Calcualtion
#' @examples
#' peptideSetFake <- PepFakeCode(peptideSet)

PepFakeCode <- function(peptideSet,
                        coreN = parallel::detectCores(logical=F)){
  from <- c("A","F","I","L","M","P","V","W","C","D","E","G","H","K","N","Q","R","S","T","Y")
  to <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
  PepFakeCode.single <- function(peptide){
    pepSplit <- stringr::str_split(peptide,pattern = "")[[1]]
    pepSplitCov <- plyr::mapvalues(x = pepSplit,
                                   from = from,
                                   to = to)
    pepSplitCovPaste <- paste0(pepSplitCov,collapse ="")
    return(pepSplitCovPaste)
  }
  cl <- parallel::makeCluster(coreN,type = "PSOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=length(peptideSet), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  peptideSetFake <- foreach::foreach(i=1:length(peptideSet), .inorder=F, .options.snow=opts)%dopar%{
    PepFakeCode.single(peptideSet[i])
  } %>% unlist()  #Attention, this could not be unique
  close(pb)
  parallel::stopCluster(cl)
  gc()
  # Clear logs
  rm(list=setdiff(ls(), c("peptideSetFake")))
  return(peptideSetFake)
}



#### Create all fragment of peptide sets  ####
#' @title Create fragment of peptide set
#' @return A vector with all fragment of peptide set
#' @export
#' @examples
#' PepFragResult <- PepFrag(peptideSet = peptides.allmouse,fragLenSet = 3:8)

PepFrag <- function(
  peptideSet = peptides,
  coreN=parallel::detectCores(logical=F),
  fragLenSet = 3:8, #make sure the length of peptide longer than 8
  CutStart = NULL, # a numeric of cut peptides from default is NULL
  Cutend = NULL
){
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(tidyverse)
  #start calculation
  set.seed(12)
  time.start <- proc.time()

  if(is.numeric(CutStart) & is.numeric(Cutend)){
    pepCut = sapply(1:length(peptideSet),function(i){stringr::str_sub(peptideSet[i],CutStart,Cutend)})
    peptideSet = pepCut
  }
  peptide.Frag.single <- function(peptide,fragLen){
    f <- sapply (1:max(nchar(peptide) - fragLen + 1),function(i){stringr::str_sub(peptide,i,i+fragLen-1)})
    return(f)
  }
  parameterDT <- data.table::CJ(peptideSet,fragLenSet) %>%
    magrittr::set_colnames(c("Peptide","FragLen"))
  cl <- parallel::makeCluster(coreN,type = "PSOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  dt_peptdesc <- foreach::foreach(i=1:nrow(parameterDT), .inorder=F, .options.snow=opts)%dopar%{
    peptide.Frag.single(parameterDT$"Peptide"[[i]], parameterDT$"FragLen"[[i]])
  } %>% unlist()  #Attention, this could not be unique
  close(pb)
  parallel::stopCluster(cl)
  gc();gc()
  # Finish the timer
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")

  # Clear logs
  rm(list=setdiff(ls(), c("dt_peptdesc")))
  return(dt_peptdesc)
}



#### filtering all significant fragment ####
#' Title PepFraSign
#' @export
#' @return d datatable with significant enrichment fragment
#' @examples
#' pepJudgeMatrix <- deredundancy.allmouse %>% column_to_rownames(var = "description")
#' pepFragSignResult <- PepFragSign()
PepFragSign <- function(
  PepFrag = PepFragResult,
  peptides = peptides.allmouse,
  JudgeMatrix = pepJudgeMatrix,
  prop = 0.01,
  coreN =parallel::detectCores(logical=F),
  Fake = FALSE,
  FakeMatrix = NULL
){
  library(tidyverse)
  PropPepFragTable <- table(PepFrag)[table(PepFrag)/length(peptides) > prop]
  PropPepFrag <- names(PropPepFragTable)
  PropPepFrag.single.analysis <- function(Frag,FakeMatrix.single
  ){
    require(tidyverse)
    FragMatch <- grep(pattern = Frag,x = peptides,value = T)
    if(Fake == TRUE){
      colnames(FakeMatrix.single) <- "Fake"
      FragMatch <- FakeMatrix.single %>%
        filter(Fake %in% FragMatch) %>%
        rownames()
      FragMatchAll <- JudgeMatrix[FragMatch,]
    }else{
      FragMatchAll <- JudgeMatrix[FragMatch,]
    }
    N <- length(peptides)
    M <- length(which(JudgeMatrix[,2] == 1))
    n <- length(FragMatch)
    k <- length(which(FragMatchAll[,2] == 1))
    PValue <- phyper(k-1,M,N-M,n,lower.tail = FALSE)
    FragMatchAll.num <- table(FragMatchAll$judge.final)
    dt <- data.table::data.table(Frag = Frag,
                                 FragNum = PropPepFragTable[Frag],
                     FragMatch = FragMatch,
                     FragMatchJudge = FragMatchAll[,2],
                     PValue = PValue)
    return(dt)
  }
  cl <- parallel::makeCluster(coreN,type = "PSOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=length(PropPepFrag), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  dt_peptdesc <- foreach::foreach(i=1:length(PropPepFrag), .inorder=F,.options.snow=opts)%dopar%{
    PropPepFrag.single.analysis(PropPepFrag[i],FakeMatrix)
  } %>% data.table::rbindlist() #Attention, this could not be unique
  parallel::stopCluster(cl)
  gc()
  return(dt_peptdesc)
}


#### peptide description ####
#' Title peptide Feature Description
#' @export
#' @example
PepFeatureDesc <- function(
  peptideSet,
  fragLenSet = 3:8,
  featureSet = NULL,
  coreN = parallel::detectCores(logical=F),
  CutStart = NULL, # a numeric of cut peptides from default is NULL
  Cutend = NULL
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
      Peptides::instaIndex(peptide),
      Peptides::kideraFactors(peptide),
      Peptides::mswhimScores(peptide),
      Peptides::pI(peptide),
      Peptides::protFP(peptide),
      Peptides::vhseScales(peptide),
      Peptides::zScales(peptide)
    ))
  }
  peptideDescriptor.NameSet <- c("AliphaticIndex","BLOSUM1","BLOSUM2","BLOSUM3","BLOSUM4","BLOSUM5","BLOSUM6","BLOSUM7","BLOSUM8","BLOSUM9","BLOSUM10","Boman","Charge","PP1","PP2","PP3","F1","F2","F3","F4","F5","F6","HydrophobicMoment","Hydrophobicity","Instability","KF1","KF2","KF3","KF4","KF5","KF6","KF7","KF8","KF9","KF10","MSWHIM1","MSWHIM2","MSWHIM3","pI","ProtFP1","ProtFP2","ProtFP3","ProtFP4","ProtFP5","ProtFP6","ProtFP7","ProtFP8","VHSE1","VHSE2","VHSE3","VHSE4","VHSE5","VHSE6","VHSE7","VHSE8","Z1","Z2","Z3","Z4","Z5")
  peptideDescriptor.FragStat.Single <- function(peptide, fragLen){
    f <- sapply(1:(max(nchar(peptide))-fragLen+1), function(i){stringr::str_sub(peptide, i, i+fragLen-1)})
    d <- sapply(f[nchar(f)==fragLen], function(s){peptideDescriptor.Batch(s)})
    data.table::data.table(
      "Peptide"=peptide,
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
    parameterDT <- data.table::CJ(peptideSet, fragLenSet) %>%
      magrittr::set_colnames(c("Peptide", "FragLen"))
    cl <- parallel::makeCluster(coreN, type="PSOCK")
    doSNOW::registerDoSNOW(cl)
    sink(tempfile())
    pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
    sink()
    opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
    dt_peptdesc <- foreach::foreach(i=1:nrow(parameterDT), .inorder=F, .options.snow=opts)%dopar%{
      peptideDescriptor.FragStat.Single(parameterDT$"Peptide"[[i]], parameterDT$"FragLen"[[i]])
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

  # Basic peptide characteristics
  dt_basic <- data.table::data.table("Peptide"=peptideSet, "Peptide_Length"=nchar(peptideSet))
  for(aa in Biostrings::AA_STANDARD){
    dt_basic[,paste0("Peptide_Contain", aa):=as.numeric(stringr::str_detect(peptideSet, aa))]
  }
  dt_peptdesc <- merge(dt_basic, dt_peptdesc, by="Peptide")
  data.table::setorder(dt_peptdesc, Peptide)

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




