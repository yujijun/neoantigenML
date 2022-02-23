# BASIC INFO ---------------------------
##
## Script name: ml_kmer_human_pepdesc.R
##
## Purpose of script: This is main script for NeoantigenML running in Human_PepDesc.RData
##
## Author: JijunYu
##
## Date Created: 2021-11-27
## Update Date:Wed Jan 19 07:28:59 2022
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##

#************************SET WORKING DIRECTORY ***************************####
#print(.libPaths())

require(tidyverse)
require(foreach)
require(doSNOW)
library(parallel)
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd = "/public/home/yujijun01/01_NeoantigenML"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")

input_path <- "./result_v3/01_mutate_wholepeptide_single_diffrelation_result/"
output_base_path <- "./result_v3/01_mutate_wholepeptide_single_diffrelation_result/"
All_files <- list.files(input_path)
All_files <- All_files[6]
coreN = 8
parameterDT <- data.table::CJ(All_files,output_base_path) %>%
magrittr::set_colnames(c("files","output_path"))
cl <- parallel::makeCluster(coreN,type = "PSOCK")
doSNOW::registerDoSNOW(cl)
sink(tempfile())
pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
sink()
opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
ML_fun <- function(i,output_base_path){
  #.libPaths("/public/home/yujijun01/R/x86_64-conda-linux-gnu-library/4.1")
  #system("conda info")

  ## SET OPTIONS ---------------------------
  cat("SETTING OPTIONS... \n\n", sep = "")
  options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
  options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
  Sys.setenv(TZ='Asia/Chongqing')
  # INSTALL PACKAGES & LOAD LIBRARIES -----------------
  cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
  packages <- c("tidyverse",
                "stringr",
                "devtools",
                "mlr3verse",
                "praznik",
                "kknn",
                "xgboost",
                "forestmodel",
                "rpart",
                "stats",
                "FSelectorRcpp",
                "care",
                "iml",
                "DALEX",
                "apcluster",
                "mlr3tuningspaces",
                "precrec",
                "janitor") # list of packages to load
  n_packages <- length(packages) # count how many packages are required
  for(n in 1:n_packages){
    cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
    lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
    eval(parse(text = lib_load)) # evaluate the string to load the library
  }
  #--------------Main -----------------------
  load_all()
  filesname <- basename(i)
  filesname <- stringr::str_remove(filesname,pattern = ".RData")
  #output_folder <- filesname
  output_prefix <- filesname
  load(paste0(input_path,i))
  MLmatrix <- mutateML
  if(ncol(MLmatrix) <=2){
    next
  }
  if(dir.exists(paths = output_base_path) == FALSE){
    dir.create(path = paste0(output_base_path),recursive = TRUE)
  }
  #### Pretraining of Model Tuning[evaluation by train and test dataset] --- ----####
  cat("Machine Learning is running...\n\n", sep = "")
  cat(c("basename is",filesname,"\n\n"))
  Model.glmnet <- neoML.all(MLmatrix)
  Model.kknn <- neoML.all(MLmatrix,taskid = "neoML.kknn",
                          tunningkey = "classif.kknn.rbv2",
                          measures = "classif.auc")
  Model.ranger <- neoML.all(MLmatrix,taskid = "neoML.ranger",
                            tunningkey = "classif.ranger.rbv2",
                            measures = "classif.auc")
  Model.svm <- neoML.all(MLmatrix,taskid = "neoML.svm",
                         tunningkey = "classif.svm.rbv2",
                         measures = "classif.auc")
  # Model.xgboost <- neoML.all(MLmatrix,taskid = "neoML.xgboost",
  #                            tunningkey = "classif.xgboost.default",
  #                            measures = "classif.auc")
  MLmatrix_Modellist <- list(
    Model.glmnet = Model.glmnet,
    Model.kknn = Model.kknn,
    Model.ranger = Model.ranger,
    Model.svm = Model.svm
  )
  save(MLmatrix_Modellist,
       file = paste0(output_base_path,output_prefix,"_Modellist.rda"))

  cat("Greate, All models have been pre tested and evaluated...\n\n", sep = "")
  #load(paste0(output_path,"MLmatrix_Modellist.rda"))
  outer.list <- list()
  inner.list <- list()
  for(i in 1:length(MLmatrix_Modellist)){
    model.name <- names(MLmatrix_Modellist)
    i_inner_tmp <-extract_inner_tuning_results(MLmatrix_Modellist[[i]]$rr)
    i_inner_tmp <- i_inner_tmp %>%
      select(iteration,classif.auc,task_id)
    inner.list[[model.name[i]]] <- i_inner_tmp
    i_outer_tmp <-MLmatrix_Modellist[[i]]$rr$score(msr("classif.auc"))
    i_outer_tmp <- i_outer_tmp %>%
      select(iteration,classif.auc,task_id)
    colnames(i_outer_tmp) <- paste0("outer_",colnames(i_outer_tmp))
    outer.list[[model.name[i]]] <- i_outer_tmp
  }
  inner.all <- do.call(rbind.data.frame,inner.list)
  outer.all <- do.call(rbind.data.frame,outer.list)
  evaluation.all <- cbind(inner.all,outer.all)
  write_tsv(evaluation.all,file = paste0(output_base_path,output_prefix,"_model_evaluation.tsv"))
}
foreach::foreach(i=1:nrow(parameterDT), .inorder=F, .options.snow=opts)%dopar%{
ML_fun(parameterDT$"files"[[i]],parameterDT$"output_path"[[i]])}
close(pb)
parallel::stopCluster(cl)
gc();gc()







