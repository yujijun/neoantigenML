# BASIC INFO ---------------------------
##
## Script name:main.r
##
## Purpose of script: This is main script for NeoantigenML running
##
## Author: JijunYu
##
## Date Created: 2021-11-27
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##

# SET WORKING DIRECTORY -----------------------------
print(.libPaths())
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd = "/public/home/yujijun01/NeoantigenML"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
output_path <- "./result/yjj_20211211/"
.libPaths("/public/home/yujijun01/R/x86_64-conda-linux-gnu-library/4.1")
system("conda info")
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
# new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# # install missing packages
# if(length(new.pkg)){
#   install.packages(new.pkg,repos = "https://mirrors.bfsu.edu.cn/CRAN/")
# }
# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}
#--------------Main -----------------------
load_all()

#### load dataset of figure matrix ####
#if datasets have been included into env,there is no need to be loaded.
#### Pretraining of Model Tuning[evaluation by train and test dataset] --- ----####
cat("Machine Learning is running...\n\n", sep = "")
Model.glmnet <- neoML.all(mousetestML)
Model.kknn <- neoML.all(mousetestML,taskid = "neoML.kknn",
                        tunningkey = "classif.kknn.rbv2",
                        measures = "classif.auc")
Model.ranger <- neoML.all(mousetestML,taskid = "neoML.ranger",
                        tunningkey = "classif.ranger.rbv2",
                        measures = "classif.auc")
Model.svm <- neoML.all(mousetestML,taskid = "neoML.svm",
                        tunningkey = "classif.svm.rbv2",
                        measures = "classif.auc")
# Model.xgboost <- neoML.all(mousetestML,taskid = "neoML.xgboost",
#                         tunningkey = "classif.xgboost.default",
#                         measures = "classif.auc")
Model.list <- list(
  Model.glmnet = Model.glmnet,
  Model.kknn = Model.kknn,
  Model.ranger = Model.ranger,
  Model.svm = Model.svm
  #Model.xgboost = Model.xgboost
)
save(Model.list,file = paste0(output_path,"Model.list_alldataset.rda"))
cat("Greate, All models have been pre tested and evaluated...\n\n", sep = "")

# #### Training of Model by training dataset to get final Model(with certain hyperparameter)-----####
# # Training model by whole datasets
# WholeTrainDataset <- mousetestML
# task = TaskClassif$new(id="mousetestML", backend=WholeTrainDataset, target = "judge")
# # train_set = sample(task$nrow,0.8*task$nrow)
# # test_set = setdiff(seq_len(task$nrow),train_set)
# glmnet.learner <- Model.glmnet$at$train(task)
# kknn.learner <- Model.kknn$at$train(task,row_ids = train_set)
# ranger.learner <- Model.ranger$at$train(task,row_ids = train_set)
# svm.learner <- Model.svm$at$train(task,row_ids = train_set)
#
# #### Performance Evaluation and Comparison by benchmark dataset####
# #verified model by new verified datasets
# Neodataset = mousetestML[test_set,]
# taskid = "neoML.all"
# target="judge"
# tasks.test = TaskClassif$new(id=taskid, backend=Neodataset, target = target)
# learners = list(glmnet.learner,
#                 kknn.learner,
#                 ranger.learner,
#                 svm.learner)
# resamplings = rsmps("cv", folds = 10)
# #resamplings = "none"
# bench <- BenchmarkEva(tasks = tasks.test,
#              learners = learners,
#              resamplings = resamplings)
# # bmr_small = bench$bmr$clone()$filter(task_id = "neoML.all")
# # autoplot(bmr_small, type = "roc")
# #### Model Interpretation ----代码编写【简单】独立于机器学习，不着急所有运行完毕之后再进行这一步----####
# #ref:https://mlr3book.mlr-org.com/interpretation.html
# #### output ####
# save(bench,file = paste0(output_path,"bench.rda"))
# cat("Congratulation,All result have been generated...\n\n", sep = "")

