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
output_path <- "./result/yjj/"
.libPaths("/public/home/yujijun01/R/x86_64-conda-linux-gnu-library/4.1")
system("conda info")
## SET OPTIONS ---------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

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

#### Original dataset input ####
#### Peptide visualization ----代码编写【简单】----####
#### Property calculation and selected manually ---完善即可----- ####
ManualPropertySelection <- c("aaComp","boman","crucianiProperties",
                             "fasgaiVectors","instaIndex","kideraFactors",
                             "protFP","stScales","tScales")
cat("Calculating property of peptides...\n\n", sep = "")
output.property <- PropertyofPepSingle(peptides = peptides,PropAll = ManualPropertySelection)
MLtestData <- output.property %>%
  mutate(judge = as.factor(Neodataset$judge)) ## Why some is NA
MLtestData <- MLtestData[which(rowSums(is.na(MLtestData)) == 0),]
MLtestData <- janitor::clean_names(MLtestData)
#### Feature Selection --------- ####
# cat("Feature selection by Filter method...\n\n", sep = "")
# FeatureSeByFilter.result <- FeatureSeByFilter(dataset = MLtestData)
# cat("Feature selection by Variant important Fileter method...\n\n", sep = "")
# FeatureSebyVariImpFil.result <- FeatureSebyVariImpFil(dataset = MLtestData)
# # To do Model Feature Selection based on trained model before.
# cat("Feature selection by byNestedresam...\n\n", sep = "")
# FeatureSebyNestedresam.result = FeatureSebyNestedresam(MLtestData,nevals = 20)
# cat("Save Feature selection Result...\n\n", sep = "")
# write.csv(output.property,file = paste0(output_path,
#                                         "output.property.csv"),
#           quote = F)
#
# save(output.property, file =  paste0(output_path,"output.property.rda"))
# save(MLtestData,file =  paste0(output_path,"MLtestData.rda"))
# save(FeatureSeByFilter.result,file = paste0(output_path,"FeatureSeByFilter.result.rda"))
# save(FeatureSebyVariImpFil.result, file = paste0(output_path,"FeatureSebyVariImpFil.result.rda"))
# save(FeatureSebyNestedresam.result,file = paste0(output_path,"FeatureSebyVariImpFil.result.rda"))

#### initial Training by pipeline: ---完善即可---####
# purpose: determination of data proprocess 判断是否需要进行参数预处理
# test parameter #
# Neodataset = MLtestData
# taskid = "neoML"
# target="judge"
# branchName = "nop"
# Training.rpart <- neoML.rpart(MLtestData)

#### Evaluation of Model Tuning --- ----####
cat("Machine Learning is running...\n\n", sep = "")
Model.glmnet <- neoML.all(MLtestData)
Model.kknn <- neoML.all(MLtestData,taskid = "neoML.kknn",
                        tunningkey = "classif.kknn.rbv2",
                        measures = "classif.auc")
Model.ranger <- neoML.all(MLtestData,taskid = "neoML.ranger",
                        tunningkey = "classif.ranger.rbv2",
                        measures = "classif.auc")
Model.svm <- neoML.all(MLtestData,taskid = "neoML.svm",
                        tunningkey = "classif.svm.rbv2",
                        measures = "classif.auc")
# Model.xgboost <- neoML.all(MLtestData,taskid = "neoML.xgboost",
#                         tunningkey = "classif.xgboost.default",
#                         measures = "classif.auc")
Model.list <- list(
  Model.glmnet = Model.glmnet,
  Model.kknn = Model.kknn,
  Model.ranger = Model.ranger,
  Model.svm = Model.svm
  #Model.xgboost = Model.xgboost
)
save(Model.list,file = paste0(output_path,"Model.list.rda"))
cat("Greate, All models have been created, benchmark evaluation is running...\n\n", sep = "")
#### Training of Model to get final Model(with certain hyperparameter)-----####
#### Performance Evaluation and Comparison by benchmark (with certain hyperparameter and parameter)####
#### load out final model ####
#load(paste0(output_path,"Model.list.rda"))
Neodataset = MLtestData
taskid = "neoML.all"
target="judge"
tasks = TaskClassif$new(id=taskid, backend=Neodataset, target = target)
learners = list(Model.glmnet$learner,
              Model.kknn$learner,
              Model.ranger$learner,
              Model.svm$learner)
resamplings = rsmps("cv", folds = 3)
bench <- BenchmarkEva(tasks = tasks,
             learners = learners,
             resamplings = resamplings)
# bmr_small = bench$bmr$clone()$filter(task_id = "neoML.all")
# autoplot(bmr_small, type = "roc")
#### Model Interpretation ----代码编写【简单】独立于机器学习，不着急所有运行完毕之后再进行这一步----####
#ref:https://mlr3book.mlr-org.com/interpretation.html
#### output ####
# purpose: Export and save all results.

save(bench,file = paste0(output_path,"bench.rda"))
cat("Congratulation,All result have been generated...\n\n", sep = "")





