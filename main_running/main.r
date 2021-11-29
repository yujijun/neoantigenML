# BASIC INFO ---------------------------
##
## Script name: main.r
##
## Purpose of script: This is main function for machine learning running
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
## Xgboost xgboost need more resource

# SET WORKING DIRECTORY -----------------------------
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "/cloud/project"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
output_path <- "./result/yjj/"
## SET OPTIONS ---------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse",
              "devtools",
              "mlr3verse") # list of packages to load
n_packages <- length(packages) # count how many packages are required
new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}
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
output.property <- PropertyofPepSingle(peptides = peptides,PropAll = ManualPropertySelection)
#### Create machine learning dataset ----need new dataset----####
MLtestData <- output.property %>%
  mutate(judge = as.factor(Neodataset$judge)) ## Why some is NA
MLtestData <- MLtestData[which(rowSums(is.na(MLtestData)) == 0),]
#### Feature Selection ---- 完善即可 ---- ####
FeatureSeByFilter.result <- FeatureSeByFilter(dataset = MLtestData)
FeatureSebyVariImpFil.result <- FeatureSebyVariImpFil(dataset = MLtestData)
#### Feature Selection by wrapper---- 完善即可----####

# To do Model Feature Selection based on trained model before.
rr = FeatureSebyNestedresam(MLtestData)

#### initial Training by pipeline: ---完善即可---####
# purpose: determination of data proprocess 判断是否需要进行参数预处理
# test parameter #
Neodataset = MLtestData
taskid = "neoML"
target="judge"
branchName = "nop"
Training.rpart <- neoML.rpart(MLtestData)

#### Model Tuning ----完善即可 ----####
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
Model.xgboost <- neoML.all(MLtestData,taskid = "neoML.xgboost",
                        tunningkey = "classif.xgboost.default",
                        measures = "classif.auc")
#### Performance Evaluation and Comparison by benchmark ####
Neodataset = MLtestData
taskid = "neoML.all"
target="judge"
tasks = TaskClassif$new(id=taskid, backend=Neodataset, target = target)
learners = list(Model.glmnet$learner,
              Model.kknn$learner,
              Model.ranger$learner,
              Model.svm$learner,
              Model.xgboost$learner)
resamplings = rsmps("cv", folds = 3)
bench <- BenchmarkEva(tasks = tasks,
             learners = learners,
             resamplings = resamplings)
bmr_small = bench[[1]]$clone()$filter(task_id = "neoML.all")
print(bmr_small)
#### Model Interpretation ----代码编写【简单】独立于机器学习，不着急所有运行完毕之后再进行这一步----####
#ref:https://mlr3book.mlr-org.com/interpretation.html
#### output ####
# purpose: Export and save all results.
write.csv(output.property,file = paste0(output_path,
                                          "output.property.csv"),
          quote = F)


