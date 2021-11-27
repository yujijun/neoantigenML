# Description: Neoantigen classifier for machine learning methods.
# Update Time: Fri Nov 26 18:20:12 2021
# This is script for mlr3 package learning


#' Title Classify by classif.rpart
#' @param Neodataset A data frame of feature selection dataset.
#' @param taskid A character of task id.
#' @param target A column name of prediction target in Neodataset
#' @param branchName A character  Branch to use, which should be one of "nop", "pca", "scale"
#'
#' @return A matrix with accurate info
#' @export
#'
#' @examples
#' neoML.rpart(MLtestData)
neoML.rpart <- function(Neodataset,
                  taskid = "neoML",
                  target="judge",
                  branchName = "nop"){
  # package ####
  library(tidyverse)
  library(mlr3verse)
  library(magrittr)
  # create task ####
  task = TaskClassif$new(id=taskid, backend=Neodataset, target = target)
  # Building PipeOps####
  ## preprocess
  opts <- list(po("nop", "nop"), po("pca"), po("scale")) # scale or not
  opt_ids <- mlr3misc::map_chr(opts, `[[`, "id")
  preprocess <- po("branch", options = opt_ids) %>>%
    gunion(opts) %>>%
    po("unbranch", options = opt_ids)
  ## Model ensembles in sample model
  single_pred <- po("subsample", frac = 0.7) %>>%
    po("learner", lrn("classif.rpart"))
  pred_set <- ppl("greplicate", single_pred, 10L)
  bagging <- pred_set %>>%
    po("classifavg", innum = 10)
  ## Chaining POs
  graph <- preprocess %>>%  bagging
  graph$param_set$values$branch.selection = branchName
  ## graph training
  learner = as_learner(graph)
  set.seed(123)
  train = sample(task$nrow, 0.8 * task$nrow)
  test = setdiff(seq_len(task$nrow), train)
  learner$train(task, row_ids=train) # training
  predictions = learner$predict(task, row_ids=test)  ## prediction
  ## resampleing and evaluation  ####
  resampling = rsmp("cv", folds = 10L)
  rr = resample(task, learner, resampling)
  acc.all <- rr$score(measure)
  acc <- rr$aggregate(measure)
  acc.all$mean = acc
  return(acc.all)
}





## multimodel
# multimodel <-
#   gunion(list(
#     po("learner", lrn("classif.rpart")),
#     po("learner",lrn("classif.xgboost")),
#     po("learner",lrn("classif.kknn")),
#     po("learner",lrn("classif.svm")),
#     po("learner",lrn("classif.glmnet")),
#     po("learner",lrn("classif.naive_bayes"))
#   ))
#### reference:####
# mlr3verse learn resource:https://www.zhihu.com/question/403123109
# mlr3 case:https:/
#/zhuanlan.zhihu.com/p/112845336
