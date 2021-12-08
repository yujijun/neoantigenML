# BASIC INFO ---------------------------
##
## Script name: Neoantigen classifier for machine learning methods.
##
## Purpose of script:This is script for mlr3 package learning
##
## Author: JijunYu
##
## Date Created: 2021-11-29
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
## https://github.com/mlr-org/mlr3tuningspaces
## https://mlr3book.mlr-org.com/appendix.html?q=measure#list-measures
#--------------Main -----------------------

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
  future::plan("multisession")
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

#' Title glmnet/kknn/range/svm/xgboost training model
#'
#' @param Neodataset A data frame of feature selection dataset.
#' @param taskid A character of task id.
#' @param target A column name of prediction target in Neodataset
#' @param tunningkey A character of mlr_tuning_spaces c(classif.glmnet.rbv2,classif.kknn.rbv2,classif.ranger.rbv2,classif.svm.rbv2,classif.xgboost.default)
#' @param evalsnum  the numbers of iterations evaluation
#' @param innercv The number of inner sampling in Nested sampling methods, default is 4
#' @param outercv The number of outer sampling in Nested sampling methods, default is 4
#' @param innerparall The number of inner parallelizations,default is 4
#' @param measures A character for Performance Measures detail could be found in https://mlr3book.mlr-org.com/appendix.html?q=measure#list-measures
#' @param outerparall The number of outer parallelizations,default is 4
#'
#' @return A resampling model with model.store
#' @export
#'
#' @examples
#' neoML.glmnet(MLtestData)
#' The evaluation of Tunning Space
neoML.all <- function(Neodataset = MLtestData,
                         taskid = "neoML.glmnet",
                         target="judge",
                         tunningkey = "classif.glmnet.rbv2",
                         evalsnum = 20,
                         innerparall = 4,
                         outerparall = 3,
                         innercv = 4,
                         outercv = 3,
                         measures = "classif.auc") {
  # package ####
  library(tidyverse)
  library(mlr3verse)
  library(magrittr)
  library(mlr3tuningspaces)
  library("glmnet")
  set.seed(123)
  # create task ####
  task = TaskClassif$new(id=taskid, backend=Neodataset, target = target)
  # write auto tuning learning function by parallelization methods
  tuning_space = lts(tunningkey)
  learner = tuning_space$get_learner()
  if(tunningkey == "classif.svm.rbv2"){
    learner$param_set$values$type = "C-classification"
  }
  learner$predict_type = "prob"
  learner$predict_sets = c("train", "test")
  # write Nested resampling function by Explicit parallelization methods.
  # Runs both loops in parallel
  # Ref:https://mlr3book.mlr-org.com/optimization.html#nested-resampling
  # https://github.com/mlr-org/paradox/issues/259
  future::plan(list(future::tweak("multisession", workers = innerparall),
                    future::tweak("multisession", workers = outerparall)))
  resampling = rsmp("cv", folds = innercv)
  measure = msr(measures)
  terminator = trm("evals", n_evals = evalsnum)
  tuner = tnr("grid_search", resolution = 10)
  at = AutoTuner$new(learner, resampling, measure, terminator, tuner)
  outer_resampling = rsmp("cv", folds = outercv)
  rr = resample(task, at, outer_resampling, store_models = TRUE)
  return(rr)
}

