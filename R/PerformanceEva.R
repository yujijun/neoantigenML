# BASIC INFO ---------------------------
##
## Script name:BenchmarkEva.R
##
## Purpose of script: Evaluation different methods by benchmark methods.
##
## Author: JijunYu
##
## Date Created: 2021-11-28
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main:Performance Evaluation --------------------
#' Title Performance evaluation of single machine learning model.
#'
#' @param task A task object created by mlr3
#' @param leaner A learner object created by mlr3
#' @param resamplemethod A resample object created by mlr3
#' @param bmr_small A single benchmark object such as  bmr_small = bmr$clone()$filter(task_id = tasks[1])
#'
#' @return A list with all performance results.
#' @export
#'
#' @examples
#' data("Sonar", package = "mlbench")
#' task = as_task_classif(Sonar, target = "Class", positive = "M")
#' learner = lrn("classif.rpart", predict_type = "prob")
#' PerformanceEva(task = task,learner = learner)
PerformanceEva <- function(task,
                           learner,
                           resamplemethod = rsmp("cv",folds =3),
                           bmr_small){
  require(mlr3verse)
  future::plan("multisession")
  if(is.empty(bmr_small) == TRUE){
    rr = resample(task, learner, resamplemethod, store_models = TRUE)
  }else{
    rr = bmr_small
  }
  autoplot(rr, measure = msr("classif.auc"))
  roc = autoplot(rr,type = "roc")
  roc.best = autoplot(roc$prediction(), type = "roc")
  prc = autoplot(rr,type = "prc")
  prediction = autoplot(rr, type = "prediction")
  Allperformance = list(roc,roc.best,prc,prediction)
  return(Allperformance)
}

#--------------Main:benchmark calculation-----------------------
#' Title Benchmark calculation
#' In particular unique combination of Task,Learner and Resampling.

#' @param tasks A vector of task
#'
#' @param resamplings A resamplings methods
#' @param learners A list of learners
#'
#' @return
#' @export
#'
#' @examples
#' design = benchmark_grid(
#' tasks = tsks(c("spam", "german_credit", "sonar")),
#' learners = lrns(c("classif.glmnet", "classif.rpart", "classif.featureless"),
#'                predict_type = "prob", predict_sets = c("train", "test")),
#' resamplings = rsmps("cv", folds = 3)
#' classif.glmnet.rbv2,classif.kknn.rbv2,classif.ranger.rbv2,classif.svm.rbv2,classif.xgboost.rbv2,
BenchmarkEva <- function(tasks,
                         learners,
                         resamplings){
  require("mlr3verse")
  require("data.table")
  future::plan("multisession")
  set.seed(123)
  design = benchmark_grid(
    tasks = tasks,
    learners = learners,
    resamplings = resamplings
  )
  bmr = benchmark(design)
  measures = list(
    msr("classif.auc", predict_sets = "train", id = "auc_train"),
    msr("classif.auc", predict_sets = "test",id = "auc_test")
  )
  tab = bmr$aggregate(measures)
  #ranks = tab[, .(learner_id, rank_train = rank(-auc_train), rank_test = rank(-auc_test)), by = task_id]
  #ranks = ranks[, .(mrank_train = mean(rank_train), mrank_test = mean(rank_test)), by = learner_id]
  bmr.barplot = autoplot(bmr) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # bmr_small = bmr$clone()$filter(task_id = tasks[1])
  # bmr_roc = autoplot(bmr_small, type = "roc")
  # Extracting ResampleResults:https://mlr3book.mlr-org.com/perf-eval-cmp.html
  Allresult <- list(bmr = bmr,tab = tab,barplot = bmr.barplot)
  return(Allresult)
}

