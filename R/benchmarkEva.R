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
#--------------Main -----------------------
#' Title
#'
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param learnerlist A vector of learners which would like to be compared
#'
#' @return
#' @export
#'
#' @examples
BenchmarkEva <- function(dataset = MLtestData,
                         target = "judge",
                         taskid = "MLtest",
                         learnerlist = list(lrn("classif.rpart"),lrn("classif.featureless"))){
  require("mlr3verse")
  grid = benchmark_grid(
    task = TaskClassif$new(id=taskid, backend=dataset, target = target),
    learner = learnerlist,
    resampling = rsmp("cv",folds = 3)
  )
  bmr = benchmark(grid, store_models = TRUE)
  return(bmr)
}

