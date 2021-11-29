# BASIC INFO ---------------------------
##
## Script name: HyperparameterTuning.R
##
## Purpose of script: Model tuning
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
## Different machine learning method have different tunning parameter
#--------------Main -----------------------
# Should write auto model
Tunning <- function(dataset = MLtestData,
                    target = "judge",
                    taskid = "MLtest",
                    learnmethod = "classif.rpart"){
  require("mlr3verse")
  require("FSelectorRcpp")
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  learner = lrn(learnmethod)
  search_space = ps(
    cp = p_dbl(lower = 0.001, upper = 0.1),
    minsplit = p_int(lower = 1, upper = 10))
  hout = rsmp("holdout")
  measure = msr("classif.ce")
  evals20 = trm("evals", n_evals = 20)
  instance = TuningInstanceSingleCrit$new(
    task = task,
    learner = learner,
    resampling = hout,
    measure = measure,
    search_space = search_space, # different model need different search_space
    terminator = evals20
  )
  tuner = tnr("grid_search", resolution = 5)
  tuner$optimize(instance)
  return(tuner)
}
