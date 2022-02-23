# BASIC INFO ---------------------------
##
## Script name:PepFeatureSelect.R
##
## Purpose of script: Feature Selection of peptides' properties.
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
## Ref:https://mlr3book.mlr-org.com/optimization.html#fs
#--------------Main -----------------------
# Feature Selection by filter methods. ####
#' Title
#'
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param parallele A logical TRUE or FALSE， Whether to parallelize. default is FALSE in mlr3 package
#' @param nthreads if parallele is TRUE, the nthreads will be worked.
#' @param filtermethod A vector of filtering methods, each character should be one of mlr_filters$keys,which filter methods should be used,all filters method could be found by mlr_filters_,
#'
#' @return
#' @export
#'
#' @examples
FeatureSeByFilter <- function(dataset = MLtestData,
                              target = "judge",
                              taskid = "FeatureSeByFilter",
                              filtermethod = c("anova","auc",
                                               "cmim","disr","find_correlation",
                                               "importance","information_gain",
                                               "jmi","jmim","kruskal_test",
                                               "mim","mrmr","njmim",
                                               "performance","permutation","relief",
                                               "variance"),
                              parallele = FALSE,
                              nthreads = 4){
  require("mlr3verse")
  require("FSelectorRcpp")
  set.seed(123)
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  Fillist <- list()
  for (i in filtermethod) {
    print(paste0("Filtering by ",i))
    Fil <- flt(i)
    if(parallele == TRUE){
      set_threads(Fil,n=nthreads)
    }
    Fil$calculate(task)
    Fil.result <- as.data.table(Fil)
    Fillist[[i]] <- Fil.result
  }
    Filall <- do.call(cbind.data.frame,Fillist)
    return(Filall)
}
# Feature Selection by Variable Importance Filters ####
## Learners With Embedded Filter Methods
##  [1] "classif.featureless" "classif.ranger"      "classif.rpart"
##  [4] "classif.xgboost"     "regr.featureless"    "regr.ranger"
##  [7] "regr.rpart"          "regr.xgboost"        "surv.ranger"
## [10] "surv.rpart"          "surv.xgboost"
## https://mlr3book.mlr-org.com/appendix.html#fs-filter-embedded-list

#' Title Feature Selection by Embedded Filter Methods
#'
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param filtermethodVI A vectors of Learners With Embedded Filter Methods
#'
#' @return
#' @export
#'
#' @examples
FeatureSebyVariImpFil <- function(dataset = MLtestData,
                                  target = "judge",
                                  taskid = "FeatureSebyVariImpFil",
                                  filtermethodVI = c("classif.ranger",
                                                     "classif.rpart",
                                                     "classif.xgboost")){
  require("mlr3verse")
  require("FSelectorRcpp")
  require("ranger")
  future::plan("multisession")
  set.seed(123)
  importantlist <- list()
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  for(i in filtermethodVI){
    print(i)
    if(i == "classif.ranger"){
      lrn = lrn(i, importance = "impurity")
    }else{
      lrn = lrn(i)
    }
    Fil = flt("importance", learner = lrn)
    Fil$calculate(task)
    Filresult <- as.data.table(Fil)
    importantlist[[i]] <- Filresult
  }
  importanttable <- do.call(cbind.data.frame,importantlist)
  return(importanttable)
}
# Feature Selection by Wrapper Methods ####

#' Title Feature Selection by Wrapper Methods
#' Use this function in specific ML methods
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param filtermethodlrn A filter selection by which lrn
#' @param nevals A number of nevals
#'
#' @return
#' @export
#'
#' @examples
FeatureWrapperMethod <- function(dataset = MLtestData,
                                 target = "judge",
                                 taskid = "FeatureWrapperMethod",
                                 filtermethodlrn = "FeatureWrapperMethod",
                                 nevals = 20){
  require("mlr3verse")
  require("data.table")
  set.seed(123)
  future::plan("multisession")
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  terminator = trm("evals", n_evals = nevals)
  learner = lrn("classif.rpart")
  hout = rsmp("holdout")
  measure = msr("classif.ce")
  instance = FSelectInstanceSingleCrit$new(
    task = task,
    learner = learner,
    resampling = hout,
    measure = measure,
    terminator = terminator
  )
  fselector = fs("random_search")
  # reduce logging output
  lgr::get_logger("bbotk")$set_threshold("warn")
  fselector$optimize(instance)
  instanceresult <- as.data.table(instance$archive)
  return(instanceresult)
  #Now the optimized feature subset can be used to subset the task and fit the model on all observations.


}

# Feature Selection by Automating Method ####
#' Title Automating the Feature Selection and benchmark judgement
#'
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param filtermethodlrn A filter selection by which lrn
#' @param nevals A number of nevals
#' @param fselector A character of fselector
#'
#' @return
#' @export
#'
#' @examples
FeatureSebyAuto <- function(dataset = MLtestData,
                            target = "judge",
                            taskid = "FeatureSebyAuto",
                            filtermethodlrn = "classif.rpart",
                            nevals = 20,
                            fselector = "random_search"){
  require("mlr3verse")
  require("data.table")
  future::plan("multisession")
  set.seed(123)
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  learner = lrn(filtermethodlrn)
  terminator = trm("evals", n_evals = nevals)
  fselector = fs(fselector)
  at = AutoFSelector$new(
    learner = learner,
    resampling = rsmp("holdout"),
    measure = msr("classif.ce"),
    terminator = terminator,
    fselector = fselector
  )
  grid = benchmark_grid(
    task = task,
    learner = list(at, lrn(filtermethodlrn)),
    resampling = rsmp("cv", folds = 3)
  )

  bmr = benchmark(grid, store_models = TRUE)
  bmr$aggregate(msrs(c("classif.ce", "time_train")))
  return(bmr)
}


#' Title
#'
#' @param dataset A dataframe of Machine Learing test dataset
#' @param target A character of column name in dataset
#' @param taskid A character of taskid
#' @param filtermethodlrn A filter selection by which lrn
#' @param nevals A number of nevals
#' @param fselector A character of fselector
#'
#' @return
#' @export
#'
#' @examples
#' rr = FeatureSebyNestedresam(MLtestData)
#' extract_inner_fselect_archives(rr) #https://mlr3fselect.mlr-org.com/
FeatureSebyNestedresam <- function(dataset = MLtestData,
                                  target = "judge",
                                  taskid = "FeatureSebyNestedresam",
                                  filtermethodlrn = "classif.rpart",
                                  nevals = 20,
                                  fselector = "random_search"){
  require("mlr3verse")
  require("data.table")
  require("mlr3fselect")
  future::plan("multisession")
  set.seed(123)
  task = TaskClassif$new(id=taskid, backend=dataset, target = target)
  learner = lrn(filtermethodlrn)
  # nested resampling
  rr = fselect_nested(
    method = fselector,
    task =  task,
    learner = lrn(filtermethodlrn),
    inner_resampling = rsmp("holdout"),
    outer_resampling = rsmp("cv", folds = 3),
    measure = msr("classif.ce"),
    term_evals = nevals,
    batch_size = 5
  )
 return(rr)
}
