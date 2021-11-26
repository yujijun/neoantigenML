# Description: Neoantigen classifier for machine learning methods.
# Update Time: Thu Nov 11 00:17:57 2021
# This is script for mlr3 package learning

#' Title
#'
#' @param Neodataset A data frame of feature selection dataset
#'
#' @return
#' @export
#'
#' @examples
#' neoML(Neodataset)
neoML <- function(Neodataset){
  #### package ####
  library(tidyverse)
  library(mlr3)
  library(mlr3verse)
  ## task： 封装了数据及额外信息 如预测目标 ####
  task = TaskClassif$new(id="neodata", backend=Neodataset, target="judge")
  ## 选择学习器：注意超参数 ####
  learner = lrn("classif.rpart" , cp=0.1, minsplit=10)
  ## 划分训练集和测试集 ####
  set.seed(123)
  train = sample(task$nrow, 0.8 * task$nrow)  # 行的指标
  test = setdiff(seq_len(task$nrow), train)
  ## 模型训练 ####
  learner$train(task, row_ids=train)
  ## 模型预测 ####
  predictions = learner$predict(task, row_ids=test)
  ## 模型评估 ####
  predictions$confusion
  measure = msr("classif.acc")
  predictions$score(measure)     # 预测精度
  ## 重采样 ####
  resampling = rsmp("cv", folds = 3L)
  rr = resample(task, learner, resampling)
  rr$score(measure)
  rr$aggregate(measure)
}



#### reference:####
# mlr3verse learn resource:https://www.zhihu.com/question/403123109
# mlr3 case:https:/
#/zhuanlan.zhihu.com/p/112845336
