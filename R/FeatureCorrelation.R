# BASIC INFO ---------------------------
##
## Script name: FeatureCorrelation.R
##
## Purpose of script: This is a series of function for Feature correlation calculation.
##
## Author: JijunYu
##
## Date Created: 2021-12-08
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main -----------------------

library(corrplot)
#' Title
#'
#' @param FeaMatrix a dataframe of feature matrix
#' @param outputname a character of output name
#' @param width a numeric of width
#' @param height a numeric of height
#' @param corvismethod a character of visulization method
#' @param order order to visualization,detail see corrplot
#'
#' @return a figure with png format
#' @export
#'
#' @examples
CorPlot <- function(FeaMatrix,outputname,
                    width = 3000,
                    height = 3000,
                    corvismethod = "circle",
                    order = "original"){
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(FeaMatrix)
  M <- cor(FeaMatrix)
  pdf(file = outputname,
      width = width,
      height = height)
  corrplot(M, method=corvismethod,
           order = order,p.mat = p.mat,
           sig.level = 0.01,col = rev(COL2('RdBu',200)),
          number.cex = 0.75,tl.cex = 1.2,cl.cex = 1.5,pch.cex = 1.6)
  dev.off()
}


#### correlation between each other ####
#INPUT DATA
#' Title
#'
#' @param FeaMatrix a dataframe of feature matrix
#' @param outputname a character of output name
#' @param width a numeric of width
#' @param height a numeric of height
#'
#' @return a dataframe with all correlation between each other which is bigger than 0.8;
#' and a hist plot
#' @export
#'
#' @examples
CorBetween <- function(FeaMatrix,
                       outputname,
                       width = 3000,
                       height = 3000,
                       cutoff = 0.75){
  metagenomics <- FeaMatrix
  cor <- 0
  d<-c("0","0","0")
  end=dim(metagenomics)[2]-1
  for (i in 1:end) {
    x = metagenomics[,i]
    s=i+1
    for (j in s:dim(metagenomics)[2]) {
      cor_n <- cor(x, y = metagenomics[,j], method = "pearson")
      cor<-append(cor,cor_n)
      if(as.numeric(as.character(cor_n))>-2){
        metab1<-colnames(metagenomics)[i]
        metab2<-colnames(metagenomics)[j]
        re<-c(metab1,metab2,cor_n)
        d<-rbind(d,re)
      }
    }
  }

  d<-as.data.frame(d)
  colnames(d)<-c("metab1","metab2","pearson_cor")
  d<-d[-1,]
  cor<-cor[-1]
  pdf(file = outputname,
      width = width,
      height = height)
  hist(cor,breaks=100,xlim=c(-1,1))
  abline(v = cutoff, col = "red")
  dev.off()
  return(d)
}


# https://zhuanlan.zhihu.com/p/94070722
#Point-biserial
#' Title
#'
#' @param FeaMatrix a dataframe of Featrue matrix
#' @param judge a vector of judge
#'
#' @return
#' @export
#'
#' @examples
Pointbiserial <- function(FeaMatrix,judge){
  judge <- as.numeric(judge)
  cor.all <- c()
  pvalue.all <- c()
  for(i in 1:ncol(FeaMatrix)){
    x = FeaMatrix[[i]]
    y = judge
    test_tmp <- cor.test(x,y)
    cor.all <- c(cor.all,test_tmp$estimate)
    pvalue.all <- c(pvalue.all,test_tmp$p.value)
  }
  cor.result <- data.frame(
    feature = colnames(FeaMatrix),
      cor = cor.all,
      pvalue = pvalue.all)
  cor.result <- cor.result %>%
    arrange(cor)
  return(cor.result)
}

#### permutation test #####
#https://blog.csdn.net/zhouyijun888/article/details/69524200
#http://www.360doc.com/content/19/0704/18/40208890_846710655.shtml
#https://blog.csdn.net/wukong1981/article/details/72820049
# Permutation <- function(FeaMatrix,judge){
#   library(coin)
#   judge <- as.numeric(judge)
#   cor.all <- c()
#   pvalue.all <- c()
#   for(i in 1:ncol(FeaMatrix)){
#     colnames(FeaMatrix[,i]
#     y = judge
#     test_tmp <- oneway_test(x~y,data = )
#     cor.all <- c(cor.all,test_tmp$estimate)
#     pvalue.all <- c(pvalue.all,test_tmp$p.value)
#   }
# }

# library(Deducer)
# x<-c(20,34,67,53,12,13,55,89)
#
# y<-c(23,45,12,56,23,67,22,66)
#
# perm.t.test(x,y,alternative = "greater",midp = TRUE,B = 1000) #进行1000次置换检验

