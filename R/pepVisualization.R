## ---------------------------
##
## Script name: pepVisualization.R
##
## Purpose of script: This is a function for peptide visualization
##
## Author: Dr. JijunYu
##
## Date Created: 2021-12-22
##
## Copyright (c) Timothy Farewell, 2021
## Email: jijunyuedu@outlook.com
##
## ---------------------------
##
## Notes:
## Input is from PepFragSign function which is in pepFeatureNew.R
##
## ---------------------------

#' Title visualization for fragment of peptides
#'
#' @param PepFragChoose a vector with all signaficant enriched fragment of peptides
#' @param peptideSet A vector of all peptides with uesed to generate fragment
#' @param pepJudgeMatrix a dataframe with judge info of peptides
#' @param outPath output path
#' @param prefix a character of prefix of output name
#'
#' @return
#' @export
#'
#' @examples
pepFragVis <- function(PepFragChoose,
                       peptideSet = peptides.allhuman,
                       pepJudgeMatrix = PepjudgeMatrixexam,
                       outPath = "./result/Kmer_result/",
                       prefix = "Kmer100_"){
  library("treeio")
  library("Biostrings")
  library("ggtree")
  library("ggmsa")
  library("seqmagick")
  library("cowplot")
  library("ggplot2")
  pepFragVis.singleFrag <- function(Frag){

    FragMatchPep <- grep(pattern = Frag,x = peptideSet,value = T)
    names(FragMatchPep) <- FragMatchPep
    FragMatchPepAAstring <- AAStringSet(x = FragMatchPep)
    FragMatchPepalign <- msa(FragMatchPepAAstring)
    FragMatchPepalign.BstringSet <- as(FragMatchPepalign,"BStringSet")
    ### tree plot ##
    d <- as.dist(stringDist(FragMatchPepalign.BstringSet,method = "hamming")/width(FragMatchPepalign.BstringSet)[1])
    library(ape)
    tree <- bionj(d)
    library(ggtree)
    FragMatchPepJudge <- pepJudgeMatrix[FragMatchPep,]
    groupInfo <- split(FragMatchPep,FragMatchPepJudge)
    if(length(groupInfo) >1){
      names(groupInfo) <- c("Non-immune","Immune")
    }else{
      names(groupInfo) <- c("Immune")
    }
    tree <- groupOTU(tree,groupInfo)
    p <- ggtree(tree,layout = "rectangular",size=0.1) +
      geom_tiplab(aes(col = group),align = TRUE) + geom_tippoint(aes(color = group))
    data = tidy_msa(FragMatchPepalign.BstringSet, 1, 20)
    ptree <- p + geom_facet(geom = geom_msa, data =data,  panel = 'msa', color= "Chemistry_AA") +xlim_tree(1) +
      theme(legend.position = "top")
    outName1 = paste0(outPath,prefix,"TreePlot_",Frag,".pdf")
    ggsave(ptree,filename = outName1,device = "pdf",
           width = 14,height = 7)

    ####figure with logo
    strings <- c()
    for(i in 1:length(FragMatchPepalign.BstringSet)){
      strings <- c(strings,as.character(FragMatchPepalign@unmasked[[i]]))
    }
    outName2 = paste0(outPath,prefix,"LogoPlot_",Frag,".pdf")
    pdf(outName2)
    p1 <- ggmsa(FragMatchPepalign.BstringSet,1,20, color = "Chemistry_AA")
    library(ggseqlogo)
    p2 <- axis_canvas(p1, axis='x')+geom_logo(strings, 'probability')
    pp <- insert_xaxis_grob(p1, p2,position="top", grid::unit(.05, "null"))
    print(ggdraw(pp))
    dev.off()
  }
  for(i in PepFragChoose){
    pepFragVis.singleFrag(i)
  }

}


