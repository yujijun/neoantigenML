#' @title Visualization of length statistics in peptides.
#' @param pepSeq a vector of sequence.
#' @return ggbar A bar plot from ggpplot.
#' @export
#' @examples
#' pepSeq <- Neodataset$wild_Peptide
#' pepLength(pepSeq = pepSeq)
#'
pepVisLength <- function(pepSeq){
  pepSeq %>% as.data.frame() %>%
    mutate(pep_length = stringr::str_length(pepSeq)) %>%
    ggplot(aes(x = pep_length)) +
    geom_bar(position = "stack") +
    ggtitle("Length statistics of all peptides") +
    theme(plot.title = element_text(hjust = 0.5,face = "bold")) ->
    ggbar
  return(ggbar)
}

pepVisImportant <- function(Important_matrix,output_path,prefix = "allimportant_"){
  data <- Important_matrix %>%
    tidyr::separate(col = FeatureID,into = c("resource","character","stat","fragmentLen"))
  data_fragmentLensum <- data %>% group_by(fragmentLen,stat) %>% summarise(count=n(),culmunative_importance=sum(Importance))
  data_fragmentLensum <- data_fragmentLensum %>%
    na.omit()

  data_statsum <- data %>% group_by(stat,fragmentLen) %>% summarise(count=n(),culmunative_importance=sum(Importance))
  data_statsum <- data_statsum %>%
    na.omit()

  ##绘制data_fragmentLensum的barplot
  library(ggplot2)
  plot2=ggplot(data=data_fragmentLensum, aes(x=fragmentLen,y=culmunative_importance)) + geom_bar(stat="identity",width = 0.5) + labs(x="fragmentLen",y="Culmunative importance") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.6))
  output_name <- paste0(output_path,prefix,"fragmentlen_boxplot.pdf")
  pdf(output_name,width = 10,height = 10)
  print(plot2)
  dev.off()

  ##绘制data_statsum的barplot
  plot4=ggplot(data=data_statsum, aes(x=stat,y=culmunative_importance)) + geom_bar(stat="identity",width = 0.5) + labs(x="stat",y="Culmunative importance") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.6))
  output_name <- paste0(output_path,prefix,"stat_boxplot.pdf")
  pdf(output_name,width = 10,height = 10)
  print(plot4)
  dev.off()

}

