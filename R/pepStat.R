#' @title Length Statistics
#' @param pepSeq a vector of sequence,which could get from by 
#' @return ggbar A bar plot from ggpplot
#' @export
#' @examples
#' pepSeq <- Neodataset$wild_Peptide
#' pepLength(pepSeq = pepSeq)
#' 
pepLength <- function(pepSeq){
  pepSeq %>% as.data.frame() %>% 
    mutate(pep_length = stringr::str_length(pepSeq)) %>% 
    ggplot(aes(x = pep_length)) + 
    geom_bar(position = "stack") + 
    ggtitle("Length statistics of all peptides") + 
    theme(plot.title = element_text(hjust = 0.5,face = "bold")) -> 
    ggbar 
  return(ggbar)
}

