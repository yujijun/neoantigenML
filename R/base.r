#### basic function ####
matrixintovector <- function(mat){
  df <- as.data.frame(mat) %>%
    rownames_to_column() %>%
    tidyr::pivot_longer(cols = setdiff(colnames(.),"rowname"),
                        names_to = "name",values_to = "value") %>%
    mutate(vector_name = paste0(rowname,"_",name)) %>%
    select(c(value,vector_name))
  vec <- df$value
  names(vec) <- df$vector_name
  return(vec)
}
listtodf <- function(datalist, name = peptides){
  if(is.matrix(datalist[[1]]) == TRUE){
    #convert matrix into vector with name
    datalist <- map(datalist,matrixintovector)
  }
  names(datalist) <- name
  df <- do.call(rbind.data.frame,datalist)
  colnames(df) <- names(datalist[[1]])
  return(df)
}
