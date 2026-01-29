#' Check if input is a model
#' 
#' @param input object of class \code{MArrayLM} (see \code{limma}), 
#' \code{DGELRT} (see \code{edgeR}) or \code{DESeqResults} (see \code{DESeq2})
#' @return logical
isModel <- function(input) {
  inherits(input, "DGELRT") || inherits(input, "MArrayLM") || 
    inherits(input, "DESeqResults")
}

#' Check if input is a topTable containing colums with logFC, 
#' p-value and adjusted p-value
#' 
#' @param input a top table (see \code{limma::topTable})
#' @return logical
isTopTable <- function(input) {
  cols <- c("logFC", "P.Value", "adj.P.Val")
  check <- is.data.frame(input) && all(cols %in% colnames(input))
  return(check)
}

#' Check input
#' @param error logical whether to return error message
#' @inheritParams daVis-common-args
#' @return (invisible) logical indicating if \code{input} passes the check,
#' and error if input is not of the desired class (if \code{error} is TRUE).
checkInput <- function(input, error = TRUE) {
  
  check <- isModel(input) || isTopTable(input)
  
  if(!check){
    if(is.list(input)){
      check <- all(vapply(input, checkInput, error = FALSE, logical(1)))
    }else{
      check <- FALSE
    }
  }
  
  if(error && !check){
    stop("'input' must be an object of class 'DGELRT', 'MArrayLM', ", 
         "'DESeqResults', top table (see limma::topTable), ", 
         "or a list of those objects.")
  }
  
  return(invisible(check))
  
}
