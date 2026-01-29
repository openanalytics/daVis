#' Check input fields
#' @inheritParams daVis-common-args
#' @param featuresIdVar column name with unique feature identifiers
#' @param coef coefficient name
#' @param cols columns with feature annotation
#' @param error Logical, if TRUE (by default) an error is returned
#' if required columns are not present
#' @return (invisibly) columns not available in the input and
#' error if required columns not present (if 
#' \code{error} is TRUE)
#' @author Katarzyna Gorczak
checkColumns <- function(
  input, 
  featuresIdVar, 
  coef, 
  cols,
  error = TRUE
) {
  
  if (length(featuresIdVar) > 0 && length(featuresIdVar) != 1) {
    stop("'featuresIdVar' must indicate only one column.")
  }
  
  if (isModel(input) || isTopTable(input)) {
    
    if(isModel(input)){
      input <- getTopTableFromModel(input = input, coef = coef)
    }
    colsModel <- colnames(input)
    notAvailableCols <- cols[which(!cols %in% colsModel)]
    
  } else {
    
    if(is.list(input)){
      notAvailableCols <- lapply(input, checkColumns, 
        featuresIdVar = featuresIdVar, 
        coef = coef, cols = cols,
        error = FALSE
      )
      notAvailableCols <- Reduce(union, notAvailableCols)
    }

  }
  
  if (error && length(notAvailableCols) != 0) {
    idCols <- unique(paste(notAvailableCols, collapse = "', '"))
    stop("The column(s): '", idCols, "'",
         " must be present in 'input'.")
  }
  
  return(invisible(notAvailableCols))
  
}
