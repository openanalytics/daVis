#' Get coefficient label
#'
#' @param coef character vector with coefficients
#' @param coefLabel character vector of coefficient labels or a 
#' function to transform the coefficients
#' @return labels of the coefficients; either a vector or a list 
#' with nested labels
getCoefLabel <- function(coef, coefLabel) {
  
  extractLabels <- function(coef, coefLabel) {
    if (!is.list(coefLabel)) {
      if (is.function(coefLabel)) {
        coefLabel <- vapply(coef, coefLabel, character(1))
      } else {
        if (is.null(coefLabel)) {
          coefLabel <- coef
        } else {
          if (length(coef) != length(coefLabel)) 
            stop("The length of 'coef' must be the same as length of ",
				"'coefLabel'.")
        }
      }
    } else  coefLabel <- lapply(coefLabel, extractLabels, coef = coef)
    coefLabel
  }
  
  coefLabel <- extractLabels(coef = coef, coefLabel = coefLabel)
  
  coefLabelCheck <- t(do.call(
    rbind.data.frame, if(is.list(coefLabel)) coefLabel else list(coefLabel)))
  # check if duplicated labels
  if(any(duplicated(coefLabelCheck)))
    stop("The specified coefficient labels are not unique.")
  
  coefLabel
  
}

#' Get coefficients in input
#' @inheritParams isModel
#' @param coef (optional) coefficient of interest
#' @return Character vector with coefficients
getModelCoefs <- function(input, coef = NULL){
    
  method <- attr(class(input), "package")
    
  tmp <- switch(
    method, 
    'limma' = {
      colnames(input[["coefficients"]])
    },
    'edgeR' = {
      cols <- colnames(input[["table"]])[
        grep("logFC", colnames(input[["table"]]))]
      if (length(cols) > 1) {
        cols <- gsub("logFC.", "", cols)
      }
      cols
    },
    'DESeq2' = {
      cols <- grep("log2FoldChange", colnames(input), value = TRUE)
      cols
    }
  )
    
  allCoefs <- if (length(tmp) == 1 & 
                  (method == "edgeR" || method == "DESeq2"))  coef  else  tmp

  return(allCoefs)
  
}

#' Check coefficients
#' @inheritParams daVis-common-args
#' @param coef character vector with coefficients
#' @return no returned value, error if any of the coefficients not present
checkCoef <- function(input, coef = NULL) {
  
  if (is.null(coef)) stop("Please specify 'coef'.")
  
  if(isModel(input)){
    allCoefs <- getModelCoefs(input, coef = coef)
    
  }else  if(isTopTable(input)) {
    # cannot be checked
    allCoefs <- coef
  }else{
    if(is.list(input)){
      isModelInput <- vapply(input, isModel, logical(1))
      coefsModel <- lapply(input[isModelInput], getModelCoefs)
      coefsTopTable <- names(input[!isModelInput])
      allCoefs <- Reduce(union, c(coefsModel, list(coefsTopTable)))
    }
  }
  
  if(any(!coef %in% allCoefs))
    stop("At least one of the coefficients is not present in 'input', ", 
         "Please specify correct 'coef'.")
  
}

#' Get ID of the input \code{list} matching to the specific coefficient
#' @param input input list
#' @param coef character vector of length 1 with coefficient
#' @return numeric of length 1 specifying which element of the list
#' contains the coefficient
getInputIdCoef <- function(input, coef){
  
  # is coef in top table...
  isCoefInTopTable <- (coef %in% names(input))
  
  # ... or in a model?
  isInputModel <- which(vapply(input, isModel, logical(1)))
  coefsModel <- lapply(input[isInputModel], getModelCoefs, coef = coef)
  isCoefInModel <- (coef %in% unlist(coefsModel))
  
  # extract results from the top table
  if(isCoefInTopTable){
    
    if(isCoefInModel){
      warning(
        "Coefficient: ", coef, " is available in both a top table",
        " and a model of the 'input' list, the top table is considered."
      )
    }
    iModel <- match(coef, names(input))
    
    # or from the model
  }else if(isCoefInModel){
    
    iModel <- isInputModel[
      which(vapply(coefsModel, function(coefs) coef %in% coefs, logical(1)))
    ]
    if(length(iModel) > 1){
      warning("Coefficient: ", coef, " is available in multiple elements",
              " of the 'input' list (", toString(iModel), "),",
              " the first one is considered.")
      iModel <- iModel[1]
    }
    
  }else{
    stop("Coefficient: ", coef, " is not available in the 'input' list.")
  }
  
  return(iModel)
  
}
