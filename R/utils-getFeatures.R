#' Subset input with specific features
#' @param features features of interest
#' @param featuresIdVar column with unique feature identifiers
#' @inheritParams daVis-common-args
#' @return input with subset of features
#' @author Katarzyna Gorczak
subsetFeatures <- function(
  input, 
  features, 
  featuresIdVar = character()
) {
  
  if (isModel(input) || isTopTable(input)) {
  
    if (isModel(input)) {
      if ("genes" %in% names(input) & length(featuresIdVar) > 0) {
        # keep only selected features
        input <- input[match(features, input[["genes"]][, featuresIdVar]), ]
		# rm features with NA feature IDs
        input <- input[!is.na(input[["genes"]][[featuresIdVar]]), ]
        if (nrow(input) == 0)
          stop("No 'features' matching 'featuresIdVar'.")
      } else {
        features <- features[which(features %in% rownames(input))]
        input <- input[features, ]
      }
    } else {
      if (length(featuresIdVar) == 0)
        features <- rownames(input[features[which(features %in% 
                                                    rownames(input))], ])
      
      if (length(featuresIdVar) > 0) {
        input <- input[match(features, input[, featuresIdVar]), ]
        input <- input[!is.na(input[[featuresIdVar]]), ]
        if (nrow(input) == 0)
          stop("No 'features' matching 'featuresIdVar'.")
      } else {
        input <- input[features, ]
      }
    }
    
  }else{
    
    if(is.list(input)){
      input <- lapply(input, subsetFeatures,
        features = features, 
        featuresIdVar = featuresIdVar
      )
    }
  }
  
  return(input)
}


#' Get top features
#' @inheritParams daVis-common-args
#' @param coef character vector with coefficient
#' @param featuresIdVar column with unique feature identifiers
#' @param n numeric, number of top features
#' @return top n features based on their adjusted p-value
#' @author Laure Cougnaud, Katarzyna Gorczak
getTopFeatures <- function(input, coef, featuresIdVar = character(), n = 10) {
  checkInput(input = input)
  checkCoef(input = input, coef = coef)
  
  if (isModel(input)) {
    method <- attr(class(input), "package")
    features <- switch(
      method, 
      'limma' = {
        topFeaturesLimma(
          input = input, coef = coef, n = n, featuresIdVar = featuresIdVar)
      },
      'edgeR' = {
        topFeaturesEdger(input = input, n = n, featuresIdVar = featuresIdVar)
      },
      'DESeq2' = {
        topFeaturesDeseq(input = input, n = n, featuresIdVar = featuresIdVar) 
      }
    )
  } else {
    if(isTopTable(input)) {
      tbl <- input[order(input[, "adj.P.Val"]), ]
      if (n == Inf) n <- nrow(tbl)
      if (n > nrow(tbl)) n <- nrow(tbl)
      features <- if (length(featuresIdVar) > 0) {
        tbl[, featuresIdVar][seq_len(n)] 
      } else rownames(tbl)[seq_len(n)]
    } else {
      if(is.list(input)) {
        iModel <- getInputIdCoef(input = input, coef = coef)
        features <- getTopFeatures(
          input = input[[iModel]], coef = coef, 
          featuresIdVar = featuresIdVar, n = n
        )
      }
    }
  }
  if (any(duplicated(features)) & length(featuresIdVar) > 0)
    stop("'featuresIdVar' must indicate unique values.")
  if (any(is.na(features)) & length(featuresIdVar) > 0)
    stop("'featuresIdVar' cannot contain 'NA' values.")
  
  features
}

#' Extract top features from \code{limma} output
#' @inheritParams getTopFeatures
#' @return character vector with feature names
#' @author Katarzyna Gorczak
topFeaturesLimma <- function(input, coef, n, featuresIdVar) {
  requireNamespace("limma")
  tbl <- limma::topTable(input, coef = coef, n = n)
  if (length(featuresIdVar) > 0) tbl[, featuresIdVar] else rownames(tbl)
}

#' Extract top features from \code{edgeR} output
#' @inheritParams getTopFeatures
#' @return character vector with feature names
#' @author Katarzyna Gorczak
topFeaturesEdger <- function(input, n, featuresIdVar) {
  requireNamespace("edgeR")
  tbl <- edgeR::topTags(input, n = n)
  tbl <- tbl[["table"]]
  if (length(featuresIdVar) > 0) tbl[, featuresIdVar] else rownames(tbl)
}

#' Extract top features from \code{DESeq2} output
#' @inheritParams getTopFeatures
#' @return character vector with feature names
#' @author Katarzyna Gorczak
topFeaturesDeseq <- function(input, n, featuresIdVar) {
  tbl <- as.data.frame(input)
  tbl <- tbl[order(tbl[, "padj"]), ]
  if (n == Inf) n <- nrow(tbl)
  if (n > nrow(tbl)) n <- nrow(tbl)
  if (length(featuresIdVar) > 0) {
    tbl[seq_len(n), featuresIdVar] 
  } else rownames(tbl)[seq_len(n)]
}

#' Get features from input
#' @inheritParams daVis-common-args
#' @param featuresIdVar column with unique feature identifiers
#' @return character with all features from input
#' @author Katarzyna Gorczak
getFeatures <- function(
  input, 
  featuresIdVar = character()
) {
  
  if (isModel(input)) {
    features <- if (length(featuresIdVar) > 0) input[["genes"]][[
      featuresIdVar]] else rownames(input)
  }else if(isTopTable(input)){
    features <- if (length(featuresIdVar) > 0) input[[
      featuresIdVar]] else rownames(input)
  }else if(is.list(input)){
    features <- getFeatures(input = input[[1]], featuresIdVar = featuresIdVar)  
  }
  
  if (any(duplicated(features)) & length(featuresIdVar) > 0)
    stop("'featuresIdVar' must indicate unique values.")
  if (any(is.na(features)) & length(featuresIdVar) > 0)
    stop("'featuresIdVar' cannot contain 'NA' values.")
  
  return(features)
  
}
