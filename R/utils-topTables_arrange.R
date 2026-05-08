#' Arrange top tables 
#' 
#' @param input list of top tables or long table
#' @param featuresIdVar column with feature identifier (must be unique)
#' @param logFCrange numeric vector with two values 
#' (low- and high-threshold for log-foldchange)
#' @param commonFeatures logical whether to keep the same set of features
#'  per coefficients
#' @param fdr fdr threshold
#' @param dir direction for the logFC ("pos" or "neg"; NULL by default). 
#' The features with "pos" diraction will be extracted with logFC above 0. 
#' The features with "neg" diraction will be extracted with logFC below 0. 
#' @param output whether to output list of top tables or long table 
#' @return arranged list or long data.frame with top tables for coefficients 
#' of interest
#' @importFrom plyr rbind.fill
#' @author Katarzyna Gorczak
arrangeTopTables <- function(
  input, 
  featuresIdVar, 
  logFCrange = NULL, 
  commonFeatures = FALSE, 
  fdr = NULL, 
  dir = NULL, 
  output = c("table", "list")
) {
  
  output <- match.arg(output)
  
  if (!inherits(input, "list")) input <- reshapeTable(input, timevar = "coef")
  
  # subset logFC (logFCrange)
  if (!is.null(logFCrange)) input <- filterLogFC(input, logFCrange)
  
  if (!is.null(fdr)) {
    input <- lapply(input, function(table){
      # add 'which' otherwise NA rows (if present) are included
      # NA values can be present in 'padj' when DESeq2 is used 
      # (see ?DESeq2::results)
      table[which(table[, "adj.P.Val"] < fdr), ] 
    })
  }
  
  if (!is.null(dir)) {
    input <- lapply(input, function(table){ 
      table[which(if(dir == "pos") { 
        table[, "logFC"] > 0
      } else { table[, "logFC"] < 0 }), ]
    })
  }
  
  # this check should be at the end, after all filtering (if necessary)
  if (commonFeatures) input <- filterCommonFeatures(input, featuresIdVar)
  
  if (output == "table") {
    df <- do.call(plyr::rbind.fill, input)
    rownames(df) <- NULL
  } else {
    input <- lapply(input, function(x) {rownames(x) <- NULL; x})
    df <- input
  }
  
  df
}

#' Filter top tables based on logFC threshold
#' 
#' @param input list of top tables 
#' @param logFCrange numeric vector with two values 
#' (low- and high-threshold for log-foldchange)
#' @return list of top tables with subset of features based on logFC range
#' @author Katarzyna Gorczak
filterLogFC <- function(input, logFCrange) {
  
  logFCrange <- sort(logFCrange)
  input <- lapply(input, function(x) {
    x[which(x[, "logFC"] >= logFCrange[1] & x[, "logFC"] <= logFCrange[2]), ]
  })
  
  dfNrow <- unlist(lapply(input, function(x) nrow(x)))
  if (any(dfNrow == 0)) {
    idZero <- paste(names(which(dfNrow == 0)), collapse = ", ")
    stop("There are no genes with logFC in range (", 
         logFCrange[1], "; ", logFCrange[2], ") for coef: ", idZero, ".")
  }
  
  input
}

#' Filter top tables based on common features across all tables
#' 
#' @param input list of top tables 
#' @param featuresIdVar column with feature identifier (must be unique) 
#' @return list of top tables with common features
#' @author Katarzyna Gorczak
filterCommonFeatures <- function(input, featuresIdVar) {
  
  # make sure that all the top tables contain the same set of features
  commonFeatures <- Reduce(
    intersect, lapply(input, function(x) x[, featuresIdVar]))
  input <- lapply(input, function(x) 
    x[match(commonFeatures, x[, featuresIdVar]), ])
  
  dfNrow <- unlist(lapply(input, function(x) nrow(x)))
  if (any(dfNrow == 0)) {
    idZero <- paste(names(which(dfNrow == 0)), collapse = ", ")
    stop("There are no common genes for ", idZero, ".")
  }
  
  input
  
}

#' Make elements unique
#' 
#' @param x vector with values
#' @return vector with unique values, 
#' for duplicated values a suffix '_' is added
makeElementsUnique <- function(x){
  x <- as.character(x)
  duplElements <- names(which(table(x) > 1))
  for (i in duplElements) {
    tmp <- x[x == i]
    x[x == i] <- paste(tmp, seq_along(tmp), sep = "_")
  }
  x
}

#' Reshape table
#' 
#' @param table long- or wide-format table
#' @param timevar column name in long-format table that differentiates multiple
#'  records from the same group
#' @return a list of tables 
#' @author Katarzyna Gorczak
reshapeTable <- function(table, timevar) {
  x <- lapply(unique(table[, timevar]), function(x)
    table[which(table[, timevar] == x), ])
  names(x) <- unique(table[, timevar])
  x
}
