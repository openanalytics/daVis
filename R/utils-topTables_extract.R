#' Extract top tables for coefficients of interest
#' @param coef character vector with coefficients
#' @param coefLabel character vector with coefficient labels or list with
#' such vector (in case of multiple labels)
#' @param features character vector with features to subset \code{input}. 
#' If NULL, top table for the top \code{n} features (for the first coefficient) 
#' is extracted.
#' @param n top features to extract. If set to Inf, all features are extracted
#' @param columns columns to extract 
#' @param errorBars logical whether to add standard errors. 
#' Allowed only for \code{input} of class \code{"MArrayLM"} (or list of those).
#' @param mean logical (FALSE by default) whether to extract column with 
#' averaged expression or logCPM ("AveExpr" or "logCPM") depending on the 
#' \code{input}
#' @param stat logical (FALSE by default) whether to extract column with 
#' statistic: 't' or 'F'
#' @param featuresIdVar column name with unique feature ids (empty by default)
#' @param hoverText logical whether to add point label for interactive plot
#'  (FALSE by default)
#' @param text (optional) String with name of a column, or function to extract
#' a text from the columns of a top table, to display as text.\cr
#' See the available columns in the top table in section:
#' \link[=daVis-common-doc]{Top table format}.
#' @param output output can be one of "table" or "list"
#' @return data.frame with top tables or list of those 
#' (if \code{output} is 'list') for coefficients of interest
#' @inheritParams daVis-common-args
#' @importFrom plyr rbind.fill
#' @author Katarzyna Gorczak, Laure Cougnaud
extractTopTables <- function(
  input, coef, coefLabel = NULL, features = NULL, n = 20, columns, 
  errorBars = FALSE, mean = FALSE, stat = FALSE, featuresIdVar = character(),
  hoverText = FALSE, text = NULL, output = c("table", "list")
) {
  
  output <- match.arg(output)
  checkInput(input = input)
  checkCoef(input = input, coef = coef)
  coefLabel <- getCoefLabel(coef = coef, coefLabel = coefLabel)

  if (is.null(features)) features <- getTopFeatures(
    input = input, coef = coef[1], featuresIdVar = featuresIdVar, n = n)
  
  dfList <- lapply(seq_len(length(coef)), function(i) {
    
    coefI <- coef[i]
    topTableCoef <- extractTTcoef(input, coefI, errorBars, text)
    columns <- extractColsOfInterest(
      mean, topTableCoef, columns, stat, errorBars, coefI, text)
    checkColumns(
      input = topTableCoef, featuresIdVar = featuresIdVar, coef = coefI, 
      cols = c(if(length(featuresIdVar) > 0) featuresIdVar, columns)
    )
    topTableCoef <- subsetFeatures(
      input = topTableCoef, features = features, featuresIdVar = featuresIdVar)
    df <- formatTTcoef(
      topTableCoef, columns, mean, coefI, coef, featuresIdVar, hoverText)

    if (!is.null(coefLabel) && is.list(coefLabel)){
      df[, paste0("comparison", seq_along(coefLabel))] <-
        lapply(coefLabel, `[`, i)
    }else{
      df[, "comparison"] <- coefLabel[i]
    }
    
    df
  })
  names(dfList) <- coef
  # sort the comparison[X] columns according to the coefficients
  dfList <- lapply(dfList, orderComparison, coefLabel)
  
  if (output == "table") {
    res <- do.call(plyr::rbind.fill, dfList)
    rownames(res) <- NULL
  } else res <- dfList
  
  return(res)
}

#' Add column with hover text to top table 
#' 
#' @param input top table
#' @param columns columns to extract from top table for hover info
#' @return top table with additional 'hoverText' column
#' @author Katarzyna Gorczak
addHoverText <- function(input, columns) {
  
  df <- input
  df[, "hoverText"] <- paste0(
    "logFC: ", round(df[, "logFC"], digits = 4), "<br>",
    "p-value: ", round(df[, "P.Value"], digits = 4), "<br>",
    "adj. p-value: ", round(df[, "adj.P.Val"], digits = 4), "<br>",
    apply(
      vapply(unique(setdiff(
        columns, c("logFC", "P.Value", "adj.P.Val"))), function(iVar) {
        if (is.numeric(df[, iVar])) {
          paste0(iVar, ": ", round(df[, iVar], digits = 4), "<br>")
        } else {
          paste0(iVar, ": ", df[, iVar], "<br>")
        }
      }, character(nrow(df))), 
      1, paste, collapse = "")
  )
  
  df
}

#' Extract top tables from a model
#' @inheritParams isModel
#' @param coef coefficient of interest
#' @param se Logical (FALSE by default), should standard errors be extracted?
#' @return data.frame with top tables for coefficient of interest.\cr
#' The features are ordered in original (unsorted) order in the model.\cr
#' Standard errors, if requested are in the column: 'SE'.
#' @author Katarzyna Gorczak
getTopTableFromModel <- function(input, coef, se = FALSE) {
  
  method <- attr(class(input), "package")
  
  topTable <- switch(
    method, 
    'limma' = {
      requireNamespace("limma")
      limma::topTable(input, coef = coef, n = Inf, sort.by = "none")
    },
    
    'edgeR' = {
      requireNamespace("edgeR")
      tbl <- edgeR::topTags(input, n = Inf, sort.by = "none")
      tbl <- tbl[["table"]]
      logFCcolumns <- grep("logFC", colnames(tbl))
      logFCcolumnCoef <- grep(paste0("[.]", coef, collapse="|"), colnames(tbl))
      if (length(logFCcolumns) > 1 & length(logFCcolumns[which(
        !logFCcolumns %in% logFCcolumnCoef)]) > 0) 
        tbl <- tbl[, -logFCcolumns[which(!logFCcolumns %in% logFCcolumnCoef)]]
      if (length(logFCcolumns[which(logFCcolumns %in% logFCcolumnCoef)]) > 1) {
        colnames(tbl) <- gsub("logFC.", "", colnames(tbl))
      } else colnames(tbl)[grep("logFC", colnames(tbl))] <- "logFC"
      colnames(tbl)[which(colnames(tbl) == "PValue")] <- "P.Value"
      colnames(tbl)[which(colnames(tbl) == "FDR")] <- "adj.P.Val"
      colnames(tbl)[which(colnames(tbl) == "logCPM")] <- "AveExpr"
      tbl
    },
    
    'DESeq2' = {
      tbl <- as.data.frame(input)
      colnames(tbl)[which(colnames(tbl) == "pvalue")] <- "P.Value"
      colnames(tbl)[which(colnames(tbl) == "padj")] <- "adj.P.Val"
      colnames(tbl)[which(colnames(tbl) == "baseMean")] <- "AveExpr"
      colnames(tbl)[which(colnames(tbl) == "log2FoldChange")] <- "logFC"
      colnames(tbl)[which(colnames(tbl) == "lfcSE")] <- "se"
      tbl
    }
  )
  
  if(se){
    if("se" %in% colnames(topTable)){
      stop("Standard errors ('se' column) are already available in the",
		" top table.")
    }
    topTable[, "se"] <- getSEModel(input = input, coef = coef)
  }
  
  topTable
}

#' Get standard error from a limma model.
#'
#' The standard errors are computed, according to BioC support 
#' \href{https://support.bioconductor.org/p/70175/#70185}{question 70175}, 
#' as: \deqn{stdev.unscaled * \sqrt{s2.post}}
#' with: 
#' \itemize{
#' \item{stdev.unscaled: }{unscaled standard deviations}
#' \item{s2.post: }{posterior values for sigma^2}
#' }
#' @param input object of class \code{MArrayLM} (see \code{limma})
#' @param coef character vector with one single coefficient
#' @return If \code{input} is \code{MArrayLM}, a numeric vector with SE named
#' by feature ID, NA otherwise.
#' @author Laure Cougnaud
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' getSEModel(input = model, coef = "B.LvsP")
#' @export
getSEModel <- function(input, coef){
  
  if(length(coef) != 1)
    stop("One unique coefficient should be specified.")
    
  if(inherits(input, "MArrayLM")){
    se <- c(input$stdev.unscaled[, coef] * sqrt(input$s2.post))
  }else{
    warning("Error bars not implemented for: ", class(input), ".")
    se <- NA_real_
  }
  
  return(se)
  
}

#' Extract top table for a specific coef
#' @inheritParams extractTopTables
#' @param coefI coefficient of interest
#' @return data.frame with statistics for a specific coef
#' @author Laure Cougnaud, Katarzyna Gorczak
extractTTcoef <- function(input, coefI, errorBars, text) {
  # extract top table for the specific coefficient
  if (isModel(input)) {
    # modelI <- input
    topTableCoef <- getTopTableFromModel(input, coef = coefI, se = errorBars)
  } else if(isTopTable(input)){
    topTableCoef <- input
  }else if(is.list(input)){
    
    # which element of the list contains the coefficient?
    iModel <- getInputIdCoef(input = input, coef = coefI)
    iInput <- input[[iModel]]
    
    # extract the top table
    if(isModel(iInput)){
      # modelI <- iInput
      topTableCoef <- getTopTableFromModel(
        input = iInput, coef = coefI, se = errorBars
      )
    }else{topTableCoef <- iInput}
    
  }else{stop("Input is not of the correct format")}
  
  if(!is.null(text)){
    if(is.character(text)){
      if(!text %in% colnames(topTableCoef))
        stop("'text' is not available in the column of the top table.")
      topTableCoef[, "text"] <- topTableCoef[, text]
    }else if(is.function(text)){
      topTableCoef[, "text"] <- text(topTableCoef)
    }else stop("'text' should be a character vector or a function.")
    # columns <- c(columns, "text")
  }
  
  topTableCoef
}

#' Extract columns of interest from top table
#' @inheritParams extractTopTables
#' @param coefI coefficient of interest
#' @param topTableCoef top table for one coefficient
#' @return character vector
#' @author Laure Cougnaud, Katarzyna Gorczak
extractColsOfInterest <- function(
    mean, topTableCoef, columns, stat, errorBars, coefI, text
) {
  if (mean){
    colMeanExpr <- ifelse(
      "AveExpr" %in% colnames(topTableCoef), "AveExpr", 
      ifelse(
        "logCPM" %in% colnames(topTableCoef), "logCPM", 
        stop("'input' must contain either 'AveExpr' or 'logCPM' column.")))
    columns <- unique(c(columns, colMeanExpr))
  }
  
  if(stat){
    colStat <- intersect(c("t", "F", "stat"), colnames(topTableCoef))
    if(length(colStat) == 0)
      stop("The 'F', 't' or 'stat' column with test statistics is not", 
           " available in the top table output.")
    colStat <- colStat[1]
    columns <- unique(c(columns, colStat))
  }
  
  if(errorBars){
    if("se" %in% colnames(topTableCoef)){
      columns <- unique(c(columns, "se"))
    }else{
      warning("Standard errors ('se') are not available in the ",
              "top table of coefficient: ", coefI, ", so", 
              " error bars are not included for this coefficient.")
    }
  }
  
  if(!is.null(text))
    columns <- c(columns, "text")
  
  columns
}

#' Format top table for a specific coefficient
#' @inheritParams extractTopTables
#' @param coefI coefficient of interest
#' @param topTableCoef top table for one coefficient
#' @return data.frame
#' @author Laure Cougnaud, Katarzyna Gorczak
formatTTcoef <- function(
    topTableCoef, columns, mean, coefI, coef, featuresIdVar, hoverText
) {
  
  if (mean)
    colMeanExpr <- ifelse(
      "AveExpr" %in% colnames(topTableCoef), "AveExpr", 
      ifelse(
        "logCPM" %in% colnames(topTableCoef), "logCPM", 
        stop("'input' must contain either 'AveExpr' or 'logCPM' column.")))
  
  # format top table
  df <- data.frame(
    topTableCoef[, columns, drop = FALSE], stringsAsFactors = FALSE)
  
  if (mean)  colnames(df)[grep(colMeanExpr, colnames(df))] <- "mean"
  
  # and coefficient
  df[, "coef"] <- factor(coefI, levels = coef)
  
  # if featuresIdVar is not present in input (list of topTables) or 
  # input is a model and does not contain 'genes' slot, then topTable(input) 
  # does not contain gene annotation and
  # a column 'featureID' is created to store the row names of input 
  # these are unique feature ids
  if (length(featuresIdVar) == 0) columns <- c(columns, "featureID")
  df[, "featureID"] <- if (length(featuresIdVar) > 0) 
    df[, featuresIdVar] else rownames(df)
  
  if (hoverText) df <- addHoverText(input = df, columns = columns)
  
  rownames(df) <- NULL
  df
}

#' Sort the comparison columns according to the coefficients
#' @inheritParams extractTopTables
#' @param df top table for one coefficient
#' @return data.frame
#' @author Laure Cougnaud, Katarzyna Gorczak
orderComparison <- function(df, coefLabel){
  if(!is.null(coefLabel) && is.list(coefLabel)){
    varsComp <- paste0("comparison", seq_along(coefLabel))
    df[, varsComp] <-
      lapply(seq_along(coefLabel), function(i) 
        factor(df[, varsComp[i]], levels = unique(coefLabel[[i]]))
      )
  }else{
    df[, "comparison"] <- factor(df[, "comparison"], levels = coefLabel)
  }
  return(df)
}