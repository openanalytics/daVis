#' Get top genes with highest significance and top genes with highest logFC
#' 
#' @param input long-format table
#' @param topGenes number of top genes to extract
#' @param featuresIdVar column with feature identifier (must be unique)
#' @return data.frame with top genes
#' @author Laure Cougnaud, Katarzyna Gorczak
getTopGenes <- function(
  input, 
  topGenes, 
  featuresIdVar = character()
) {
  
  seqTopGenes <- if(topGenes > 0) seq_len(topGenes)	else 0
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  
  requireNamespace("plyr")
  tbl <- plyr::ddply(input, "coef", function(x) {
    dt <- x
    if (any(grepl(".compCoef", colnames(dt)))) {
      x <- rbind(
        x[, c(featuresIdVar, "P.Value", "logFC")], 
        setNames(x[, c(featuresIdVar, "P.Value.compCoef", "logFC.compCoef")], 
                 c(featuresIdVar, "P.Value", "logFC"))
      )
    }
    x <- x[rank(x[, "P.Value"], ties.method = "random") %in% seqTopGenes | 
             (nrow(x) - rank(abs(x[, "logFC"]), ties.method = "random") + 
                1) %in% seqTopGenes, ]
    
    dt[dt[, featuresIdVar] %in% x[, featuresIdVar], ]
  })
  
  tbl
}



#' Create table for genes of interest
#' 
#' @param input long-format table with top tables for coefficients of interest
#' @param featuresIdVar column name with feature ids, 'ENTREZID' by default.
#' @param genesToHighlight string with identifiers of the genes to highlight.
#' The gene identifiers should correspond to the variable specified in 
#' \code{genesToHighlightVar}.
#' @param genesToHighlightVar column of the data to which the genes to 
#' highlight should be map
#' @param genesToHighlightThresholdPValue numeric, if specified (1 by default)
#' keep among the genes to highlight, the genes which have a raw p-value lower 
#' (strict) than this threshold for at least one of the coefficient considered
#' @param genesToHighlightThresholdLogFC numeric, if specified (NULL by default)
#' keep among the genes to highlight, the genes which have an absolute log FC 
#' higher (strict) than this threshold for at least one of the coefficient 
#' considered
#' @return long-format table
#' @author Laure Cougnaud, Kirsten Van Hoorde
createTopTableGenesOfInterest <- function(
  input, 
  featuresIdVar = character(),
  genesToHighlight, 
  genesToHighlightVar, 
  genesToHighlightThresholdPValue = 1, 
  genesToHighlightThresholdLogFC = NULL
) {
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  
  tblGenesOfInterest <- input[which(
    input[, featuresIdVar] %in% genesToHighlight), ]
  if (nrow(tblGenesOfInterest) == 0) {
    warning("No features of interest to highlight. Make sure that ",
            "'genesToHighlight' are correct feature identifiers.")
  }
  
  outputList <- filterGenesOfInterest(
    input = tblGenesOfInterest, 
    featuresIdVar = featuresIdVar,
    genesToHighlightVar = genesToHighlightVar,
    genesToHighlightThresholdPValue = genesToHighlightThresholdPValue, 
    genesToHighlightThresholdLogFC = genesToHighlightThresholdLogFC)
  
  outputList
}


#' Filter genes of interest
#' 
#' @param input long-format table
#' @param featuresIdVar column name with feature ids, 'ENTREZID' by default.
#' @param genesToHighlightVar column name with genes to highlight
#' @param genesToHighlightThresholdPValue numeric, if specified 
#' keep among the genes to highlight, the genes which have a raw p-value lower 
#' (strict) than this threshold for at least one of the coefficient considered
#' @param genesToHighlightThresholdLogFC numeric, if specified 
#' keep among the genes to highlight, the genes which have an absolute log FC 
#' higher (strict) than this threshold for at least one of the coefficient 
#' considered
#' @return long-format table
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
filterGenesOfInterest <- function(
    input, featuresIdVar = character(), genesToHighlightVar, 
    genesToHighlightThresholdPValue = 1, genesToHighlightThresholdLogFC = NULL
) {
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  
  pval <- genesToHighlightThresholdPValue
  logfc <- genesToHighlightThresholdLogFC
  requireNamespace("plyr")
  tblGenesOfInterest <- plyr::ddply(
    input, genesToHighlightVar, function(x) {
      dt <- x
      if (any(grepl(".compCoef", colnames(dt)))) 
        x <- rbind(
          x[, c(featuresIdVar, "P.Value", "logFC")], 
          setNames(
            x[, c(featuresIdVar, "P.Value.compCoef", "logFC.compCoef")], 
            c(featuresIdVar, "P.Value", "logFC"))
        )
      if (!is.null(pval) & !is.null(logfc)){ 
        if(sum(x[, "P.Value"] < pval) > 0 && 
           sum(abs(x[, "logFC"]) > logfc) > 0) dt
      } else {
        if (!is.null(pval)) {
          if(sum(x[, "P.Value"] < pval) > 0) dt
        } else if (!is.null(logfc)) if (sum(abs(x[, "logFC"]) > logfc) > 0) dt
      }
    }    
  )
  if (nrow(tblGenesOfInterest) > 0) {
    includeTableGenesOfInterest <- TRUE
  } else {
    text <- paste0(
      "No features of interest have",
      if(!is.null(pval)) paste0(" a raw p-value smaller than the specified ",
                                "threshold (i.e. ", pval, ")"),
      if(!is.null(pval) & !is.null(logfc)) " or", 
      if(!is.null(logfc)) paste0(" an absolute logFC higher than the ", 
                                 "specified threshold (i.e. ", logfc, ")"), 
      " for any of the coefficients.")
    warning(text)
    includeTableGenesOfInterest <- FALSE
  }
  
  list(
    topTableOutputGenesOfInterest = tblGenesOfInterest, 
    includeTableGenesOfInterest = includeTableGenesOfInterest
  )
}

#' Calculate correlation 
#' 
#' @param input long-format table
#' @return text with information from \code{cor.test}
#' @author Katarzyna Gorczak
calcCorrelation <- function(input) {
  plyr::ddply(input, "comparison", function(x){
    corTest <- stats::cor.test(x[, "logFC"], x[, "logFC.compCoef"])
    paste0("Pearson cor=", round(corTest$estimate, 2), " (", "p=", 
           format(corTest$p.value, scientific = TRUE, digits = 2), ")")
    
  })
}

#' Color-blind palette
#' 
#' @param n number of colors to return
#' @param grey whether '#999999' color should be included
#' @return n colors 
colorBlindPalette <- function(n, grey = TRUE) {
  cbp <- c(if(grey) "#999999",
           c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
             "#0072B2", "#D55E00", "#CC79A7"))
  requireNamespace("grDevices")
  palette <- grDevices::colorRampPalette(cbp)
  colours <- palette(n)
  colours
}

#' Get number of significant genes for coefficient
#' 
#' @param input a list of top tables.
#' @param fdr threshold for adjusted p-value.
#' @param logFCrange numeric, upper and lower bound for logFC
#' @return named numeric vector with number of significant genes 
#' @author Katarzyna Gorczak
getNumberOfSignificantGenes <- function(
  input, 
  logFCrange = NULL,
  fdr = 0.05
) {
  
  tbl <- arrangeTopTables(
    input = input, 
    commonFeatures = FALSE, 
    fdr = fdr,
    logFCrange = logFCrange,
    output = "list"
  )
  
  vapply(tbl, function(x) length(which(x[, "adj.P.Val"] <= fdr)), numeric(1))
}

#' Get number of significant genes for coefficient
#' 
#' @param input a list of top tables.
#' @param fdr threshold for adjusted p-value.
#' @param dir direction, select either up- or down-regulated genes 
#' ('up' - genes with logFC larger than 0, 
#' 'down' - genes with logFC smaller than 0).
#' @param logFCrange numeric, upper and lower bound for logFC
#' @return named numeric vector with number of up- or down-regulated genes 
#' @author Katarzyna Gorczak
getNumberOfRegulatedGenes <- function(
  input,
  fdr = 0.05,
  logFCrange = NULL, 
  dir = c("up", "down")
) {
  
  dir <- match.arg(dir)
  
  tbl <- arrangeTopTables(
    input = input, 
    commonFeatures = FALSE, 
    fdr = fdr,
    logFCrange = logFCrange,
    dir = ifelse(dir == "up", "pos", "neg"),
    output = "list"
  )
  vapply(tbl, function(x) nrow(x), numeric(1))
}

#' check if the aesthetic is fixed (e.g. color, shape, size 'palette')
#' @param typeVar name of variable for aesthetic
#' @param valVar fixed value of variable of aesthetic
#' @return logical, if TRUE the element is fixed
#' @author Laure Cougnaud
setFixElement <- function(typeVar, valVar) {
  (length(typeVar) == 0) & (length(valVar) > 0)
}

#' check if manual aesthetic should be set
#' 
#' This is the case only if \code{typeVar} and \code{valVar} are specified,
#' and if the variable is not numeric or integer (doesn't work with ggplot2)
#' @param x data.frame with \code{typeVar}
#' @param typeVar name of variable for aesthetic
#' @param valVar fixed value of variable of aesthetic
#' @return logical, if TRUE the manual scale should be set
#' @author Laure Cougnaud
setManualScale <- function(x, typeVar, valVar) {
  (length(typeVar)) > 0 & (length(valVar) > 0) & 
  !class(x[, typeVar]) %in% c("numeric", "integer")
}

#' extend manual scale values if required
#' @param x data.frame with \code{nameVar}
#' @param valVar fixed value of variable of aesthetic
#' @param nameVar name of variable for aesthetic
#' @return vector of manual scales
#' @author Laure Cougnaud
formatManualScale <- function(x, valVar, nameVar) {
  if(is.null(names(valVar))){
    values <- rep(valVar, length.out = nlevels(factor(x[, nameVar])))
    names(values) <- NULL #cannot provide named argument for colors
  }else{values <- valVar}
  return(values)
}

#' check if manual aesthetic for the gradient should be set
#' 
#' This is the case only if \code{typeVar} and \code{valVar} are specified,
#' and if the variable is numeric or integer 
#' @param x data.frame with \code{typeVar}
#' @param typeVar name of variable for aesthetic
#' @param valVar fixed value of variable of aesthetic
#' @return logical, if TRUE the manual scale should be set
#' @author Katarzyna Gorczak
setGradientScale <- function(x, typeVar, valVar) {
  (length(typeVar)) > 0 & (length(valVar) > 0) & 
  class(x[, typeVar]) %in% c("numeric", "integer")
}

#' remove space in variable name
#' @param var string
#' @return String with formatted variable
formatVariableSpace <- function(var) {
  if(grepl(" ", var))	paste0("`", var, "`")	else var
}

#' shorten string length
#' @param x data.frame with \code{var}
#' @param var column name 
#' @param length max number of characters in each string
#' @return string of length \code{length}
#' @author Katarzyna Gorczak
#' @return Character vector with shortened variable.
formatVariableLength <- function(x, var, length) {
  substr(x[, var], 1, length)
}

#' check if variable is not numeric (or integer)
#' @param x data.frame
#' @param typeVar column name in \code{x}
#' @return Logical, TRUE if the variable is not a numeric or integer
setCategoricalElement <- function(x, typeVar)	{
  (length(typeVar)) > 0 & !class(x[, typeVar]) %in% c("numeric", "integer")
}

#' Create colors for feature labels
#' @inheritParams daLogRatioPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param order logical whether to order features 
#' (used only when plot is 'logRatio')
#' @param plot plot 
#' @return named color palette for feature labels 
#' @author Katarzyna Gorczak
getFeatureColor <- function(
    topTableOutput, features, featuresIdVar, featuresColor, order = FALSE,
    plot = c("logRatio", "heatmap", "waterfall")
) {
  
  plot <- match.arg(plot)
  providedFeatures <- features
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  
  features <- unique(topTableOutput[, featuresIdVar])
  # order features based on provided ones -> necessary when 
  # featuresOrder = NULL ('daLogRatioPlot')
  if (length(providedFeatures) > 0) {
    if (length(providedFeatures) != length(features))
      warning("Some 'features' are not present in 'input'.")
    providedFeatures <- providedFeatures[which(providedFeatures %in% features)] 
  } else providedFeatures <- features
  if (!order & plot == "logRatio") 
    features <- features[match(providedFeatures, features)]
  nFeatures <- length(features)
  if (length(featuresColor) == 1) 
    featuresColor <- setNames(rep_len(featuresColor, nFeatures), features)
  if (!is.null(names(featuresColor))) {
    if (!all(features %in% names(featuresColor))) {
      stop("'featuresColor' is not specified for all 'features'.")
    }
  } else {
    if (length(featuresColor) == nFeatures) {
      names(featuresColor) <- providedFeatures
    } else stop("'featuresColor' should be of length 1 or ", 
                "the same length as 'features'.")
  }
  featuresColor <- featuresColor[match(features, names(featuresColor))]
  if (!order & plot == "logRatio") featuresColor <- rev(featuresColor)  
  
  featuresColor
}

#' Create colors for coefficient labels
#' 
#' @inheritParams daHeatmapLogFC
#' @return named character (named color palette)
#' @author Katarzyna Gorczak
getCoefColor <- function(coef, coefLabel, coefColor) {
  
  nCoefs <- length(coef)
  if (!is.list(coefLabel)) {
    if (length(coefColor) == 1) 
      coefColor <- setNames(rep_len(coefColor, nCoefs), coef)
    if (!is.null(names(coefColor))) {
      if (!all(coef %in% names(coefColor))) {
        stop("'coefColor' is not specified for all coefficients ('coef').")
      } else coefColor <- coefColor[match(coef, names(coefColor))]
    } else {
      if (length(coefColor) == nCoefs) {
        names(coefColor) <- coef
      } else stop("'coefColor' should be of length 1 or the", 
                  " same length as 'coef'.")
    }
  } else {
    if (length(coefColor) != 1) 
      warning("'coefColor' is not supported when 'coefLabel' ", 
              "is provided as a list.")
    coefColor <- setNames(rep_len("black", nCoefs), coef)
  }
  
  coefColor
}

#' Combine data.frames by columns allowing for different number of rows.
#' @param list list with data frames which should be combined by columns
#' @param featuresIdVar columns with feature ids by which data frames should
#' be combined
#' @param sort logical, whether the rows should be sorted
#' @return data.frame
#' @author Katarzyna Gorczak
cbindFill <- function(list, featuresIdVar, sort = FALSE) {
  # use suffixes ncol(x)+1 and ncol(x)+2 to avoid duplicated col names
  # use merge with all = TRUE to keep features which are not present
  Reduce(function(x, y) merge(
    x, y, by = featuresIdVar, all = TRUE, sort = sort,
    suffixes = c(ncol(x)+1, ncol(x)+2)), list)
}

#' Concatenate feature variables 
#' @param tbl long-form top table for all coefficients
#' @param vars column names to concatenate
#' @param nChar maximum number of characters to truncate the 
#' feature labels to
#' @author Katarzyna Gorczak
#' @return character vector with concatenated with ' | ' columns
concatenateVars <- function(tbl, vars, nChar = NULL) {
  
  tbl[, vars] <- lapply(tbl[, vars, drop = FALSE], as.character)
  vars <- do.call(paste, c(tbl[, vars, drop=FALSE], sep=" | "))
  if(!is.null(nChar)) substr(vars, 1, nChar) else vars
  
}

#' Process table for visualizations. 
#' 
#' Some can be missing across different coefficients (either different subsets
#' of to tables are provided or provided 'features' are not present in each top
#' table). The table will be processed for each coefficient: 
#' \itemize{
#' \item{1: missing features will be added (if a missing feature in coef X is
#' present in coef Y, this label will be used)}
#' \item{2: feature missing in all coefs will obviously not be plotted}
#' \item{3: the same 'featuresVar' can happen for different 'featuresIdVar' - 
#' add index to duplicates as 'featuresVar' is used for plotting}
#' }
#' @inheritParams daLogRatioPlot
#' @param tbl logn-form table (top tables for specified coefs)
#' @param order character vector with feature ids in a specific order
#' @importFrom stats na.omit
#' @author Katarzyna Gorczak
#' @return data.frame with processed features
processFeatures <- function(tbl, coef, featuresIdVar, order) {
  
  allNames <- tbl[!duplicated(tbl[, c(featuresIdVar, "featuresVar")]), c(
    featuresIdVar, "featuresVar")]
  allNames <- setNames(allNames[, "featuresVar"], allNames[, featuresIdVar])
  
  if (any(table(tbl[, "featuresVar"]) > length(coef))) {
    warning("Truncated feature annotations are not unique. ", 
            "A suffix '_' is added to make them unique.")
    duplicates <- TRUE
  } else duplicates <- FALSE
  
  compVars <- grep("^comparison([[:digit:]]{1,})?$", colnames(tbl), value=TRUE)
  tbl <- plyr::ddply(tbl, "coef", function(x) {
    x <- x[match(order, x[, featuresIdVar]), ]
    x[, featuresIdVar] <- order # replace NA (if any) with correct feature id
    rownames(x) <- NULL
    x <- x[rev(rownames(x)), ]
    idx <- which(is.na(x$coef)) 
    # replace NA (if any) with correct coef, coefLabel and featureVar
    if (length(idx) > 0) {
      idFeat <- paste(x[idx, featuresIdVar, drop = TRUE], collapse = ", ")
      warning(length(idx), " features (", idFeat, 
              ") are not present for the coefficient ", na.omit(unique(x$coef)))
      x[idx, c("coef", compVars)] <- na.omit(unique(x[, c("coef", compVars)]))
      x[idx, "featuresVar"] <- allNames[match(
        x[idx, featuresIdVar], names(allNames))]
      if (featuresIdVar != "featureID") 
        x[idx, "featureID"] <- x[idx, featuresIdVar]
    }
    if (duplicates) x[, "featuresVar"] <- makeElementsUnique(x[, "featuresVar"])
    x
  })
  
  if (any(table(tbl$featuresVar) != length(coef)) ||
      any(table(tbl[, featuresIdVar]) != length(coef)))
    stop("Comparison data is not correct.")
  
  tbl
  
}

#' Add facets to ggplot object
#' @param facetNCol number of columns in facets, by default the function 
#' \code{n2mfrow} is used
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object
#' @param scales see \code{facet_wrap} or \code{facet_nested_wrap}
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggh4x facet_nested_wrap
#' @return ggplot object
#' @author Laure Cougnaud, Katarzyna Gorczak
facet <- function(g, topTableOutput, facetNCol, scales = "fixed") {
  facetVars <- grep("^comparison([[:digit:]]{1,})?$", 
                    colnames(topTableOutput), value = TRUE)
  fm <- stats::as.formula(paste("~", paste(facetVars, collapse = " + ")))
  if (is.null(facetNCol)){
    datFacet <- unique(topTableOutput[, facetVars, drop = FALSE])
    requireNamespace("grDevices")
    facetNCol <- grDevices::n2mfrow(nrow(datFacet))[2]
  }
  fctFacet <- if(length(facetVars) > 1){
    facet_nested_wrap
  }else{
    facet_wrap
  }
  g <- g + do.call(fctFacet, 
                   list(facets = fm, ncol = facetNCol, scales = scales))
  
  g
}

#' Add point labels to ggplot object
#' @param topTableOutputTopGenes data.frame with top genes
#' @param includeTableGenesOfInterest whether to label \code{genesToHighlight}
#' @param topTableOutputGenesOfInterest data.frame with \code{genesToHighlight}
#' @param colorVar string with column name containing variable used 
#' for coloring
#' @param topGenesCex cex for topGenes labels
#' @param g ggplot object
#' @param ... Extra parameters passed to \code{geom_text_repel} to customize
#' the position of the gene labels
#' @importFrom ggplot2 aes
#' @importFrom rlang sym
#' @return ggplot object
#' @author Laure Cougnaud, Katarzyna Gorczak
labelTopGenes <- function(
    g, includeTableGenesOfInterest, topTableOutputTopGenes, 
    topTableOutputGenesOfInterest, colorVar, topGenesCex, ...
) {
  if (includeTableGenesOfInterest) {
    topTableOutputTopGenes <- topTableOutputTopGenes[
      ! topTableOutputTopGenes[, "topGenesVar"] %in% 
        topTableOutputGenesOfInterest[, "topGenesVar"], ]
  } 
  
  mainArgsRepel <- c(
    list(label = 'topGenesVar'),
    if(length(colorVar) > 0)  list(color = formatVariableSpace(colorVar))
  )
  mainArgsRepel <- lapply(mainArgsRepel, sym)
  aesStringRepel <- c(
    list(data = topTableOutputTopGenes, 
         mapping = do.call(aes, mainArgsRepel),
         size = topGenesCex,
         show.legend = FALSE),
    list(...)
  )
  
  g <- g + do.call(getFromNamespace(
    "geom_text_repel", ns = "ggrepel"), aesStringRepel)
  
  g
}


#' Add layer with genes of interest to ggplot object
#' @inheritParams daVolcanoPlot
#' @param topTableOutputGenesOfInterest data.frame with \code{genesToHighlight}
#' @param g ggplot object
#' @importFrom ggplot2 geom_point aes
#' @importFrom rlang .data
#' @return ggplot object
#' @author Laure Cougnaud, Katarzyna Gorczak
labelGenesOfInterest <- function(
    g, topTableOutputGenesOfInterest, typePlot, color = "red",
    sizeVar = character(), size = if(length(sizeVar) > 0) numeric()  else 2, 
    genesToHighlightCex, ...
) {
  # to avoid a warning from plotly, use geom_point
  # Warning message:
  # In geom2trace.default(dots[[1L]][[1L]],dots[[2L]][[1L]],dots[[3L]][[1L]]):
  #   geom_GeomTextRepel() has yet to be implemented in plotly.
  
  # geom_point does not have 'label' param 
  # instead of using geom_text, use geom_text_repel
  argsPointGoI <- c(
    list(data = topTableOutputGenesOfInterest),
    # added because the size does not change otherwise
    # (without it, the size only changes if sizeVar specified)
    if(setFixElement(sizeVar, size))  list(size = size),
    list(color = color, pch = 21, show.legend = FALSE)
  )
  
  g <- g + do.call(geom_point, argsPointGoI)
  
  if(typePlot == "static") {
    requireNamespace("ggrepel")
    g <- g +
      ggrepel::geom_text_repel(data = topTableOutputGenesOfInterest, 
                               aes(label = .data[['genesToHighlightVar']]),
                               color = color,
                               size = genesToHighlightCex,
                               show.legend = FALSE, ...)
  }
  
  g
}

