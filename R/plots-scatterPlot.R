#' Scatter plot
#' 
#' This is a function to create a scatter plot. When several coefficients are 
#' used, multiple plots side by side are returned. 
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' existing labels.
#' @param featuresIdVar column name with feature ids, empty by default. 
#' If specified and input is a model, featuresIdVar should be a column 
#' name in 'genes' slot.
#' @param featuresVar column name with feature ids to label, the same as 
#' 'featuresIdVar' by default.
#' @param fdr threshold considered for significance, NULL by default.
#' @param logFCrange numeric, two values (upper and lower bounds for logFC).
#' @param xlab x-axis title, NULL by default.
#' @param ylab y-axis title, NULL by default.
#' @param axesCex cex for the axis text.
#' @param axesTitleCex cex for the axis title.
#' @param title plot title, NULL by default.
#' @param titleCex cex for the plot title.
#' @param legendPosition legend position ("right", "bottom", "none"). 
#' If 'none', no legend is shown.
#' @param legendTitleCex cex for the legend title.
#' @param legendCex cex for the legend text.
#' @param facetCex cex for the facets if multiple plot are used.
#' @param facetColor color for the text of the facets
#' @param facetNCol number of columns in facets, by default the function 
#' \code{n2mfrow} is used to.
#' @param color color palette to distinguish significance groups. 
#' Four colors must be specified.
#' @param topGenes numeric, number of top genes with highest logFC or p-value 
#' to highlight in the plot for each considered coefficient, 0 by default.
#' @param topGenesVar column name with feature identifier to label topGenes
#' @param topGenesCex cex for topGenes labels
#' @param returnTopGenes logical, if TRUE (FALSE by default), return a list 
#' with the top genes highlighted in the plot.
#' @param genesToHighlight string with identifiers of the genes to highlight, 
#' NULL by default. The gene identifiers should correspond to the variable 
#' specified in \code{genesToHighlightVar}, and be contained among the column 
#' names of the output of the \code{topTable} function from \code{limma}.
#' @param genesToHighlightVar column name with the genes to highlight.
#' Same as \code{featuresIdVar} by default.
#' @param genesToHighlightCex cex for genesToHighlight
#' @param genesToHighlightThresholdPValue numeric, if specified (1 by default)
#' keep among the genes to highlight, the genes which have a raw p-value lower 
#' (strict) than this threshold for at least one of the coefficient considered.
#' @param genesToHighlightThresholdLogFC numeric, if specified (NULL by default)
#' keep among the genes to highlight, the genes which have an absolute log FC
#' higher (strict) than this threshold for at least one of the coefficient 
#' considered.
#' @param alpha transparency level for the points, 0.4 by default.
#' @param pointSize point size, 2 by default.
#' @param correlation logical whether to calculate the correlation and add 
#' text to the plot.
#' @param correlationCex cex for the text with correlation value
#' @param typePlot plot type, can be one of "static" or "interactive"
#' @param ... Extra parameters passed to \code{geom_text_repel} to customize
#' the position of the gene labels.
#' @inheritParams createDataScatterPlot
#' @return ggplot object or a list with ggplot object and top genes 
#' highlighted in the scatter plot 
#' (top 10 genes with highest significance and/or highest logFC)
#' \code{featuresVar} with names \code{featuresIdVar} 
#' (if \code{returnTopGenes} is set to TRUE)
#' @author Katarzyna Gorczak
#' @importFrom utils getFromNamespace
#' @importFrom plyr dlply
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' 
#' # Simple scatter plot
#' daScatterPlot(input = model, coef = coefs[c(1,2)])
#' 
#' # More coefficients
#' daScatterPlot(input = model, coef = coefs, coefLabel = c("A", "B", "C", "D"),
#' facetNCol = 3)
#' 
#' # LogFC range
#' daScatterPlot(input = model, coef = coefs, coefLabel = c("A", "B", "C", "D"),
#' facetNCol = 3, logFCrange = c(-2, 2))
#' 
#' # see vignette for other examples
#' 
#' @export
daScatterPlot <- function(
  input, 
  coef = NULL, 
  coefLabel = coef,
  featuresIdVar = character(), 
  featuresVar = featuresIdVar,
  fdr = 0.05,
  logFCrange = NULL,
  xlab = NULL, 
  ylab = NULL, 
  axesCex = 1, 
  axesTitleCex = 1, 
  title = NULL,
  titleCex = 1,
  legendPosition = c("right", "bottom", "none"),
  legendTitleCex = 1, 
  legendCex = 0.8,
  facetCex = 1,
  facetColor = "black",
  facetNCol = grDevices::n2mfrow(length(coef))[2], 
  color = c("gray90", "darkgoldenrod1", "dodgerblue", "darkgreen"), 
  topGenes = 0, 
  topGenesVar = featuresIdVar,
  topGenesCex = 2.5,
  returnTopGenes = FALSE,
  genesToHighlight = NULL, 
  genesToHighlightVar = featuresIdVar,
  genesToHighlightCex = 2.5,
  genesToHighlightThresholdPValue = 1,
  genesToHighlightThresholdLogFC = NULL,
  alpha = 0.4, 
  pointSize = 2,
  correlation = FALSE, 
  correlationCex = 3,
  typePlot = c("static", "interactive"),
  ...
) {
  
  # Check plot type and required columns
  typePlot <- match.arg(typePlot)
  legendPosition <- match.arg(legendPosition)
  
  if (!is.null(genesToHighlight) & length(genesToHighlightVar) > 1) 
    stop("'genesToHighlightVar' must indicate only one column name.")
  
  if (topGenes > 0 & length(topGenesVar) > 1) 
    stop("'topGenesVar' must indicate only one column name.")
  
  if (typePlot == "interactive") {
    featuresVar <- unique(
      c(featuresIdVar, featuresVar, topGenesVar, genesToHighlightVar))
  }
  
  # Check color palette (4 colors: non-sign, sign in both, 2x sign in each)
  if (length(color) != 4)  stop("'color' must contain 4 colors.")
  multiplePlot <- length(coef) > 2
  if (length(coef) < 2) stop("At least 2 'coef' must be provided.")
  
  tblList <- createDataScatterPlot(
    input = input, coef = coef, coefLabel = coefLabel,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar,
    genesToHighlight = genesToHighlight, topGenes = topGenes, 
    topGenesVar = topGenesVar, genesToHighlightVar = genesToHighlightVar, 
    genesToHighlightThresholdPValue = genesToHighlightThresholdPValue, 
    genesToHighlightThresholdLogFC = genesToHighlightThresholdLogFC,
    fdr = fdr, logFCrange = logFCrange, typePlot = typePlot
  )
  
  coefs <- tblList[["coef"]]
  topTableOutput <- tblList[["topTableOutput"]]
  topTableOutputTopGenes <- tblList[["topTableOutputTopGenes"]]
  includeTableGenesOfInterest <- tblList[["includeTableGenesOfInterest"]]
  
  if(includeTableGenesOfInterest) {
    topTableOutputGenesOfInterest <- tblList[["topTableOutputGenesOfInterest"]]
  } else topTableOutputGenesOfInterest <- NULL
  
  ggPlot <- callScatterPlot(
    topTableOutput = topTableOutput, coef = coefs,
    topTableOutputTopGenes = topTableOutputTopGenes, 
    includeTableGenesOfInterest = includeTableGenesOfInterest,
    topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
    legendPosition = legendPosition, legendTitleCex = legendTitleCex, 
    legendCex = legendCex, xlab = xlab, ylab = ylab, title = title, 
    axesCex = axesCex, axesTitleCex = axesTitleCex, titleCex = titleCex, 
    facetNCol = facetNCol, facetColor = facetColor, facetCex = facetCex, 
    topGenesCex = topGenesCex, genesToHighlightCex = genesToHighlightCex,
    alpha = alpha, pointSize = pointSize, fdr = fdr, color = color, 
    typePlot = typePlot, multiplePlot = multiplePlot,
    correlation = correlation, correlationCex = correlationCex, ...
  )
  
  if (returnTopGenes) {
    if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
    if (length(topGenesVar) == 0) topGenesVar <- "featureID"
    
    if (topGenes > 0) {
      topGenesEachTP <- dlply(topTableOutputTopGenes, "comparison", function(x){
        res <- as.character(x[, topGenesVar])
        names(res) <- x[, featuresIdVar]
        factor(res)
      })
      
    } else topGenesEachTP <- NULL
    
    list("scatterPlot" = ggPlot, "topGenesEachTP" = topGenesEachTP)
  } else ggPlot
  
}


#' Create table for the scatter plot
#' @inheritParams extractTopTables
#' @inheritParams daScatterPlot
#' @return data.frame
#' @author Katarzyna Gorczak
createDataScatterPlot <- function(
  input, 
  coef, coefLabel = NULL,
  featuresIdVar,
  featuresVar,
  fdr, 
  topGenes, 
  topGenesVar, 
  genesToHighlight, 
  genesToHighlightVar, 
  genesToHighlightThresholdPValue, 
  genesToHighlightThresholdLogFC,
  logFCrange,
  typePlot
) {
  
  columns <- unique(c(featuresIdVar, featuresVar, topGenesVar, 
                      genesToHighlightVar, "logFC", "P.Value", "adj.P.Val"))
  tblList <- extractTopTables(
    input = input, 
    coef = coef, coefLabel = coefLabel,
    columns = columns,
    features = NULL, n = Inf,
    featuresIdVar = featuresIdVar,
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), 
    output = "list"
  )
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(topGenesVar) == 0) topGenesVar <- featuresIdVar
  if (length(featuresVar) == 0) featuresVar <- featuresIdVar
  columns <- unique(c(columns, featuresIdVar))
  
  tbl <- arrangeTopTables(
    input = tblList, 
    logFCrange = logFCrange, 
    commonFeatures = TRUE, 
    featuresIdVar = featuresIdVar, 
    output = "list"
  )
  
  res <- createPairData(tbl = tbl, featuresIdVar = featuresIdVar, fdr = fdr, 
                        typePlot = typePlot, columns = columns)
  tbl <- res[["tbl"]]
  tbl[, "coef"] <- tbl[, "comparison"]
  tbl[, "x"] <- tbl[, "logFC"]
  tbl[, "y"] <- tbl[, "logFC.compCoef"]
  tbl[, "topGenesVar"] <- tbl[, topGenesVar]
  
  # Subset with top genes: top 10 genes (by default) with highest significance,
  # and top 10 genes with highest logFC
  tblTopGenes <- getTopGenes(
    input = tbl, topGenes = topGenes, featuresIdVar = featuresIdVar
  )
  
  # Select genes of interest (if specified)
  if (!is.null(genesToHighlight)) {
    includeTableGenesOfInterest <- TRUE
    if (length(genesToHighlightVar) == 0) genesToHighlightVar <- featuresIdVar
    tbl[, "genesToHighlightVar"] <- tbl[, genesToHighlightVar]
  } else {includeTableGenesOfInterest <- FALSE}
  
  if (includeTableGenesOfInterest) {
    tblGenesOfInterestList <- createTopTableGenesOfInterest(
      input = tbl, featuresIdVar = featuresIdVar,
      genesToHighlight = genesToHighlight, 
      genesToHighlightVar = genesToHighlightVar, 
      genesToHighlightThresholdPValue = genesToHighlightThresholdPValue, 
      genesToHighlightThresholdLogFC = genesToHighlightThresholdLogFC
    )
    
    tblGenesOfInterest <- 
		tblGenesOfInterestList[["topTableOutputGenesOfInterest"]]
    includeTableGenesOfInterest <- 
		tblGenesOfInterestList[["includeTableGenesOfInterest"]]
  }
  
  output <- list("coef" = res[["coef"]],
                 "topTableOutput" = tbl, 
                 "topTableOutputTopGenes" = tblTopGenes, 
                 "includeTableGenesOfInterest" = includeTableGenesOfInterest)
  
  if (includeTableGenesOfInterest) {
    tblGenesOfInterest[, "genesToHighlightVar"] <- 
		tblGenesOfInterest[, genesToHighlightVar]
    output <- c(output, 
		list("topTableOutputGenesOfInterest" = tblGenesOfInterest))
  }	
  
  output
  
}

#' Create pairwise data
#' @param tbl list of top tables
#' @param columns columns of interest from top table
#' @inheritParams daScatterPlot
#' @importFrom limma decideTests
#' @return data.frame with reference and compared contrasts
#' @author Katarzyna Gorczak
createPairData <- function(tbl, featuresIdVar, fdr, typePlot, columns) {
  
  pval <- lapply(tbl, function(x) x[, "adj.P.Val"])
  pvalTbl <- do.call(cbind, pval)
  rownames(pvalTbl) <- tbl[[1]][, featuresIdVar]
  
  decideTestsOutput <- decideTests(pvalTbl, adjust.method = "none", p.value=fdr)
  multiplePlot <- length(tbl) > 2
  
  coefLabels <- setNames(
    unname(vapply(tbl, function(x) 
      unique(as.character(x[, "coef"])), character(1))), 
    unname(vapply(tbl, function(x) 
      unique(as.character(x[, "comparison"])), character(1))) 
  )
  refCoef <- coefLabels[1]
  compCoef <- coefLabels[-1]
  dfPairs <- lapply(compCoef, function(iCoef) {
    extractPairs(
      decideTestsOutput = decideTestsOutput, refCoef = refCoef, iCoef = iCoef, 
      multiplePlot = multiplePlot, compCoef = compCoef, tbl = tbl, 
      typePlot = typePlot, columns = columns)
  })
  
  tbl <- do.call("rbind", dfPairs)
  rownames(tbl) <- NULL
  tbl[, "comparison"] <- factor(tbl[, "comparison"], levels = names(compCoef))
  
  list(tbl = tbl, coef = coefLabels)
}

#' Create ggplot object with scatter plot
#' 
#' @inheritParams daScatterPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param topTableOutputTopGenes data.frame with top genes
#' @param includeTableGenesOfInterest whether to label \code{genesToHighlight}
#' @param topTableOutputGenesOfInterest data.frame with \code{genesToHighlight}
#' @param multiplePlot logical whether to use facet on coefficients
#' @importFrom ggplot2 facet_wrap
#' @return ggplot object
#' @author Katarzyna Gorczak
callScatterPlot <- function(
    topTableOutput, topTableOutputTopGenes, includeTableGenesOfInterest, coef,
    topTableOutputGenesOfInterest, facetNCol, xlab, ylab, title, fdr, color,
    axesCex, axesTitleCex, titleCex, legendPosition, legendTitleCex, legendCex,
    facetCex, facetColor, topGenesCex, genesToHighlightCex, alpha, pointSize,
    typePlot, multiplePlot, correlation, correlationCex, ...
) {
  g <- mainSP(
    typePlot = typePlot, color = color, topTableOutput = topTableOutput, 
    xlab = xlab, coef = coef, ylab = ylab, multiplePlot = multiplePlot, 
    pointSize = pointSize, alpha = alpha)
  
  if (includeTableGenesOfInterest) 
    g <- labelGenesOfInterest(
      g = g, topTableOutputGenesOfInterest = topTableOutputGenesOfInterest,
      size = pointSize, typePlot = typePlot, 
      genesToHighlightCex = genesToHighlightCex, ...) 
  
  if (nrow(topTableOutputTopGenes) > 0 & typePlot == "static") 
    g <- labelTopGenes( 
      g = g, includeTableGenesOfInterest = includeTableGenesOfInterest, 
      topTableOutputTopGenes = topTableOutputTopGenes, colorVar = "signif", 
      topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
      topGenesCex = topGenesCex, ...) 
   
  if(multiplePlot)	
    g <- g + facet_wrap(
      facets = stats::as.formula(paste("~", "comparison")),
      ncol = if(is.null(facetNCol)) grDevices::n2mfrow(
        nlevels(factor(topTableOutput[, "comparison"])) - 1)[2] else facetNCol)
  
  xlab <- if(!is.null(xlab)) xlab else paste0("logFC ", names(coef[1]))
  ylab <- if(!is.null(ylab)) ylab else if(multiplePlot) {
    paste0("logFC compared contrast")} else {paste0("logFC ", names(coef)[2])}
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, xTextSize = axesCex,
    facetLabelSize = facetCex, facetLabelColor = facetColor,yTextSize = axesCex,
    xTitle = xlab, xTitleSize = axesTitleCex, legendTextSize = legendCex,
    yTitle = ylab, yTitleSize = axesTitleCex, legendPosition = legendPosition,
    legendTitleSize = legendTitleCex)
  
  if (correlation)
    g <- labelCorr(
      topTableOutput = topTableOutput, g = g, correlationCex = correlationCex)
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
}

#' Create main plot object with scatter plot
#' @inheritParams daScatterPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param multiplePlot logical whether to use facet on coefficients
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline 
#' @importFrom ggplot2 geom_segment scale_colour_manual guides guide_legend
#' @importFrom rlang sym .data
#' @return ggplot object
#' @author Katarzyna Gorczak
mainSP <- function(
    typePlot, color, topTableOutput, xlab, coef, ylab, multiplePlot, 
    pointSize, alpha
) {
  addHoverText <- typePlot == "interactive"
  names(color) <- levels(topTableOutput[, "signif"])
  range <- c(min(c(topTableOutput[, "x"], topTableOutput[, "y"])) - 0.2, 
             max(c(topTableOutput[, "x"], topTableOutput[, "y"])) + 0.2) 
  dfLine <- data.frame(matrix(c(3*range/4, 3*range/4), nrow = 1))
  
  mainArgs <- list(x = 'x', y = 'y', color = 'signif')
  mainArgs <- c(mainArgs, if(addHoverText) list(text = "hoverText"))
  mainArgs <- lapply(mainArgs, sym)
  
  aesString <- list(data = topTableOutput, mapping = do.call(aes, mainArgs))
  g <- do.call(ggplot, aesString)
  
  g <- g + 
    geom_point(size = pointSize, alpha = alpha, show.legend = TRUE) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_segment(data = dfLine, 
                 mapping = aes(x = .data[["X1"]], 
                               xend = .data[["X2"]], 
                               y = .data[["X3"]], 
                               yend = .data[["X4"]]), 
                 linewidth = 0.2, 
                 color = "grey", 
                 inherit.aes = FALSE) 
  g <- g + scale_colour_manual(values = color, drop = FALSE)
  g <- g + guides(colour = guide_legend(title = "Significance"))
  g
}

#' Add label with correlation to scatter plot
#' @inheritParams daScatterPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @return ggplot object
#' @author Katarzyna Gorczak
labelCorr <- function(topTableOutput, g, correlationCex) {
  corr <- calcCorrelation(input = topTableOutput)
  requireNamespace("ggrepel")
  g <- g +
    ggrepel::geom_text_repel(data = corr,
                             aes(x = -Inf, 
                                 y = Inf,
                                 label = .data[["V1"]]), 
                             size = correlationCex,
                             show.legend = FALSE,
                             inherit.aes = FALSE)
  g
}


#' Get top table per pair (ref and comp coef) for scatter plot
#' @inheritParams daScatterPlot
#' @param decideTestsOutput output of \code{decideTests}
#' @param refCoef coefficient on the x-axis
#' @param iCoef coefficients on the y-axis
#' @param multiplePlot whether to facet on coefs to compare
#' @param compCoef coefficients to compare with \code{refCoef}
#' @inheritParams createPairData
#' @return ggplot object
#' @author Katarzyna Gorczak
extractPairs <- function(
    decideTestsOutput, refCoef, iCoef, multiplePlot, compCoef, tbl, 
    typePlot, columns
) {
  # find significant genes common for ref 'coef' and comp 'coef', 
  # and the ones significant only in one coef
  signGenes <- decideTestsOutput[, c(refCoef, iCoef)]
  sign <- factor(
    ifelse(rowSums(signGenes) == 2, 
           "Both", 
           ifelse(signGenes[, 1] == 1 & signGenes[, 2] == 0,
                  paste0("Only ", names(refCoef)),
                  ifelse(signGenes[, 1] == 0 & signGenes[, 2] == 1,
                         ifelse(multiplePlot,
                                "Only comp",
                                paste0("Only ", names(compCoef)[
                                  which(compCoef == iCoef)])), "None"))), 
    levels = c("None", 
               paste0("Only ", names(refCoef)),
               ifelse(multiplePlot,
                      "Only comp",
                      paste0("Only ", names(compCoef)[
                        which(compCoef == iCoef)])), 
               "Both"))
  
  coefToCompareTbl <- tbl[[iCoef]]
  coefToCompareTbl <- coefToCompareTbl[, c(
    "logFC", "P.Value", "adj.P.Val", "comparison", 
    if(typePlot == "interactive") "hoverText")]
  colnames(coefToCompareTbl)[c(1,2,3)] <- paste0(
    colnames(coefToCompareTbl)[c(1,2,3)], ".compCoef")
  dat <- cbind(tbl[[refCoef]][, c(
    columns, if(typePlot == "interactive") "hoverText")], coefToCompareTbl)
  dat$signif <- sign
  
  dat
}