#' Volcano plot 
#' 
#' This is a function to create a volcano plot. When several coefficients are
#'  used, multiple plots side by side are returned. 
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' existing labels.
#' @param featuresIdVar column name with feature ids (should be unique). 
#' If not specified, row names are used. \code{featuresIdVar} is used to 
#' label \code{topGenes} unless \code{topGenesVar} is specified. 
#' @param featuresVar column name with feature ids to label; used for extra 
#' labels if \code{typePlot} is 'interactive'.
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
#' @param colorVar (optional) string with column name containing variable used 
#' for coloring. If not specified, coloring is based on the adjusted p-values.
#' @param color color palette, only used if \code{colorVar} is specified.
#' @param shapeVar (optional) string with column name containing variable used 
#' for shaping. If not specified, only points are used.
#' @param shape shape palette, only used if \code{shapeVar} is specified.
#' @param alphaVar column name used for the transparency, empty by default.
#' @param alpha character or factor with specified transparency(s) for the 
#' points, replicated if needed. By default: '1' if alphaVar is not specified.
#' @param alphaRange transparency (alpha) range used in the plot, possible only
#' if the alphaVar is 'numeric' or 'integer'.
#' @param sizeVar column name used for the size, empty by default.
#' @param size character or factor with specified size(s) (cex) for the points, 
#' replicated if needed. This is used only if sizeVar is empty. 
#' By default: '2.5' if sizeVar is not specified and 
#' default ggplot size(s) otherwise
#' @param sizeRange size (cex) range used in the plot, 
#' possible only if the sizeVar is 'numeric' or 'integer'
#' @param additionalThresholdsAdjPValue numeric, additional adjusted p-values 
#' thresholds to use for the coloring of the points, and indicated 
#' in the legend.
#' @param typePlot plot type can be one of "static" or "interactive".
#' @param ... Extra parameters passed to \code{geom_text_repel} to customize
#' the position of the gene labels.
#' @inheritParams createDataVolcanoPlot
#' @import ggplot2
#' @importFrom grDevices n2mfrow colorRampPalette pdf png
#' @importFrom stats setNames
#' @return ggplot object or list with ggplot object and top genes highlighted 
#' in the volcano plot (top 10 genes with highest significance and/or highest 
#' logFC) \code{topGenesVar} with names \code{featuresIdVar}
#' (if \code{returnTopGenes} is set to TRUE)
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
#' @importFrom utils getFromNamespace
#' @importFrom plyr dlply
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' 
#' # Simple volcano plot 
#' daVolcanoPlot(input = model, coef = "B.LvsP")
#' 
#' # Specify logFC range
#' daVolcanoPlot(input = model, coef = "B.LvsP", logFCrange = c(-2, 2))
#' 
#' # Customized aesthetics
#' daVolcanoPlot(input = model, coef = c("B.LvsP", "L.LvsP"),
#' coefLabel = c("A", "B"), facetNCol = 2, colorVar = "adj.P.Val")
#'  
#' # Customized gene annotation
#' model$genes$group <- rep(c("gr A", "gr B", "gr C"), each = 5323)
#' daVolcanoPlot(input = model, coef = "B.LvsP", colorVar = "group",
#' color = setNames(c("orange", "red", "blue"), c("gr A", "gr B", "gr C")))
#' 
#' # Facet by variable(s)
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' coefsLabel <- list(
#' sub("(.+)\\.(.+)", "\\2", coefs),
#' sub("(.+)\\.(.+)", "\\1", coefs)
#' )
#' daVolcanoPlot(input = model, coef = coefs, coefLabel = coefsLabel, 
#' facetNCol = 4, colorVar = "adj.P.Val")
#' 
#' # see vignette for other examples
#' 
#' @export
daVolcanoPlot <- function(
  input, 
  coef = NULL, 
  coefLabel = NULL,
  featuresIdVar = character(),
  featuresVar = featuresIdVar, 
  fdr = 0.05,
  logFCrange = NULL,
  xlab = "logFC",
  ylab = "-log10(p-value)", 
  axesCex = 1, 
  axesTitleCex = 1.1, 
  title = NULL, 
  titleCex = 1.1,
  legendPosition = c("right", "bottom", "none"),
  legendTitleCex = 1, 
  legendCex = 0.8,
  facetCex = 1,
  facetColor = "black",
  facetNCol = grDevices::n2mfrow(length(coef))[2], 
  topGenes = 0, 
  topGenesVar = featuresIdVar,
  topGenesCex = 2.5,
  returnTopGenes = FALSE,
  genesToHighlight = NULL, 
  genesToHighlightVar = featuresIdVar,
  genesToHighlightCex = 2.5,
  genesToHighlightThresholdPValue = 1,
  genesToHighlightThresholdLogFC = NULL,
  colorVar = character(),
  color = if(length(colorVar) > 0) character()  else "black",
  shapeVar = character(),
  shape = if(length(shapeVar) > 0) numeric()  else 19,
  alphaVar = character(),
  alpha = if(length(alphaVar) > 0) numeric()  else 0.4,
  alphaRange = numeric(),
  sizeVar = character(),
  size = if(length(sizeVar) > 0) numeric()  else 2, 
  sizeRange = numeric(),
  additionalThresholdsAdjPValue = NULL,
  typePlot = c("static", "interactive"),
  ...
) {
  
  typePlot <- match.arg(typePlot)
  legendPosition <- match.arg(legendPosition)

  if (!is.null(genesToHighlight) & length(genesToHighlightVar) > 1) 
    stop("'genesToHighlightVar' must indicate only one column name.")
  
  if (topGenes > 0 & length(topGenesVar) > 1) 
    stop("'topGenesVar' must indicate only one column name.")
  
  if (typePlot == "interactive") featuresVar <- unique(
    c(featuresIdVar, featuresVar, topGenesVar, genesToHighlightVar))
  
  nCoefs <- length(coef)
  multiplePlot <- nCoefs > 1
  
  tblList <- createDataVolcanoPlot(
    input = input, coef = coef, coefLabel = coefLabel,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar, 
    topGenes = topGenes, genesToHighlight = genesToHighlight, 
    topGenesVar = topGenesVar, genesToHighlightVar = genesToHighlightVar, 
    genesToHighlightThresholdPValue = genesToHighlightThresholdPValue, 
    genesToHighlightThresholdLogFC = genesToHighlightThresholdLogFC,
    additionalThresholdsAdjPValue = additionalThresholdsAdjPValue,
    colorVar = colorVar, shapeVar = shapeVar, alphaVar = alphaVar,
    sizeVar = sizeVar, logFCrange = logFCrange, fdr = fdr, typePlot = typePlot
  )
  
  topTableOutput <- tblList[["topTableOutput"]]
  topTableOutputTopGenes <- tblList[["topTableOutputTopGenes"]]
  includeTableGenesOfInterest <- tblList[["includeTableGenesOfInterest"]]
  
  if(includeTableGenesOfInterest) {
    topTableOutputGenesOfInterest <- tblList[["topTableOutputGenesOfInterest"]]
  } else topTableOutputGenesOfInterest <- NULL

  # build palette; if colorVar = 'adj.P.Val, change color for adj. p-values
  if (length(colorVar) > 0 && colorVar == "adj.P.Val") {
    adjPValAvailableLevels <- levels(topTableOutput[, "adj.P.ValFct"])
    adjPValAvailableLevelsNum <- as.numeric(
      sub(".+,(.+)(\\)|\\])", "\\1", adjPValAvailableLevels))
    colorRampPaletteUsed <- c(
      colorRampPalette(c("dodgerblue3", "skyblue1"))(
        sum(adjPValAvailableLevelsNum <= fdr)),
      colorRampPalette(c("grey55", "grey10"))(
        sum(adjPValAvailableLevelsNum > fdr))
    )
    
    if (length(color) == 0) color <- colorRampPaletteUsed
    colorVar <- "adj.P.ValFct"
  }
    
  ggPlot <- callVolcanoPlot(
    topTableOutput = topTableOutput,
    topTableOutputTopGenes = topTableOutputTopGenes, 
    includeTableGenesOfInterest = includeTableGenesOfInterest,
    topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
    xlab = xlab, ylab = ylab, title = title, fdr = fdr, 
    colorVar = colorVar, color = color, shapeVar = shapeVar, shape = shape, 
    alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange,
    sizeVar = sizeVar, size = size, sizeRange = sizeRange,
    axesCex = axesCex, axesTitleCex = axesTitleCex, titleCex = titleCex, 
    facetCex = facetCex, facetColor = facetColor, facetNCol = facetNCol,
    topGenesCex = topGenesCex, genesToHighlightCex = genesToHighlightCex,
    legendPosition = legendPosition, multiplePlot = multiplePlot, 
    legendTitleCex = legendTitleCex, legendCex = legendCex, 
	typePlot = typePlot, ...)
  
  if (returnTopGenes) {
    if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
    if (length(topGenesVar) == 0) topGenesVar <- "featureID"
    
    if (topGenes > 0) {
      topGenesEachTP <- dlply(topTableOutputTopGenes, "coef", function(x){
        res <- as.character(x[, topGenesVar])
        names(res) <- x[, featuresIdVar]
        factor(res)
      })
    } else topGenesEachTP <- NULL
    
    list("volcanoPlot" = ggPlot, "topGenesEachTP" = topGenesEachTP)
  } else ggPlot
}

#' Create table for the volcano plot
#' @inheritParams extractTopTables
#' @inheritParams daVolcanoPlot
#' @return data.frame 
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak, 
#' Michela Pasetto
createDataVolcanoPlot <- function(
  input, 
  coef, coefLabel = NULL,
  featuresIdVar, 
  featuresVar, 
  colorVar, 
  shapeVar,
  alphaVar,
  sizeVar,
  fdr, 
  additionalThresholdsAdjPValue,
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
                      genesToHighlightVar, "logFC", "P.Value", "adj.P.Val", 
                      colorVar, shapeVar, alphaVar, sizeVar))
  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel, features = NULL, 
    n = Inf, columns = columns, featuresIdVar = featuresIdVar,
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), output = "list"
  )
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(topGenesVar) == 0) topGenesVar <- featuresIdVar
  if (length(featuresVar) == 0) featuresVar <- featuresIdVar
  tbl <- arrangeTopTables(
    input = tblList, logFCrange = logFCrange, output = "table"
  )
  
  # Set p-values of 0 to the smallest non-zero normalized floating-point number
  # (otherwise issue log-transformation during plot creation) in a new column 
  # to use original P.Value to extract top genes
  minVal <- .Machine$double.xmin
  tbl[, "P.ValuePlot"] <- ifelse(
    tbl[, "P.Value"] < minVal, minVal, tbl[, "P.Value"])
  tbl[, "adj.P.ValPlot"] <- ifelse(
    tbl[, "adj.P.Val"] < minVal, minVal, tbl[, "adj.P.Val"])
  tbl[, "topGenesVar"] <- tbl[, topGenesVar]
  
  # Use factor for the adjusted p-values
  thresholdsLegendVolcanoPlot <- rev(sort(
    unique(c(0, 0.001, 0.01, 0.05, 1, fdr, additionalThresholdsAdjPValue))))
  tbl[, "adj.P.ValFct"] <- cut(tbl[, "adj.P.Val"], 
                               breaks = thresholdsLegendVolcanoPlot, 
                               right = FALSE, 
                               include.lowest = TRUE)
  
  # Subset with top genes: top 10 genes (by default) with highest significance, 
  # and top 10 genes with highest logFC
  tblTopGenes <- getTopGenes(
    input = tbl, topGenes = topGenes, featuresIdVar = featuresIdVar
  )
  
  # Select genes of interest (if specified)
  if(!is.null(genesToHighlight)) {
    includeTableGenesOfInterest <- TRUE
    if (length(genesToHighlightVar) == 0) genesToHighlightVar <- featuresIdVar
    tbl[, "genesToHighlightVar"] <- tbl[, genesToHighlightVar]
  } else {includeTableGenesOfInterest <- FALSE}
  
  if (includeTableGenesOfInterest) {
    tblGenesOfInterestList <- createTopTableGenesOfInterest(
      input = tbl, featuresIdVar = NULL, genesToHighlight = genesToHighlight, 
      genesToHighlightVar = genesToHighlightVar, 
      genesToHighlightThresholdPValue = genesToHighlightThresholdPValue, 
      genesToHighlightThresholdLogFC = genesToHighlightThresholdLogFC
    )
    
    tblGenesOfInterest <- 
		tblGenesOfInterestList[["topTableOutputGenesOfInterest"]]
    includeTableGenesOfInterest <- 
		tblGenesOfInterestList[["includeTableGenesOfInterest"]]
  }
  
  output <- list("topTableOutput" = tbl, 
                 "topTableOutputTopGenes" = tblTopGenes, 
                 "includeTableGenesOfInterest" = includeTableGenesOfInterest)
  
  if (includeTableGenesOfInterest) {
    output <- c(
      output, list("topTableOutputGenesOfInterest" = tblGenesOfInterest))
  }	
  
  output
  
}

#' Create ggplot object with volcano plot
#' @inheritParams daVolcanoPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param topTableOutputTopGenes data.frame with top genes
#' @param includeTableGenesOfInterest whether to label \code{genesToHighlight}
#' @param topTableOutputGenesOfInterest data.frame with \code{genesToHighlight}
#' @param multiplePlot whether to use facet_wrap on coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
callVolcanoPlot <- function(
  topTableOutput, topTableOutputTopGenes, includeTableGenesOfInterest,
  topTableOutputGenesOfInterest, xlab, ylab, title, fdr, colorVar, color, 
  shapeVar, shape, alphaVar, alpha, alphaRange, sizeVar, size, sizeRange,
  axesCex, axesTitleCex, titleCex, legendPosition, legendTitleCex, legendCex,
  facetCex, facetColor, facetNCol, topGenesCex, genesToHighlightCex, typePlot,
  multiplePlot, ...) {
  
  g <- mainVP(
    topTableOutput = topTableOutput, colorVar = colorVar, color = color,
    shapeVar = shapeVar, shape = shape, alphaVar = alphaVar, alpha = alpha, 
    sizeVar = sizeVar, typePlot = typePlot, size = size)
  
  g <- formatAesVP(
    g = g, topTableOutput = topTableOutput, colorVar = colorVar, color = color, 
    shapeVar = shapeVar, shape = shape, sizeVar = sizeVar, sizeRange =sizeRange,
    size = size, alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange)
  
  if (includeTableGenesOfInterest) { 
    g <- labelGenesOfInterest(
      g = g, topTableOutputGenesOfInterest = topTableOutputGenesOfInterest,
      sizeVar = sizeVar, size = size, typePlot = typePlot, 
      genesToHighlightCex = genesToHighlightCex, ...) 
  }
  
  if (nrow(topTableOutputTopGenes) > 0 & typePlot == "static") {
    g <- labelTopGenes( 
      g = g, includeTableGenesOfInterest = includeTableGenesOfInterest, 
      topTableOutputTopGenes = topTableOutputTopGenes, colorVar = colorVar, 
      topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
      topGenesCex = topGenesCex, ...) 
  }
  
  if (multiplePlot) 
    g <- facet(g = g, topTableOutput = topTableOutput, facetNCol = facetNCol)

  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, xTextSize = axesCex, 
    facetLabelSize = facetCex, facetLabelColor = facetColor,yTextSize = axesCex,
    xTitle = xlab, xTitleSize = axesTitleCex, legendTextSize = legendCex, 
    yTitle = ylab, yTitleSize = axesTitleCex, legendPosition = legendPosition,
    legendTitleSize = legendTitleCex
  )
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
  
}

#' Create main plot object with volcano plot
#' @inheritParams daVolcanoPlot
#' @param topTableOutput combined topTables for all coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
mainVP <- function(
    colorVar, shapeVar, alphaVar, sizeVar, topTableOutput, 
    typePlot, color, shape, alpha, size
) {
  # Transform p-values to -log10 scale
  transData <- scales::trans_new(
    name = "log10Reverse",
    transform = function(x) -log10(x), 
    inverse  = function(x) 10^(-x),
    breaks = function(x) scales::log_breaks(base = 10)(x),
    domain = c(-Inf, Inf))
  
  mainArgs <- c(
    list(x = 'logFC', y = 'P.ValuePlot'),
    if(length(colorVar) > 0)	  list(color = formatVariableSpace(colorVar)),
    if(length(shapeVar) > 0)	  list(shape = formatVariableSpace(shapeVar)),
    if(length(alphaVar) > 0)	  list(alpha = formatVariableSpace(alphaVar)),
    if(length(sizeVar) > 0)	    list(size = formatVariableSpace(sizeVar))
  )
  aesString <- c(    
    list(data = topTableOutput, 
         mapping = do.call(ggplot2::aes_string, c(
           mainArgs, if(typePlot == "interactive") list(text = "hoverText"))))
  )
  
  g <- do.call(getFromNamespace("ggplot", ns = "ggplot2"), aesString)
  g <- g + scale_y_continuous(trans = transData)
  
  aesPointArgs <- c(
    if(setFixElement(colorVar, color))    list(color = color),
    if(setFixElement(shapeVar, shape))    list(shape = shape),
    if(setFixElement(alphaVar, alpha))    list(alpha = alpha),
    if(setFixElement(sizeVar, size))      list(size = size),
    list(show.legend = TRUE)
  )
  g <- g + do.call(getFromNamespace("geom_point", ns = "ggplot2"), aesPointArgs)
  
  g
}

#' Format aesthetics for the main plot with volcano plot
#' @inheritParams daVolcanoPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object with volcano plot
#' @import ggplot2
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
formatAesVP <- function(
    g, topTableOutput, colorVar, color, shapeVar, shape, 
    alphaVar, alpha, alphaRange, sizeVar, size, sizeRange
) {
  # manual specifications: custom scales
  # only if variable, values are specified and if the variable is not numeric 
  # or integer (doesn't work with ggplot2)
  setManualScaleStatic <- function(typeVar, nameVar, valVar){
    values <- if(nameVar == "adj.P.ValFct") valVar else formatManualScale(
      x = topTableOutput, valVar, nameVar)
    do.call(getFromNamespace(
      paste("scale", typeVar, "manual", sep = "_"), ns = "ggplot2"),
      c(list(values = values),
        if(nameVar == "adj.P.ValFct") list(drop = FALSE)))
  }
  
  if (setManualScale(topTableOutput, colorVar, color)) 
    g <- g + setManualScaleStatic("color", colorVar, color)
  
  if (setManualScale(topTableOutput, shapeVar, shape))
    g <- g + setManualScaleStatic("shape", shapeVar, shape)
  
  if (setManualScale(topTableOutput, alphaVar, alpha))
    g <- g + setManualScaleStatic("alpha", alphaVar, alpha)
  
  if (setManualScale(topTableOutput, sizeVar, size))
    g <- g + setManualScaleStatic("size", sizeVar, size)
  
  # custom transparency range, works only if alpha var is numeric, or integer
  if (length(alphaVar) > 0 && 
      class(topTableOutput[, alphaVar]) %in% c("numeric", "integer") & 
      length(alphaRange) > 0)
    g <- g + ggplot2::scale_alpha(
      range = alphaRange,
      guide = ggplot2::guide_legend(override.aes = list(fill = "black"))
    )
  
  # custom size range, works only if size variable is numeric, or integer
  if(length(sizeVar) > 0 &&
     class(topTableOutput[, sizeVar]) %in% c("numeric", "integer") & 
     length(sizeRange) > 0)	
    g <- g + ggplot2::scale_size(range = sizeRange)
  
  if (length(colorVar) > 0 && colorVar == "adj.P.ValFct") 
    g <- g + guides(colour = guide_legend(title = "adjusted p-value"))
  
  g
}

