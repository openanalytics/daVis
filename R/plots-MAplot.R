#' MA plot 
#' 
#' This is a function to create a MA plot. When several coefficients are used, 
#' multiple plots side by side are returned. 
#' MAplot visualizes mean expression versus logFC for each gene.
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' existing labels.
#' @param featuresIdVar column name with feature ids.
#' @param featuresVar column name with feature ids to label.
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
#' specified in \code{featuresIdVar} or row names of \code{input}
#' @param genesToHighlightVar column name for the labels of genesToHighlight.
#' Same as \code{featuresIdVar} by default.
#' @param genesToHighlightCex cex for genesToHighlight
#' @param direction logical whether to color significant up- and 
#' down-regulated genes
#' @param fdr threshold considered for direction (up- or down- 
#' significant features), 0.05 by default.
#' @param color colors for points indicating direction (should be three colors:
#' significant up- and down-regulated and non-significant)
#' @param alpha transparency level for the points, 0.4 by default.
#' @param sizeVar column name used for the size, empty by default.
#' @param size character or factor with specified size(s) (cex) for the points, 
#' replicated if needed. This is used only if sizeVar is empty. 
#' By default: '2.5' if sizeVar is not specified and 
#' default ggplot size(s) otherwise
#' @param typePlot plot type can be one of "static" or "interactive".
#' @param ... Extra parameters passed to \code{geom_text_repel} to customize
#' the position of the gene labels.
#' @inheritParams createDataMAplot
#' @import ggplot2
#' @importFrom grDevices n2mfrow colorRampPalette pdf png
#' @importFrom stats setNames
#' @return ggplot object or a list with ggplot object and top genes highlighted 
#' in the MA plot (top 10 genes with highest significance and/or highest logFC)
#' \code{featuresVar} with names \code{featuresIdVar} 
#' (if \code{returnTopGenes} is set to TRUE)
#' @author Katarzyna Gorczak
#' @importFrom utils getFromNamespace
#' @importFrom plyr dlply
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' coefs <- c("B.LvsP", "L.LvsP")
#' 
#' # Simple MA plot
#' daMAplot(input = model, coef = coefs[1])
#' 
#' # Color by significant direction
#' daMAplot(input = model, coef = coefs, coefLabel = c("A", "B"),
#' direction = TRUE, color = c("steelblue", "firebrick", "grey"), facetNCol = 2)
#' 
#' # see vignette for other examples
#' 
#' @export
daMAplot <- function(
  input, 
  coef = NULL, 
  coefLabel = NULL,
  featuresIdVar = character(),
  featuresVar = featuresIdVar, 
  logFCrange = NULL,
  xlab = "log2 mean expression",
  ylab = "logFC", 
  axesCex = 1, 
  axesTitleCex = 1, 
  title = NULL, 
  titleCex = 1,
  legendPosition = c("right", "bottom", "none"),
  legendTitleCex = 1, 
  legendCex = 0.8,
  facetCex = 1.1,
  facetColor = "black",
  facetNCol = grDevices::n2mfrow(length(coef))[2], 
  topGenes = 0,
  topGenesVar = featuresIdVar,
  topGenesCex = 2.5,
  returnTopGenes = FALSE,
  genesToHighlight = NULL, 
  genesToHighlightVar = featuresIdVar,
  genesToHighlightCex = 2.5,
  direction = FALSE,
  fdr = 0.05,
  color = if(direction) character()  else "grey",
  alpha = 0.5,
  sizeVar = character(),
  size = if(length(sizeVar) > 0) numeric()  else 2, 
  typePlot = c("static", "interactive"),
  ...
) {
  
  typePlot <- match.arg(typePlot)
  legendPosition <- match.arg(legendPosition)
  
  if (length(genesToHighlightVar) > 1) 
    stop("'genesToHighlightVar' must indicate only one column name.")
  
  if (topGenes > 0 & length(topGenesVar) > 1) 
    stop("'topGenesVar' must indicate only one column name.")
  
  if (typePlot == "interactive") featuresVar <- unique(c(
    featuresIdVar, featuresVar, topGenesVar, genesToHighlightVar))
  
  multiplePlot <- length(coef) > 1
  
  tblList <- createDataMAplot(
    input = input, coef = coef, coefLabel = coefLabel, typePlot = typePlot,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar, 
    genesToHighlight = genesToHighlight, topGenes = topGenes, sizeVar = sizeVar,
    topGenesVar = topGenesVar, genesToHighlightVar = genesToHighlightVar,
    logFCrange = logFCrange, fdr = fdr, direction = direction
  )
  
  topTableOutput <- tblList[["topTableOutput"]]
  topTableOutputTopGenes <- tblList[["topTableOutputTopGenes"]]
  includeTableGenesOfInterest <- tblList[["includeTableGenesOfInterest"]]
  
  if(includeTableGenesOfInterest) {
    topTableOutputGenesOfInterest <- tblList[["topTableOutputGenesOfInterest"]]
  } else {
    topTableOutputGenesOfInterest <- NULL
  }
  
  ggPlot <- callMAplot(
    topTableOutput = topTableOutput,
    topTableOutputTopGenes = topTableOutputTopGenes, 
    includeTableGenesOfInterest = includeTableGenesOfInterest,
    topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
    xlab = xlab, ylab = ylab, title = title, size = size, sizeVar = sizeVar,
    axesCex = axesCex, axesTitleCex = axesTitleCex, titleCex = titleCex, 
    facetCex = facetCex, facetColor = facetColor, facetNCol = facetNCol, 
    topGenesCex = topGenesCex, genesToHighlightCex = genesToHighlightCex,
    color = color, alpha = alpha, legendPosition = legendPosition, 
    legendCex = legendCex, legendTitleCex = legendTitleCex, 
    direction = direction, multiplePlot = multiplePlot, typePlot = typePlot, ...
  )
  
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
    
    list("MAplot" = ggPlot, "topGenesEachTP" = topGenesEachTP)
  } else ggPlot
}

#' Create table for MA plot
#' @inheritParams extractTopTables
#' @inheritParams daMAplot
#' @importFrom plyr ddply
#' @return data.frame 
#' @author Katarzyna Gorczak
createDataMAplot <- function(
  input, coef, coefLabel = NULL, featuresIdVar, featuresVar, topGenes,
  topGenesVar, genesToHighlight, genesToHighlightVar, logFCrange, direction,
  sizeVar, fdr, typePlot
) {
  
  columns <- c(featuresIdVar, featuresVar, topGenesVar, genesToHighlightVar, 
               "logFC", "P.Value", "adj.P.Val", sizeVar)
  
  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel, features = NULL, 
    featuresIdVar = featuresIdVar, n = Inf, columns = columns, mean = TRUE,
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), output = "list"
  )
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(topGenesVar) == 0) topGenesVar <- featuresIdVar
  if (length(featuresVar) == 0) featuresVar <- featuresIdVar
  
  tbl <- arrangeTopTables(
    input = tblList, logFCrange = logFCrange, featuresIdVar = featuresIdVar, 
    output = "table"
  )
  
  tbl[, "topGenesVar"] <- tbl[, topGenesVar]
  
  tbl <- ddply(tbl, "coef", function(x) {
    x$direction <- "NS" 
    x$direction[x[, "adj.P.Val"] < fdr & x[, "logFC"] > 0] <- "Up"
    x$direction[x[, "adj.P.Val"] < fdr & x[, "logFC"] < 0] <- "Down"
    x
  })
  tbl[, "direction"] <- factor(
    tbl[, "direction"], levels = c("Up", "Down", "NS"))

  # Subset with top genes: top 10 genes (by default) with highest significance, 
  # and top 10 genes with highest logFC
  tblTopGenes <- getTopGenes(
    input = tbl, topGenes = topGenes, featuresIdVar = featuresIdVar)
  
  # Select genes of interest (if specified)
  if(!is.null(genesToHighlight)) {
    includeTableGenesOfInterest <- TRUE
    if (length(genesToHighlightVar) == 0) genesToHighlightVar <- featuresIdVar
    tbl[, "genesToHighlightVar"] <- tbl[, genesToHighlightVar]
  } else {includeTableGenesOfInterest <- FALSE}
  
  if (includeTableGenesOfInterest) {
    tblGenesOfInterest <- tbl[
		which(tbl[, featuresIdVar] %in% genesToHighlight), 
	]
    if (nrow(tblGenesOfInterest) == 0) {
      warning("No features of interest to highlight. Make sure that ", 
              "'genesToHighlight' are correct feature identifiers.")
      includeTableGenesOfInterest <- FALSE
    } else includeTableGenesOfInterest <- TRUE
  }
  
  output <- list("topTableOutput" = tbl, 
                 "topTableOutputTopGenes" = tblTopGenes, 
                 "includeTableGenesOfInterest" = includeTableGenesOfInterest)
  
  if (includeTableGenesOfInterest) {
    output <- c(output, 
                list("topTableOutputGenesOfInterest" = tblGenesOfInterest))
  }	
  
  output
  
}

#' Create ggplot object with MA plot
#' 
#' @inheritParams daMAplot
#' @param topTableOutput combined topTables for all coefficients
#' @param topTableOutputTopGenes data.frame with top genes
#' @param includeTableGenesOfInterest whether to label \code{genesToHighlight}
#' @param topTableOutputGenesOfInterest data.frame with \code{genesToHighlight}
#' @param multiplePlot whether to use facet_wrap on coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
callMAplot <- function(
    topTableOutput, topTableOutputTopGenes, includeTableGenesOfInterest, title, 
    topTableOutputGenesOfInterest, xlab, ylab, axesCex, axesTitleCex, titleCex, 
    facetCex, facetColor, facetNCol, genesToHighlightCex, topGenesCex, color,
    direction, alpha, sizeVar, size, legendPosition, legendTitleCex, legendCex,
    typePlot, multiplePlot, ...
) {
  
  g <- mainMA(
    typePlot = typePlot, direction = direction, color = color, size = size,
    sizeVar = sizeVar, topTableOutput = topTableOutput, alpha = alpha)
  g <- formatAesMA(
    g = g, topTableOutput = topTableOutput, direction = direction, size = size,
    color = color, sizeVar = sizeVar
  )
  
  if (includeTableGenesOfInterest) { 
    g <- labelGenesOfInterest(
      g = g, topTableOutputGenesOfInterest = topTableOutputGenesOfInterest,
      sizeVar = sizeVar, size = size, typePlot = typePlot, color = "#0e4209",
      genesToHighlightCex = genesToHighlightCex, ...) 
  }
  
  if (nrow(topTableOutputTopGenes) > 0 & typePlot == "static") {
    colorVar <- if(direction) "direction" else character()
    g <- labelTopGenes( 
      g = g, includeTableGenesOfInterest = includeTableGenesOfInterest, 
      topTableOutputTopGenes = topTableOutputTopGenes, colorVar = colorVar, 
      topTableOutputGenesOfInterest = topTableOutputGenesOfInterest, 
      topGenesCex = topGenesCex, ...) 
  }
  
  if (multiplePlot) 
    g <- facet(g = g, topTableOutput = topTableOutput, facetNCol = facetNCol)
  if (direction) g <- g + guides(colour = guide_legend(title = ""))
  
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, xTextSize = axesCex,
    facetLabelSize = facetCex, facetLabelColor = facetColor,yTextSize = axesCex,
    xTitle = xlab, xTitleSize = axesTitleCex, legendTextSize = legendCex,
    yTitle = ylab, yTitleSize = axesTitleCex, legendPosition = legendPosition, 
    legendTitleSize = legendTitleCex)
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
}

#' Create main plot object with MA plot
#' @inheritParams daMAplot
#' @param topTableOutput combined topTables for all coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
mainMA <- function(
    typePlot, direction, color, sizeVar, topTableOutput, alpha, size
) {
  addHoverText <- typePlot == "interactive"
  if(direction & length(color) == 0) 
    color <- c("firebrick", "dodgerblue4", "grey")
  
  mainArgs <- c(
    list(x = 'mean', y = 'logFC'),
    if(direction)           list(color = 'direction'),
    if(length(sizeVar) > 0)	list(size = formatVariableSpace(sizeVar))
  )
  
  aesString <- c(    
    list(data = topTableOutput, 
         mapping = do.call(ggplot2::aes_string, c(
           mainArgs,if(addHoverText) list(text = "hoverText"))))
  )
  
  g <- do.call(getFromNamespace("ggplot", ns = "ggplot2"), aesString)
  
  aesPoint <- c(
    list(alpha = alpha),
    if(!direction) list(color = color),
    if(setFixElement(sizeVar, size)) list(size = size),
    list(show.legend = TRUE)
  )
  
  g <- g + do.call(getFromNamespace("geom_point", ns = "ggplot2"), aesPoint)
  g <- g + geom_hline(yintercept = 0, size = 0.2, color = "red")
  
  g
}

#' Format aesthetics for the main plot with MA plot
#' @inheritParams daMAplot
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object with MA plot
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
formatAesMA <- function(g, topTableOutput, direction, color, sizeVar, size) {
  # manual specifications: custom scales
  # only if variable, values are specified and if the variable is not numeric 
  # or integer (doesn't work with ggplot2)
  setManualScaleStatic <- function(typeVar, nameVar, valVar){
    # always keep three levels for colors (up, down and NS)
    values <- rep(valVar, length.out = 3)
    names(values) <- NULL
    # add number of up- and down-regulated genes when coloring by direction
    if (length(unique(topTableOutput[, "coef"])) == 1) {
      nsign <- c(length(which(topTableOutput$direction == "Up")),
                 length(which(topTableOutput$direction == "Down")))
      nsign[3] <- nrow(topTableOutput) - sum(nsign)
      labels <- paste0(c("Up", "Down", "NS"), ": ", nsign)
    } else labels <- c("Up", "Down", "NS")
    do.call(getFromNamespace(
      paste("scale", typeVar, "manual", sep = "_"), ns = "ggplot2"),
      list(values = values, labels = labels, drop = FALSE) )
  }
  
  if (direction & setManualScale(topTableOutput, 'direction', color))	
    g <- g + setManualScaleStatic("color", 'direction', color)
  if (setManualScale(topTableOutput, sizeVar, size))
    g <- g + setManualScaleStatic("size", sizeVar, size)
  
  g
}