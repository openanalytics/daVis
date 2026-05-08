#' Plot log ratios
#'
#' This is a function to create barplot with log-foldchanges.
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or function to transform 
#' the existing labels, or list with such labels.
#' \cr If a list is specified, nested facets are used. In that case, 
#' the coefficients are ordered based on the groups defined by the coefficient 
#' labels.
#' @param features feature ids or numeric. If NULL (by default), 
#' show top 20 features according to the first \code{coef}. If numeric and 
#' \code{input} is a list, the tables are subsetted according to the row names
#'  of the table for the first \code{coef}. 
#' @param featuresIdVar column name with feature ids.
#' @param featuresVar column name with feature ids to label.
#' @param featuresOrder if not NULL, features in the graph are re-ordered, 
#' either based on similarity ('similarity') or significance ('significance') 
#' of the statistics. See section 'Feature ordering'.
#' @param featuresMaxNChar maximum number of characters to truncate the 
#' feature labels to.
#' @param featuresColor character vector specifying colors to use for 
#' the feature label text. Either length 1 or the same length as 
#' \code{features}, allowing particular features to be highlighted.
#' @param color string with colors for the bars, one per coefficient.
#' @param xlab x-axis title, or \code{NULL} to remove.\cr
#' By default, set to: 'logFC' or 'logFC (+- SE)' if error bar(s)
#' are available.
#' @param xexpand expansion factor for the x-axis 
#' (see \code{\link[ggplot2]{expansion}}. If \code{text} is specified, 
#' the x-axis is expanded by 20\% on each side.
#' @param axesCex relative size for the label of the axes.
#' @param axesTitleCex relative size for the title of the axes.
#' @param title plot title, NULL by default.
#' @param titleCex relative size for the title of the plot.
#' @param facetCex cex for the facets if multiple plot are used.
#' @param facetColor color for the text of the facets
#' @param facetNCol number of columns in facets, by default the function 
#' \code{n2mfrow} is used to.
#' @param errorBars logical whether to add error bars to the plot (+/- one
#' standard error). If \code{input} (or each element within a 
#' \code{input} list) is a:
#' \itemize{
#' \item{model: only supported for \code{MArrayLM}}
#' \item{top table: a 'se' column should be available}
#' }
#' @param typePlot plot type, can be one of "static" or "interactive".
#' @param textCex cex for the text next to bars
#' @inheritParams createDataLogRatioPlot
#' @inheritSection orderFeatures Feature ordering
#' @return ggplot object
#' @author Laure Cougnaud, Katarzyna Gorczak, Heather Turner, 
#' Aditya Bhagwat, Kirsten Van Hoorde
#' @importFrom ggplot2 expansion
#' @examples 
#' exampleData <- createExampleData(path = ".", output = c("limma", "topTable"))
#' model <- exampleData$limma
#' topTableList <- exampleData$topTable
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' 
#' # Simple log-ratio plot
#' daLogRatioPlot(input = model, coef = coefs, facetNCol = 4)
#' 
#' # Specify features and annotation
#' features <- daVis:::getTopFeatures(input = model, coef = "B.LvsP", 
#' featuresIdVar = "ENTREZID", n = 20)
#' daLogRatioPlot(input = model, coef = coefs, 
#' features = features, featuresIdVar = "ENTREZID", 
#' featuresVar = c("SYMBOL", "GENENAME"), featuresMaxNChar = 35, facetNCol = 4)
#' 
#' # Specify different set of features
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV", "A")
#' set.seed(123)
#' features <- sample(features, 20)
#' daLogRatioPlot(
#' input = list(model, A = topTableList[["B.LvsP"]][c(seq_len(6), 9, 10), ]),
#' featuresIdVar = "ENTREZID", features = features, coef = coefs, 
#' facetNCol = 5, errorBars = TRUE)
#' 
#' # Sort features as specified
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' daLogRatioPlot(input = model, featuresIdVar = "ENTREZID", 
#' features = features, coef = coefs, facetNCol = 4, errorBars = TRUE)
#' 
#' # Sort features based on similarity
#' daLogRatioPlot(input = model, featuresIdVar = "ENTREZID", 
#' features = features, coef = coefs, facetNCol = 4, errorBars = TRUE,
#' featuresOrder = "similarity")
#' 
#' # see vignette for other examples
#' 
#' @export
daLogRatioPlot <- function(
    input, coef = NULL, coefLabel = NULL, features = NULL, 
    featuresIdVar = character(), featuresVar = featuresIdVar, 
    featuresOrder = NULL, featuresMaxNChar = 50, featuresColor = "black",
    color = character(), text = NULL, textCex = 4, xlab, 
    xexpand = if(!is.null(text)){ggplot2::expansion(mult = 0.2)}, axesCex = 1, 
    axesTitleCex = 1, title = NULL, titleCex = 1, facetCex = 1, 
    facetColor = "black", facetNCol = NULL, errorBars = TRUE,
    typePlot = c("static", "interactive")
) {
  
  typePlot <- match.arg(typePlot)
  
  if (length(color) == 0){
    color <- colorBlindPalette(n = length(coef))
    names(color) <- coef 
  }
  if (length(color) != length(coef)) {
    stop("'color' must contain the same number of colors as number of 'coef'.")
  }
  
  if (isTRUE(!is.null(featuresOrder) & length(features) > 1)) {
    order <- TRUE 
  } else order <- FALSE
  
  multiplePlot <- length(coef) > 1 | is.list(coefLabel)
  
  topTableOutput <- createDataLogRatioPlot(
    input = input, features = features, coef = coef, coefLabel = coefLabel,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar, text = text,
    order = order, featuresOrder = featuresOrder, errorBars = errorBars,
    featuresMaxNChar = featuresMaxNChar, typePlot = typePlot)
  
  featuresColor <- getFeatureColor(
    topTableOutput = topTableOutput, features = features, plot = "logRatio",
    featuresIdVar = featuresIdVar, featuresColor = featuresColor, order = order
  )
  
  ggPlot <- callLogRatioPlot(
    topTableOutput = topTableOutput, textVar = if(!is.null(text)){"text"}, 
    textVarCex = textCex, xlab = xlab, xexpand = xexpand, axesCex = axesCex, 
    axesTitleCex = axesTitleCex, title = title, titleCex = titleCex,
    color = color, featuresColor = featuresColor, facetCex = facetCex, 
    facetColor = facetColor, facetNCol = facetNCol,
    multiplePlot = multiplePlot, typePlot = typePlot)
  
  ggPlot
  
}


#' Create table for log ratio plot
#' @inheritParams extractTopTables
#' @inheritParams daLogRatioPlot
#' @param order logical whether to order features
#' @return data.frame 
createDataLogRatioPlot <- function(
    input, coef, coefLabel = NULL, text = NULL, features = NULL, featuresIdVar, 
    featuresVar, featuresMaxNChar, order, featuresOrder, errorBars, typePlot
) {
  
  columns <- unique(c(featuresIdVar, featuresVar, "logFC", "P.Value", 
                      "adj.P.Val"))
  tblList <- extractTopTables(
    input = input, features = features, coef = coef, coefLabel = coefLabel,
    columns = columns, featuresIdVar = featuresIdVar, output = "list",
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), text = text,
    stat = identical(featuresOrder, "similarity"), errorBars = errorBars)
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(featuresVar) == 0) featuresVar <- "featureID"
  
  tbl <- arrangeTopTables(
    input = tblList, featuresIdVar = featuresIdVar, 
    commonFeatures = FALSE, output = "table"
  )
  
  tbl[, "featuresVar"] <- concatenateVars(
    tbl = tbl, vars = featuresVar, nChar= featuresMaxNChar)
  
  if (order){
    ord <- orderFeatures(
      input = tblList, 
      featuresOrder = featuresOrder, featuresIdVar = featuresIdVar)
  }else{
    if (length(features) > 0) {
      ord <- features 
    } else ord <- unique(tbl[, featuresIdVar])
  }
  
  tbl <- processFeatures(
    tbl = tbl, coef = coef, featuresIdVar = featuresIdVar, order = ord)
  
  tbl[, "featuresVar"] <- factor(
    tbl[, "featuresVar"], levels = unique(tbl[, "featuresVar"]))
  
  tbl
}

#' Create ggplot object with log ratio plot
#' @param textVar (optional) String with name of a column to display as text.
#' @param textVarCex cex for the text next to bars
#' @inheritParams daLogRatioPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param multiplePlot whether to use facet_wrap on coefficients
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
#' @importFrom ggplot2 ylab scale_y_continuous
#' @importFrom utils getFromNamespace
#' @importFrom ggh4x facet_nested_wrap
callLogRatioPlot <- function(
    topTableOutput, xlab, xexpand, axesCex, axesTitleCex, title, titleCex, 
    facetNCol, facetCex, facetColor, color, textVar = NULL, textVarCex,
    featuresColor, typePlot, multiplePlot
) {
  
  errorBars <- "se" %in% colnames(topTableOutput)
  if (errorBars) {
    topTableOutput[, "ymin"] <- topTableOutput[, "logFC"] - 
      topTableOutput[, "se"]
    topTableOutput[, "ymax"] <- topTableOutput[, "logFC"] + 
      topTableOutput[, "se"]
  }
  
  g <- mainLRP(
    typePlot = typePlot, color = color, topTableOutput = topTableOutput, 
    errorBars = errorBars)
  
  # (optional) text
  if(!is.null(textVar))
    g <- labelTextLRP(
      g = g, topTableOutput = topTableOutput, errorBars = errorBars, 
      textVar = textVar, textVarCex = textVarCex)
  
  if (multiplePlot) 
    g <- facet(g = g, topTableOutput = topTableOutput, facetNCol = facetNCol)
  
  if(missing(xlab)) xlab <- paste0("logFC", if(errorBars) " +- SE")
  
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, 
    facetLabelSize = facetCex, facetLabelColor = facetColor, xTextSize =axesCex,
    xTitle = if(!is.null(xlab)) "", xTitleSize = axesTitleCex, 
    yTitle = NULL, yTitleSize = axesTitleCex, yTextSize = axesCex, 
    yTextColor = featuresColor, yTicks = FALSE, gridMinor = FALSE
  )
  
  if (!is.null(xlab)) g <- g + ylab(xlab)
  
  if(!is.null(xexpand)){
    g <- g + scale_y_continuous(expand = xexpand)
  }
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
  
}


#' Order features based on similarity or significance
#' @param input list of top tables
#' @inheritParams daLogRatioPlot
#' @importFrom stats dist hclust 
#' @section Feature ordering:
#' The features are ordered based on:
#' \itemize{
#' \item{'similarity': a hierarchical clustering
#' of the (Euclidean) distances between the statistics. The statistics are for:
#' \itemize{
#' \item{limma: the t-statistic of each coefficient}
#' \item{edgeR: the (overall) F-statistic is considered}
#' }}
#' \item{'significance': decreasing average (-log10) p-values
#' across coefficients. The p-values are for:
#' \itemize{
#' \item{limma: the significance of each coefficient}
#' \item{edgeR: the (overall) p-value of the model}
#' }}
#' }
#' @return ordered features
orderFeatures <- function(input, featuresOrder, featuresIdVar) {
  
  featuresOrder <- match.arg(
    arg = featuresOrder, choices = c("similarity", "significance"), 
    several.ok = FALSE
  )
  
  switch(featuresOrder,
         similarity = {
           # get statistics
           stats <- lapply(input, function(x){
             cols <- c(featuresIdVar, "t", "F")
             cols <- intersect(cols, colnames(x))
             x[, cols]
           })
           # combine across coefficients into one data.frame
           stats <- cbindFill(list = stats, featuresIdVar = featuresIdVar)
           # sort statistics
           ord <- hclust(dist(stats[, -grep(
             featuresIdVar, colnames(stats))]))$order
           ord <- stats[, featuresIdVar][ord]
         },
         
         significance = {
           pValues <- lapply(input, function(x){
             x[, grep(paste0(featuresIdVar, "|P.Value"), colnames(x))]
           })
           pValues <- cbindFill(list = pValues, featuresIdVar = featuresIdVar)
           ord <- order(base::rowMeans(-log10(
             pValues[, -grep(featuresIdVar, colnames(pValues)), drop = FALSE]), 
             na.rm = TRUE), decreasing = TRUE)
           ord <- pValues[, featuresIdVar][ord]
         }
  )
  
  return(ord)
  
}

#' Create main plot object with log-ratio plot
#' @inheritParams daLogRatioPlot
#' @param topTableOutput combined topTables for all coefficients
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip scale_fill_manual 
#' @importFrom ggplot2 geom_hline geom_errorbar
#' @importFrom rlang sym expr syms
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
mainLRP <- function(typePlot, color, topTableOutput, errorBars) {
  addHoverText <- typePlot == "interactive"
  fillVar <- "coef"
  # for compatibility with <= 1.0.15
  if(!is.null(names(color)) && "comparison" %in% colnames(topTableOutput) && 
     !all(names(color) %in% topTableOutput[, "coef"])){
    warning("Since version > 1.0.15, the 'color' should be named with the ",
            "coefficients and not the coefficient labels.")
    fillVar <- "comparison"
  }
  mainArgs <- list(x = 'featuresVar', y = 'logFC', fill = fillVar)
  mainArgs <- c(mainArgs, if(addHoverText) list(text = "hoverText"))
  mainArgs <- lapply(mainArgs, sym)
  
  aesString <- list(data = topTableOutput, mapping = do.call(aes, mainArgs))
  g <- do.call(ggplot, aesString)
  
  ymax <- ifelse(
    errorBars, 
    max(ifelse(is.na(topTableOutput$se), topTableOutput$logFC, 
               topTableOutput$logFC + topTableOutput$se), na.rm = TRUE), 
    max(topTableOutput$logFC, na.rm = TRUE))
  ymin <- ifelse(
    errorBars, 
    min(ifelse(is.na(topTableOutput$se), topTableOutput$logFC,
               topTableOutput$logFC - topTableOutput$se), na.rm = TRUE), 
    min(topTableOutput$logFC, na.rm = TRUE))
  
  aesBarArgs <- c(
    list(stat = "identity", position = "identity", 
         show.legend = FALSE, na.rm = TRUE)
  )
  g <- g + do.call(geom_bar, aesBarArgs)
  g <- g +
    coord_flip(ylim = c(ymin, ymax)) +
    scale_fill_manual(values = color) +
    geom_hline(yintercept = 0)
  
  if (errorBars)
    g <- g + geom_errorbar(data = topTableOutput, 
                           aes(ymin = !!sym('ymin'), ymax = !!sym('ymax')),
                           color = "black")
  
  g
}

#' Add optional text to the log-ratio plot
#' @inheritParams daLogRatioPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object with volcano plot
#' @param textVar (optional) String with name of a column to display as text.
#' @param textVarCex cex for the text next to bars
#' @importFrom ggplot2 geom_text aes
#' @importFrom rlang sym
#' @return ggplot object
#' @author Laure Cougnaud, Kirsten Van Hoorde, Katarzyna Gorczak
labelTextLRP <- function(topTableOutput, errorBars, g, textVar, textVarCex) {
  
  # x-position
  topTableOutput[, "textX"] <- if (errorBars){
    logFC <- topTableOutput[, "logFC"]
    minEB <- topTableOutput[, "ymin"]
    maxEB <- topTableOutput[, "ymax"]
    ifelse(
      logFC < 0, 
      ifelse(!is.na(minEB), minEB, logFC),
      ifelse(!is.na(maxEB), maxEB, logFC)
    )
  }else{topTableOutput[, "logFC"]}
  # x-adjustment
  topTableOutput[, "textXJust"] <- ifelse(
    topTableOutput[, "logFC"] < 0, 1.1, -0.1)
  textArgs <- list(label = textVar, y = "textX", hjust = "textXJust")
  textArgs <- lapply(textArgs, sym)
  g <- g + geom_text(
    data = topTableOutput,
    mapping = aes(!!!textArgs),
    cex = textVarCex
  )
  
  g
}