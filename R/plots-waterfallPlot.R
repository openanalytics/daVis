#' Waterfall plot
#' 
#' This is a function to create a waterfall plot. When several coefficients are 
#' used, multiple plots side by side are returned. 
#' @param features IDs of features to show. If NULL (by default), 
#' show top 20 features.
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' existing labels.
#' @param featuresIdVar column name with feature ids, empty by default. 
#' If specified and input is a model, featuresIdVar should be a column name 
#' in 'genes' slot.
#' @param featuresVar column name with feature ids to label, the same as 
#' 'featuresIdVar' by default.
#' @param featuresColor character vector specifying colors to use for 
#' the feature label text. Either length 1 or the same length as 
#' \code{features}, allowing particular features to be highlighted.
#' @param featuresMaxNChar numeric, maximum number of characters to truncate 
#' the feature labels to.
#' @param fdr threshold considered for significance, NULL by default.
#' @param logFCrange numeric, two values (upper and lower bounds for logFC).
#' @param xlab x-axis title, 'logFC' by default.
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
#' @param fillVar name of variable (in 'input') used for filling, 
#' empty by default
#' @param fill character or factor with specified color(s) for the boxplot 
#' inside, replicated if needed. By default: 'skyblue2' if fillVar is not 
#' specified and default ggplot palette otherwise.
#' @param colorVar name of variable (in 'input') used for coloring, 
#' empty by default
#' @param color character or factor with specified color(s) for the bar border, 
#' replicated if needed. By default: 'skyblue2' if colorVar is not specified 
#' and default ggplot palette otherwise.
#' @param alphaVar column name used for the transparency, empty by default.
#' @param alpha character or factor with specified transparency(s) for the bars,
#' replicated if needed. By default: '1' if alphaVar is not specified.
#' @param alphaRange transparency (alpha) range used in the plot, possible only
#' if the alphaVar is 'numeric' or 'integer'.
#' @param typePlot plot type, can be one of "static" or "interactive".
#' @inheritParams createDataWaterfallPlot
#' @import ggplot2
#' @importFrom utils getFromNamespace
#' @return ggplot object
#' @author Katarzyna Gorczak
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' 
#' # Simple waterfall plot
#' daWaterfallPlot(input = model, coef = "B.LvsP")
#' 
#' # see vignette for other examples
#' 
#' @export
daWaterfallPlot <- function(
    input, features = NULL, coef, coefLabel = coef, featuresIdVar = character(),
    featuresVar = featuresIdVar, featuresColor = "black", featuresMaxNChar = 60,
    fdr = NULL, logFCrange = NULL, xlab = "logFC", ylab = NULL, axesCex = 1, 
    axesTitleCex = 1.1, title = NULL, titleCex = 1.1,
    legendPosition = c("right", "bottom", "none"), legendTitleCex = 1, 
    legendCex = 0.8, facetCex = 1, facetColor = "black",
    facetNCol = grDevices::n2mfrow(length(coef))[2], fillVar = character(), 
    fill = if(length(fillVar) > 0) character()  else "skyblue2",
    colorVar = character(),
    color = if(length(colorVar) > 0) character()  else "white",
    alphaVar = character(),
    alpha = if(length(alphaVar) > 0) 1  else numeric(),
    alphaRange = numeric(), typePlot = c("static", "interactive")
) {
  
  typePlot <- match.arg(typePlot)
  legendPosition <- match.arg(legendPosition)
  
  topTableOutput <- createDataWaterfallPlot(
    input = input, coef = coef, coefLabel = coefLabel, features = features,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar, fillVar = fillVar,
    featuresMaxNChar = featuresMaxNChar, colorVar = colorVar, fdr = fdr,
    alphaVar = alphaVar, logFCrange = logFCrange, typePlot = typePlot)
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(featuresVar) == 0) featuresVar <- featuresIdVar
  if (length(coef) == 1) {
    featuresColor <- getFeatureColor(
      topTableOutput = topTableOutput, features = features, plot = "waterfall",
      featuresIdVar = featuresIdVar, featuresColor = featuresColor)
  } else {
    if (length(featuresColor) != 1) 
      warning("'featuresColor' is not supported when multiple 'coef' provided.")
    features <- unique(topTableOutput[, featuresIdVar])
    featuresColor <- setNames(rep_len("black", length(features)), features)
  }
  
  ggPlot <- callWaterfallPlot(
    topTableOutput = topTableOutput, xlab = xlab, ylab = ylab, title = title, 
    axesCex = axesCex, axesTitleCex = axesTitleCex, titleCex = titleCex, 
    fillVar = fillVar, fill = fill, colorVar = colorVar, color = color,
    alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange,
    legendPosition = legendPosition, legendCex = legendCex,
    legendTitleCex = legendTitleCex, featuresColor = featuresColor,
    facetCex = facetCex, facetColor = facetColor, facetNCol = facetNCol, 
    multiplePlot = length(coef) > 1, typePlot = typePlot)
  
  ggPlot
}


#' Create data for waterfall plot
#' @inheritParams extractTopTables
#' @inheritParams daWaterfallPlot
#' @importFrom plyr ddply
#' @return data.frame 
#' @author Katarzyna Gorczak
createDataWaterfallPlot <- function(
  input, 
  coef, coefLabel = NULL,
  features = NULL, 
  featuresIdVar, 
  featuresVar, 
  featuresMaxNChar,
  fillVar,
  colorVar,
  alphaVar,
  fdr, 
  logFCrange,
  typePlot
) {
  
  columns <- unique(c(featuresIdVar, featuresVar, "logFC", "P.Value", 
                      "adj.P.Val", fillVar, colorVar, alphaVar))

  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel, features = features, 
    columns = columns, featuresIdVar = featuresIdVar, 
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), output = "list"
  )
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(featuresVar) == 0) featuresVar <- "featureID"
  tbl <- arrangeTopTables(
    input = tblList, logFCrange = logFCrange, fdr = fdr, commonFeatures = TRUE,
    featuresIdVar = featuresIdVar, output = "table"
  )
  
  tbl[, featuresVar] <- lapply(tbl[, featuresVar, drop = FALSE], as.character)
  # concatenate feature labels (if 'featuresVar' indicates multiple col names)
  tbl[, "featuresVar"] <- do.call(
    paste, c(tbl[, featuresVar, drop=FALSE], sep = " | "))
  tbl[, "featuresVar"] <- substr(tbl[, "featuresVar"], 1, featuresMaxNChar)
  
  # sort logFC within each comparison and keep the order of the features
  tbl <- ddply(tbl, "coef", function(x) {
    x <- x[order(x[, "logFC"]), ]
    x[, "grpLabel"] <- paste0(x[, "coef"], "_", x[, "featuresVar"])
    x[, "grpLabel"] <- factor(x[, "grpLabel"], levels = unique(x[, "grpLabel"]))
    x
  })
  
  tbl
}


#' Create ggplot object with waterfall plot
#' 
#' @inheritParams daWaterfallPlot
#' @param topTableOutput combined top tables for all coefficients
#' @param multiplePlot logical whether to use facet_wrap on coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
callWaterfallPlot <- function(
  topTableOutput, featuresColor, xlab, ylab, axesCex, axesTitleCex, title, 
  titleCex, facetNCol, fillVar, fill, colorVar, color, alphaVar, alpha,
  alphaRange, legendPosition, legendCex, legendTitleCex, facetCex, facetColor,
  typePlot, multiplePlot) {
  
  g <- mainWP(
    typePlot = typePlot, colorVar = colorVar, fillVar = fillVar, color = color,
    alphaVar = alphaVar, topTableOutput = topTableOutput, fill = fill, 
    alpha = alpha)

  g <- formatAesWP(
    g = g, topTableOutput = topTableOutput, colorVar = colorVar, color = color, 
    fillVar = fillVar, fill = fill, alphaVar = alphaVar, alpha = alpha, 
    alphaRange = alphaRange)
  
  if (multiplePlot) 
    g <- facet(g = g, topTableOutput = topTableOutput, facetNCol = facetNCol, 
               scales = "free")
  
  g <- g +
    ggplot2::scale_y_discrete(breaks = topTableOutput$grpLabel, 
                              labels = topTableOutput$featuresVar) 
  
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, xTitle = xlab, 
    xTitleSize = axesTitleCex, xTextSize = axesCex, legendTextSize = legendCex,
    yTitle = ylab, yTitleSize = axesTitleCex, yTextSize = axesCex, 
    yTextColor = featuresColor, yTicks = FALSE, gridMinor = FALSE, 
    facetLabelSize = facetCex, facetLabelColor = facetColor,
    legendPosition = legendPosition, legendTitleSize = legendTitleCex)
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
  
}

#' Create main plot object with waterfall plot
#' @inheritParams daWaterfallPlot
#' @param topTableOutput combined topTables for all coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
mainWP <- function(
    typePlot, colorVar, fillVar, alphaVar, topTableOutput, fill, color, alpha
) {
  addHoverText <- typePlot == "interactive"
  
  mainArgs <- c(
    list(x = 'logFC', y = 'grpLabel'),
    if(length(colorVar) > 0)	  list(color = formatVariableSpace(colorVar)),
    if(length(fillVar) > 0)	    list(fill = formatVariableSpace(fillVar)),
    if(length(alphaVar) > 0)	  list(alpha = formatVariableSpace(alphaVar))
  )
  aesString <- c(    
    list(data = topTableOutput, 
         mapping = do.call(ggplot2::aes_string, c(
           mainArgs, if(addHoverText) list(text = "hoverText"))))
  )
  
  g <- do.call(getFromNamespace("ggplot", ns = "ggplot2"), aesString)
  
  aesBarArgs <- c(
    list(stat = "identity", width = 0.7),
    if(setFixElement(fillVar, fill))            list(fill = fill),
    if(setFixElement(colorVar, color))          list(color = color),
    if(setFixElement(alphaVar, alpha))          list(alpha = alpha),
    if(setCategoricalElement(topTableOutput, fillVar))  list(position = "dodge")
  )
  g <- g + do.call(getFromNamespace("geom_bar", ns = "ggplot2"), aesBarArgs)
  g
}

#' Format aesthetics for the main plot with waterfall plot
#' @inheritParams daWaterfallPlot
#' @param topTableOutput combined topTables for all coefficients
#' @param g ggplot object with volcano plot
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
formatAesWP <- function(
    topTableOutput, colorVar, color, g, fillVar, fill, alphaVar, alpha, 
    alphaRange) {
  
  setManualScaleStatic <- function(typeVar, nameVar, valVar){
    values <- formatManualScale(x = topTableOutput, valVar, nameVar)
    do.call(getFromNamespace(
      paste("scale", typeVar, "manual", sep = "_"), ns = "ggplot2"),
      list(values = values) )
  }
  
  if (setManualScale(topTableOutput, colorVar, color))	
    g <- g + setManualScaleStatic("color", colorVar, color)
  
  if (setManualScale(topTableOutput, fillVar, fill))	
    g <- g + setManualScaleStatic("fill", fillVar, fill)
  
  if (setGradientScale(topTableOutput, colorVar, color))	
    g <- g + do.call(getFromNamespace("scale_color_gradientn", ns = "ggplot2"), 
                     list(colors = color))
  
  if (setGradientScale(topTableOutput, fillVar, fill))	
    g <- g + do.call(getFromNamespace("scale_fill_gradientn", ns = "ggplot2"), 
                     list(colors = fill))
  
  if (setManualScale(topTableOutput, alphaVar, alpha))
    g <- g + setManualScaleStatic("alpha", alphaVar, alpha)
  
  # custom transparency range, works only if alpha var is numeric, or integer
  if (length(alphaVar) > 0 && class(topTableOutput[, alphaVar]) %in% 
      c("numeric", "integer") & length(alphaRange) > 0)
    g <- g + ggplot2::scale_alpha(
      range = alphaRange,
      guide = ggplot2::guide_legend(override.aes = list(fill = "black")))
  
  g
}