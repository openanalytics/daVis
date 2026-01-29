#' Heatmap of log fold-changes
#' 
#' This is a function to create a heatmap that represents logFC values for 
#' several coefficients. 
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform the 
#' existing labels.
#' @param features IDs of features to show. If NULL (by default), 
#' show top 20 features.
#' @param featuresIdVar column name with feature ids.
#' @param featuresVar column name with feature ids to label, the same as 
#' \code{featuresIdVar} by default.
#' @param featuresMaxNChar maximum number of characters to truncate the 
#' feature labels to.
#' @param featuresColor character vector specifying colors to use for 
#' the feature label text. Either length 1 or the same length as 
#' \code{features}, allowing particular features to be highlighted.
#' @param xlab x-axis title, NULL by default.
#' @param ylab y-axis title, NULL by default.
#' @param axesCex cex for the axis text.
#' @param axesTitleCex cex for the axis title.
#' @param title plot title, NULL by default.
#' @param titleCex cex for the plot title.
#' @param legendTitleCex cex for the legend title.
#' @param legendCex cex for the legend text.
#' @param coefColor color palette for coef labels.
#' @param color palette for gradient to fill the heatmap.
#' @param colorNA color for missing values. 
#' @param typePlot plot type, can be one of "static" or "interactive".
#' @inheritParams createDataHeatmap
#' @return ggplot object
#' @author Kirsten Van Hoorde, Laure Cougnaud and Katarzyna Gorczak
#' @importFrom utils getFromNamespace
#' @examples 
#' exampleData <- createExampleData(path = ".", output = c("limma", "topTable"))
#' model <- exampleData$limma
#' topTableList <- exampleData$topTable
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' 
#' # Simple heatmap
#' daHeatmapLogFC(input = model, coef = coefs)
#' 
#' # Specify feature annotation
#' daHeatmapLogFC(input = model, coef = coefs, featuresIdVar = "ENTREZID", 
#' featuresVar = c("SYMBOL", "GENENAME"), featuresMaxNChar = 35)
#' 
#' # Color coefficient labels
#' daHeatmapLogFC(input = model, coef = coefs, 
#' coefLabel = c("A", "B", "C", "D"), 
#' coefColor = c("blue", "red", "blue", "red"))
#' 
#' # Specify different set of features
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV", "A")
#' daHeatmapLogFC(
#' input = list(model, A = topTableList[["B.LvsP"]][c(seq_len(6), 9, 10), ]), 
#' coef = coefs)
#' 
#' # see vignette for other examples
#' 
#' @export
daHeatmapLogFC <- function(
  input,
  coef = NULL, 
  coefLabel = NULL,
  features = NULL, 
  featuresIdVar = character(),
  featuresVar = featuresIdVar,
  featuresColor = "black",
  featuresMaxNChar = 50,
  xlab = NULL,
  ylab = NULL,
  axesCex = 1, 
  axesTitleCex = 1, 
  title = NULL,
  titleCex = 1,
  legendTitleCex = 1, 
  legendCex = 0.8,
  coefColor = "black",
  color = c("#0072B2", "white", "#D55E00"),
  colorNA = "grey",                  
  typePlot = c("static", "interactive")
) {
  
  if (length(coef) < 2) 
    stop("At least 2 'coef' must be provided.")
  
  typePlot <- match.arg(typePlot)
  
  topTableOutput <- createDataHeatmap(
    input = input, coef = coef, coefLabel = coefLabel, features = features,
    featuresIdVar = featuresIdVar, featuresVar = featuresVar, 
    featuresMaxNChar = featuresMaxNChar, typePlot = typePlot)
  
  featuresColor <- getFeatureColor(
    topTableOutput = topTableOutput, features = features, plot = "heatmap",
    featuresIdVar = featuresIdVar, featuresColor = featuresColor)
  
  coefColor <- getCoefColor(
    coef = coef, coefLabel = coefLabel, coefColor = coefColor)
  
  ggPlot <- callHeatmap(
    topTableOutput = topTableOutput, xlab = xlab, ylab = ylab, title = title, 
    featuresColor = featuresColor, coefColor = coefColor, color = color, 
    colorNA = colorNA, axesCex = axesCex, axesTitleCex = axesTitleCex, 
    titleCex = titleCex, legendTitleCex = legendTitleCex, legendCex = legendCex,
    typePlot = typePlot)
  
  ggPlot
}



#' Create table for heatmap
#' @inheritParams daVis-common-args
#' @inheritParams daHeatmapLogFC
#' @return data.frame 
#' @importFrom stats dist hclust reorder as.dendrogram order.dendrogram
#' @author Laure Cougnaud, Katarzyna Gorczak
createDataHeatmap <- function(
  input, 
  coef, coefLabel = NULL,
  features = NULL,
  featuresIdVar, 
  featuresVar,
  featuresMaxNChar,
  typePlot
) {
  
  columns <- unique(
    c(featuresIdVar, featuresVar, "logFC", "P.Value", "adj.P.Val"))
  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel,
    features = features, featuresIdVar = featuresIdVar, columns = columns,
    hoverText = ifelse(typePlot == "static", FALSE, TRUE), output = "list")
  
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  if (length(featuresVar) == 0) featuresVar <- "featureID"
  
  tbl <- arrangeTopTables(
    input = tblList, featuresIdVar = featuresIdVar, 
    commonFeatures = FALSE, output = "table")
  
  tbl[, "featuresVar"] <- concatenateVars(
    tbl = tbl, vars = featuresVar, nChar= featuresMaxNChar)
  
  # order features - hierarchical clustering on logFC
  logfc <- lapply(tblList, function(x){
    x[, grep(paste0(featuresIdVar, "|logFC"), colnames(x))]
  })
  logfc <- cbindFill(list = logfc, featuresIdVar = featuresIdVar)
  if (nrow(logfc) > 1) {
    rowm <- rowMeans(logfc[, -grep(
      featuresIdVar, colnames(logfc))], na.rm = FALSE)
    hcc <- hclust(dist(logfc[, -grep(
      featuresIdVar, colnames(logfc))]))
    dend <- reorder(as.dendrogram(hcc), rowm)
    ord <- order.dendrogram(dend)
    ord <- logfc[, featuresIdVar][ord]
  } else ord <- logfc[, featuresIdVar]
  
  tbl <- processFeatures(
    tbl = tbl, coef = coef, featuresIdVar = featuresIdVar, order = ord)
  tbl[, "featuresVar"] <- factor(
    tbl[, "featuresVar"], levels = unique(tbl[, "featuresVar"]))
  
  tbl
}

#' Create ggplot object with heatmap
#' 
#' @inheritParams daHeatmapLogFC
#' @param topTableOutput combined topTables for all coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
callHeatmap <- function(
  topTableOutput,
  xlab,
  ylab,
  title, 
  color,
  axesCex, 
  featuresColor,
  axesTitleCex,
  titleCex, 
  legendTitleCex, 
  legendCex,
  typePlot,
  coefColor,
  colorNA
) {
  
  g <- mainH(typePlot = typePlot, topTableOutput = topTableOutput, 
             color = color, colorNA = colorNA)
  
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex,
    xTitle = xlab, xTitleSize = axesTitleCex, xTextSize = axesCex, 
    yTitle = ylab, yTitleSize = axesTitleCex, yTextSize = axesCex, 
    xTextColor = coefColor, yTextColor = featuresColor, 
    yTicks = FALSE, xTicks = FALSE, panelBackground = NULL, 
    panelBorder = "white", gridMajor = FALSE,
    legendTitleSize = legendTitleCex, legendTextSize = legendCex)
  
  if (typePlot != "static") {
    requireNamespace("plotly")
    plotly::ggplotly(g, tooltip = "text")
  } else g
  
}

#' Create main plot object with heatmap
#' @inheritParams daHeatmapLogFC
#' @param topTableOutput combined topTables for all coefficients
#' @import ggplot2
#' @return ggplot object
#' @author Katarzyna Gorczak
mainH <- function(typePlot, topTableOutput, color, colorNA) {
  
  addHoverText <- typePlot == "interactive"
  # limits for legend
  requireNamespace("stats")
  val <- max(abs(range(stats::na.omit(topTableOutput[, "logFC"]))))
  logfcRange <- c(-val, val)
  sumAbs <- sum(abs(logfcRange))
  
  low <- color[1]
  mid <- color[2]
  high <- color[3]
  
  facetVars <- grep("^comparison([[:digit:]]{1,})?$", 
                    colnames(topTableOutput), value = TRUE)
  mainArgs <- c(
    list(x = paste0('interaction(', 
                    paste0(facetVars, collapse = ', ' ),', sep = "!")'),
         y = 'featuresVar', fill = 'logFC')
  )
  aesString <- c(    
    list(data = topTableOutput, 
         mapping = do.call(ggplot2::aes_string, c(
           mainArgs, if(addHoverText) list(text = "hoverText"))))
  )
  g <- do.call(getFromNamespace("ggplot", ns = "ggplot2"), aesString)
  
  g <- g + 
    geom_tile() +
    scale_fill_gradient2(
      low = low, mid = mid, high = high, 
      midpoint = 0, na.value = colorNA,
      limits = logfcRange, breaks = seq(
        round(-val, 0), round(val, 0), by = ceiling(sumAbs/6))) +
    scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = "!"))
  
  g
}
