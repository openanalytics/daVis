#' Barplot with the number of significant genes
#' 
#' This is a function to create a barplot indicating the number of significant 
#' genes in each coefficient. The number of up- and down-regulated genes is 
#' shown.
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' the existing labels.
#' @param fdr threshold considered for significance, NULL by default.
#' @param logFCrange numeric, two values (upper and lower bounds for logFC).
#' @param xlab x-axis title, NULL by default.
#' @param ylab y-axis title, NULL by default.
#' @param axesCex cex for the axis text.
#' @param axesTitleCex cex for the axis title.
#' @param title plot title, NULL by default.
#' @param titleCex cex for the plot title.
#' @param legendPosition legend position ("right", "bottom", "none"). If 'none',
#'  no legend is shown.
#' @param legendTitleCex cex for the legend title.
#' @param legendCex cex for the legend text.
#' @param color color palette, must contain two colors.
#' @param annotCex cex for the text displayed in the plot 
#' (indicating number of genes).
#' @param addPercentage logical whether to add percentage 
#' @inheritParams extractTopTables
#' @import ggplot2
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#'  
#' # Simple barplot
#' daSignificantGenesBarplot(input = model, coef = coefs)
#' 
#' # Add percentage of genes
#' daSignificantGenesBarplot(input = model, coef = coefs, 
#' coefLabel = c("A", "B", "C", "D"), addPercentage = TRUE)
#' 
#' # see vignette for other examples
#' 
#' @return ggplot object
#' @author Katarzyna Gorczak
#' @export
daSignificantGenesBarplot <- function(
  input, 
  coef, 
  coefLabel = coef,
  fdr = 0.05,
  logFCrange = NULL,
  xlab = NULL, 
  ylab = NULL, 
  axesCex = 1, 
  axesTitleCex = 1.1, 
  title = NULL,
  titleCex = 1.1,
  legendPosition = c("right", "bottom", "none"),
  legendTitleCex = 1, 
  legendCex = 0.8,
  color = c("#32a6d3", "#e52323"),
  annotCex = 3.5,
  addPercentage = FALSE
) {
  
  legendPosition <- match.arg(legendPosition)
  
  if (length(color) != 2)
    stop("'color' must contain two colors.")
  if (is.null(names(color))) names(color) <- c("down", "up")
  
  tbl <- createDataBarplot(
    input = input, 
    coef = coef, coefLabel = coefLabel, 
    fdr = fdr, logFCrange = logFCrange, addPercentage = addPercentage
  )
  
  ggPlot <- callBarplot(
    tbl = tbl,
    xlab = xlab, ylab = ylab, title = title, color = color, 
    axesCex = axesCex, annotCex = annotCex,
    axesTitleCex = axesTitleCex, addPercentage = addPercentage,
    titleCex = titleCex, legendPosition = legendPosition, 
    legendTitleCex = legendTitleCex, legendCex = legendCex
  )
  
  ggPlot
  
}

#' Create table for barplot
#' @inheritParams extractTopTables
#' @inheritParams daSignificantGenesBarplot
#' @return data.frame 
#' @author Katarzyna Gorczak
createDataBarplot <- function(
    input, coef, coefLabel = coef, fdr, logFCrange = NULL, addPercentage = FALSE
) {
  
  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel, n = Inf,
    columns = c("logFC", "adj.P.Val"), hoverText = FALSE, output = "list"
  )
  tbl <- do.call(plyr::rbind.fill, tblList)
  rownames(tblList) <- NULL
  
  nFeatures <- getNumberOfSignificantGenes( # all features have fdr <= 1
    input = tblList, fdr = 1, logFCrange = logFCrange) 
  nSignGenes <- getNumberOfSignificantGenes(
    input = tblList, fdr = fdr, logFCrange = logFCrange)
  nUpGenes <- getNumberOfRegulatedGenes(
    input = tblList, fdr = fdr, dir = "up", logFCrange = logFCrange)
  nDownGenes <- getNumberOfRegulatedGenes(
    input = tblList, fdr = fdr, dir = "down", logFCrange = logFCrange)
  nSignGenesPerc <- paste0(round(nSignGenes/nFeatures * 100, digits = 2), "%")
  nUpGenesPerc <- paste0(round(nUpGenes/nFeatures * 100, digits = 2), "%")
  nDownGenesPerc <- paste0(round(nDownGenes/nFeatures * 100, digits = 2), "%")
  coefInfoTbl <- cbind(nSignGenes, nUpGenes, nDownGenes, nSignGenesPerc)
  # create table with coef and coef label
  coefLabelTbl <- unique(tbl[, c("coef", grep(
      "^comparison([[:digit:]]{1,})?$", colnames(tbl), value = TRUE)
    ), drop = FALSE])
  # duplicate the table for direction ('up' and 'down')
  coefLabelTbl <- rbind(
    coefLabelTbl, coefLabelTbl[rep(seq_len(length(coef)), 1), ])
  # create table for barplot
  tbl <- data.frame(
    value = c(nUpGenes, nDownGenes),
    valueLabel = if(addPercentage) c(
      paste0(nUpGenes, " (", nUpGenesPerc, ")"), 
      paste0(nDownGenes, " (", nDownGenesPerc, ")")
    ) else c(nUpGenes, nDownGenes),
    direction = rep(c("up", "down"), each = nrow(coefInfoTbl)),
    total = rep(nSignGenes, 2),
    totalLabel = if(addPercentage) rep(
      paste0(nSignGenes, " (", nSignGenesPerc, ")"), 2
    ) else rep(nSignGenes, 2))
  tbl <- cbind(tbl, coefLabelTbl)
  tbl[, "direction"] <- factor(tbl[, "direction"], levels = c("down", "up"))
  
  tbl
}

#' Create ggplot object with barplot
#' @inheritParams daSignificantGenesBarplot
#' @param tbl combined summaries for all coefficients 
#' (i.e.: number of significantly up- and down-regulated genes)
#' @return ggplot object
#' @author Katarzyna Gorczak
#' @importFrom utils getFromNamespace
#' @importFrom ggh4x guide_axis_nested
#' @import ggplot2
callBarplot <- function(
    tbl, color, annotCex, addPercentage, title, titleCex, xlab, ylab, 
    axesTitleCex, axesCex, legendPosition, legendTitleCex, legendCex
) {
  
  labelVars <- grep("^comparison([[:digit:]]{1,})?$", colnames(tbl), value=TRUE)
  mainArgs <- c(
    list(x = paste0('interaction(', 
                    paste0(labelVars, collapse = ', ' ),', sep = "!")'), 
         y = 'value', label = 'valueLabel', fill = 'direction')
  )
  aesString <- c(    
    list(data = tbl, mapping = do.call(ggplot2::aes_string, mainArgs))
  )
  g <- do.call(getFromNamespace("ggplot", ns = "ggplot2"), aesString)
  
  tblTextTotal <- tbl[!duplicated(
    tbl[, c("total", "totalLabel", "coef", 
            grep("^comparison([[:digit:]]{1,})?$", colnames(tbl), value = TRUE)
  ), ]), ]
  
  g <- g + 
    geom_bar(stat="identity") +
    geom_text(size = annotCex, position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(values = color) + 
    geom_text(
      data = tblTextTotal,
      aes_string(x = paste0('interaction(', 
                            paste0(labelVars, collapse = ', ' ),', sep = "!")'),
                 y = "total",
                 label = "paste0('total: ', totalLabel)"),
      vjust = -0.2,
      size = annotCex) +
    scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = "!")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  g <- ggPlotTheme(
    ggObject = g, title = title, titleSize = titleCex, panelBackground = NULL,
    xTitle = xlab, xTitleSize = axesTitleCex, xTextSize = axesCex,
    yTitle = ylab, yTitleSize = axesTitleCex, yTextSize = axesCex,
    panelBorder = "grey", legendPosition = legendPosition, 
    legendTitleSize = legendTitleCex, legendTextSize = legendCex)
  
  g
}
