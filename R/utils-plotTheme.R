#' Customized ggplot theme
#' 
#' @param ggObject ggplot object
#' @param panelBackground panel background of the entire plot
#' @param panelBorder border around plotting area, by default black
#' @param title plot title
#' @param titleSize cex for plot title, relative to the parent
#' @param xTextSize cex for axis text
#' @param xTextColor color font on the x axis, by default black
#' @param xTitle title of the x axis
#' @param xTitleSize cex for the x axis title
#' @param xTextAngle angle of the x axis labels
#' @param xTextHjust numeric for horizontal adjustment, 0.5 by default
#' @param yTextSize cex for axis text
#' @param yTextColor color font on the y axis, by default black 
#' @param yTitle title of the y axis, NULL by default  
#' @param yTitleSize cex for the y axis title
#' @param facetLabelSize cex for the facet labels
#' @param facetLabelColor color of the text in facets
#' @param gridMajor logical whether to show the major grid lines
#' @param gridMinor logical whether to show the minor grid lines
#' @param yTicks logical whether to show ticks on the y axis 
#' @param xTicks logical whether to show ticks on the x axis 
#' @param legendPosition legend position, by default right
#' @param legendTitleSize cex for the legend title
#' @param legendTextSize cex for the legend 
#' @importFrom ggplot2 theme_bw theme element_rect element_blank element_text 
#' @importFrom ggplot2 rel xlab ylab labs
#' @return ggplot object with customized theme
ggPlotTheme <- function(
    ggObject,
    panelBackground = NULL, panelBorder = "black",
    title = NULL, titleSize = 2,
    xTextSize = 1.2, xTextColor = "black", xTitle = NULL, 
    xTitleSize = 1.4, xTextAngle = 0, xTextHjust = 0.5,
    yTextSize = 1.2, yTextColor = "black", yTitle = NULL, yTitleSize = 1.4,
    facetLabelSize = 1.5, facetLabelColor = "black",
    gridMajor = TRUE, gridMinor = TRUE, yTicks = TRUE, xTicks = TRUE,
    legendPosition = "right", legendTitleSize = 1.3, legendTextSize = 1.1
) {
  
  t <- theme_bw()
  if (length(xTextSize) == 2) xTextSize <- xTextSize[1]
  if (length(yTextSize) == 2) yTextSize <- yTextSize[2]
  
  if (!is.null(panelBackground)) {
    t <- t + theme(panel.background = element_rect(
      fill = panelBackground, colour = NA))
  } else t <- t + theme(panel.background = element_blank())
  
  t <- t + theme(
    panel.border = element_rect(
      linetype = "solid", fill = NA, color = panelBorder),
    legend.title = element_text(size = rel(legendTitleSize)), 
    legend.text = element_text(size = rel(legendTextSize)), 
    legend.position = legendPosition,
    axis.text.y = ggtext::element_markdown(
      hjust = 0, color = yTextColor, size = rel(yTextSize)),
    axis.text.x = ggtext::element_markdown(
      hjust = xTextHjust, angle = xTextAngle, color = xTextColor, 
      size = rel(xTextSize)),
    strip.text = element_text(
      size = rel(facetLabelSize), color = facetLabelColor)
  )
  if (!xTicks) t <- t + theme(axis.ticks.x = element_blank())
  if (!yTicks) t <- t + theme(axis.ticks.y = element_blank())
  if (!gridMajor) t <- t + theme(panel.grid.major = element_blank())
  if (!gridMinor) t <- t + theme(panel.grid.minor = element_blank())
  
  if (!is.null(xTitle)) {
    ggObject <- ggObject + xlab(xTitle)
    t <- t + theme(axis.title.x = element_text(
      size = rel(xTitleSize), vjust = -0.5))
  } else t <- t + theme(axis.title.x = element_blank())
  
  if (!is.null(yTitle)) {
    ggObject <- ggObject + ylab(yTitle)
    t <- t + theme(axis.title.y = element_text(
      size = rel(yTitleSize), vjust = 2))
  } else t <- t + theme(axis.title.y = element_blank())
  
  if (!is.null(title)) {
    ggObject <- ggObject + labs(title = title)
    t <- t +
      theme(plot.title = element_text(size = rel(titleSize)))
  }
  
  ggObject + t
}