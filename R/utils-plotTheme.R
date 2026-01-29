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
  
  t <- ggplot2::theme_bw()
  if (length(xTextSize) == 2) xTextSize <- xTextSize[1]
  if (length(yTextSize) == 2) yTextSize <- yTextSize[2]
  
  if (!is.null(panelBackground)) {
    t <- t + ggplot2::theme(panel.background = ggplot2::element_rect(
      fill = panelBackground, colour = NA))
  } else t <- t + ggplot2::theme(panel.background = ggplot2::element_blank())

  t <- t + ggplot2::theme(
    panel.border = ggplot2::element_rect(
      linetype = "solid", fill = NA, color = panelBorder),
    legend.title = ggplot2::element_text(size = ggplot2::rel(legendTitleSize)), 
    legend.text = ggplot2::element_text(size = ggplot2::rel(legendTextSize)), 
    legend.position = legendPosition,
    axis.text.y = ggtext::element_markdown(
      hjust = 0, color = yTextColor, size = ggplot2::rel(yTextSize)),
    axis.text.x = ggtext::element_markdown(
      hjust = xTextHjust, angle = xTextAngle, color = xTextColor, 
      size = ggplot2::rel(xTextSize)),
    strip.text = ggplot2::element_text(
      size = ggplot2::rel(facetLabelSize), color = facetLabelColor)
  )
  if (!xTicks) t <- t + ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
  if (!yTicks) t <- t + ggplot2::theme(axis.ticks.y = ggplot2::element_blank())
  if (!gridMajor) t <- t + ggplot2::theme(
    panel.grid.major = ggplot2::element_blank())
  if (!gridMinor) t <- t + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank())
  
  if (!is.null(xTitle)) {
    ggObject <- ggObject + ggplot2::xlab(xTitle)
    t <- t + ggplot2::theme(axis.title.x = ggplot2::element_text(
      size = ggplot2::rel(xTitleSize), vjust = -0.5))
  } else t <- t + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  if (!is.null(yTitle)) {
    ggObject <- ggObject + ggplot2::ylab(yTitle)
    t <- t + ggplot2::theme(axis.title.y = ggplot2::element_text(
      size = ggplot2::rel(yTitleSize), vjust = 2))
  } else t <- t + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  
  if (!is.null(title)) {
    ggObject <- ggObject + ggplot2::labs(title = title)
    t <- t +
      ggplot2::theme(plot.title = ggplot2::element_text(
        size = ggplot2::rel(titleSize)))
  }
  
  ggObject + t
}