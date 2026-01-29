context("Waterfall plot")

library(ggplot2)

testthat::test_that("The waterfall bars correctly visualize the logFC for a mixed input of limma and top table", {

  # create the plot
  gg <- daWaterfallPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features
  )
  
  # extract data behind the bar
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  ggDataBar <- ggplot2::layer_data(plot = gg, i = idxGeomBar)
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataBar, PANEL == 1+!!i)
        ggDataCoef <- ggDataCoef$xmin + ggDataCoef$xmax
        ggDataCoef[order(ggDataCoef)]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf)
        topTableCoef <- topTableCoef[match(features, topTableCoef$ENTREZID), "logFC"]
        topTableCoef[order(topTableCoef)]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataBar, PANEL == 1)
      ggDataCoef <- ggDataCoef$xmin + ggDataCoef$xmax
      ggDataCoef[order(ggDataCoef)]
    },
    expected = {
      topTableCoef <- inputMixed$A[match(features, inputMixed$A$ENTREZID), "logFC"]
      topTableCoef[order(topTableCoef)]
    }
  )
  
})
