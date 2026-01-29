context("Scatter plot")

library(ggplot2)

testthat::test_that("The points correctly visualize the logFC and p-value for a mixed input of limma and top table", {

  # create the plot
  gg <- daScatterPlot(
    input = inputMixed,
    coef = coefMixed
  )
  
  # extract data behind the points
  idxGeomPoint <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
  ggDataPoint <- ggplot2::layer_data(plot = gg, i = idxGeomPoint)
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataPoint, PANEL == !!i)
        ggDataCoef$y[order(ggDataCoef$y, decreasing = TRUE)]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "none")
        topTableCoef$logFC[order(topTableCoef$logFC, decreasing = TRUE)]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataPoint, PANEL == 1)
      ggDataCoef$x[order(ggDataCoef$x, decreasing = TRUE)]
    },
    expected = {
      topTableCoef <- inputMixed$A
      topTableCoef$logFC[order(topTableCoef$logFC, decreasing = TRUE)]
    }
  )
  
})

testthat::test_that("The top genes correctly labelled", {
  
  # create the plot
  gg <- daScatterPlot(
    input = inputMixed,
    coef = coefMixed,
    topGenes = 10
  )
  
  # extract data behind the points
  idxGeomText <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomTextRepel"), logical(1)))
  ggDataText <- ggplot2::layer_data(plot = gg, i = idxGeomText)
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataText, PANEL == !!i)
        ggDataCoef$label[order(ggDataCoef$y, decreasing = TRUE)]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "p")
        topTableTopGenes <- topTableCoef[rank(topTableCoef[, "P.Value"], ties.method = "random") %in% 
            seq_len(10) | (nrow(topTableCoef) - rank(abs(topTableCoef[, "logFC"]), 
                                          ties.method = "random") + 1) %in% seq_len(10), ]
        topTableTopGenes <- topTableTopGenes[order(topTableTopGenes$logFC, decreasing = TRUE), ]
        topTableTopGenes[, "ENTREZID"]
      }
    )
  }
  
})

testthat::test_that("The genes of interest correctly labelled", {
  
  # create the plot
  gg <- daScatterPlot(
    input = inputMixed,
    coef = coefMixed,
    topGenes = 10,
    genesToHighlight = features
  )
  
  # extract data behind the points 
  idxGeomText <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomTextRepel"), logical(1)))
  ggDataText <- ggplot2::layer_data(plot = gg, i = idxGeomText[1]) # first, because the second is for top genes
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataText, PANEL == !!i)
        ggDataCoef$label[order(ggDataCoef$y, decreasing = TRUE)]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "p")
        topTableCoef <- topTableCoef[order(topTableCoef$logFC, decreasing = TRUE), ]
        topTableCoef[which(topTableCoef$ENTREZID %in% features), "ENTREZID"]
      }
    )
  }
  
})
