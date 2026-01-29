context("MA plot")

library(ggplot2)

testthat::test_that("The points correctly visualize the logFC and average expression for a mixed input of limma and top table", {

  # create the plot
  gg <- daMAplot(
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
        ggDataCoef <- subset(ggDataPoint, PANEL == 1+!!i)
        # sort by multiple columns, because the same logFC for multiple AveExpr
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
        list(
          logfc = ggDataCoef$y,
          mean = ggDataCoef$x
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "none")
        # sort by multiple columns, because the same logFC for multiple AveExpr
        topTableCoef <- topTableCoef[order(topTableCoef$logFC, topTableCoef$AveExpr, decreasing = TRUE), ]
        list(
          logfc = topTableCoef$logFC,
          mean = topTableCoef$AveExpr
        )
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataPoint, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
      list(
        logfc = ggDataCoef$y,
        mean = ggDataCoef$x
      )
    },
    expected = {
      topTableCoef <- inputMixed$A[order(inputMixed$A$logFC, inputMixed$A$AveExpr, decreasing = TRUE), ]
      list(
        logfc = topTableCoef$logFC,
        mean = topTableCoef$AveExpr
      )
    }
  )
  
})

testthat::test_that("The top genes correctly labelled", {
  
  # create the plot
  gg <- daMAplot(
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
        ggDataCoef <- subset(ggDataText, PANEL == 1+!!i)
        # sort by multiple columns, because the same logFC for multiple AveExpr
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
        ggDataCoef$label
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "none")
        topTableTopGenes <- topTableCoef[rank(topTableCoef[, "P.Value"], ties.method = "random") %in% 
                                           seq_len(10) | (nrow(topTableCoef) - rank(abs(topTableCoef[, "logFC"]), 
                                                                                    ties.method = "random") + 1) %in% seq_len(10), ]
        # sort by multiple columns, because the same logFC for multiple AveExpr
        topTableTopGenes <- topTableTopGenes[order(topTableTopGenes$logFC, topTableTopGenes$AveExpr, decreasing = TRUE), ]
        topTableTopGenes[, "ENTREZID"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataText, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
      ggDataCoef$label
    },
    expected = {
      topTableCoef <- inputMixed$A
      topTableCoef <- topTableCoef[rank(topTableCoef[, "P.Value"], ties.method = "random") %in% 
                                                         seq_len(10) | (nrow(topTableCoef) - rank(abs(topTableCoef[, "logFC"]), 
                                                                                                  ties.method = "random") + 1) %in% seq_len(10), ]
      topTableTopGenes <- topTableCoef[order(topTableCoef$logFC, topTableCoef$AveExpr, decreasing = TRUE), ]
      topTableTopGenes[, "ENTREZID"]
    }
  )
  
})

testthat::test_that("The genes of interest correctly labelled", {
  
  # create the plot
  gg <- daMAplot(
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
        ggDataCoef <- subset(ggDataText, PANEL == 1+!!i)
        # sort by multiple columns, because the same logFC for multiple AveExpr
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
        ggDataCoef$label
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "none")
        # sort by multiple columns, because the same logFC for multiple AveExpr
        topTableCoef <- topTableCoef[order(topTableCoef$logFC, topTableCoef$AveExpr, decreasing = TRUE), ]
        topTableCoef[which(topTableCoef$ENTREZID %in% features), "ENTREZID"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataText, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, ggDataCoef$x, decreasing = TRUE), ]
      ggDataCoef$label
    },
    expected = {
      topTableCoef <- inputMixed$A
      topTableCoef <- topTableCoef[order(topTableCoef$logFC, topTableCoef$AveExpr, decreasing = TRUE), ]
      topTableCoef[which(topTableCoef$ENTREZID %in% features), "ENTREZID"]
    }
  )
  
})
