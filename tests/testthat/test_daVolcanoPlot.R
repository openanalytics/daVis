context("Volcano plot")

library(ggplot2)

testthat::test_that("The points correctly visualize the logFC and p-value for a mixed input of limma and top table", {

  # create the plot
  gg <- daVolcanoPlot(
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
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
        list(
          pval = ggDataCoef$y,
          logfc = ggDataCoef$x
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "p")
        minVal <- .Machine$double.xmin
        topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
        topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
        list(
          pval = topTableCoef$P.ValuePlot,
          logfc = topTableCoef$logFC
        )
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataPoint, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
      list(
        pval = ggDataCoef$y,
        logfc = ggDataCoef$x
      )
    },
    expected = {
      topTableCoef <- inputMixed$A[order(inputMixed$A$P.Value, decreasing = FALSE), ]
      minVal <- .Machine$double.xmin
      topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
      topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
      list(
        pval = topTableCoef$P.ValuePlot,
        logfc = topTableCoef$logFC
      )
    }
  )
  
})

testthat::test_that("The top genes correctly labelled", {
  
  # create the plot
  gg <- daVolcanoPlot(
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
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
        ggDataCoef$label
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "p")
        minVal <- .Machine$double.xmin
        topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
        topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
        topTableTopGenes <- topTableCoef[rank(topTableCoef[, "P.Value"], ties.method = "random") %in% 
            seq_len(10) | (nrow(topTableCoef) - rank(abs(topTableCoef[, "logFC"]), 
                                          ties.method = "random") + 1) %in% seq_len(10), ]
        topTableTopGenes[, "ENTREZID"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataText, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
      ggDataCoef$label
    },
    expected = {
      topTableCoef <- inputMixed$A[order(inputMixed$A$P.Value, decreasing = FALSE), ]
      minVal <- .Machine$double.xmin
      topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
      topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
      topTableTopGenes <- topTableCoef[rank(topTableCoef[, "P.Value"], ties.method = "random") %in% 
                                         seq_len(10) | (nrow(topTableCoef) - rank(abs(topTableCoef[, "logFC"]), 
                                                                                  ties.method = "random") + 1) %in% seq_len(10), ]
      topTableTopGenes[, "ENTREZID"]
    }
  )
  
})

testthat::test_that("The genes of interest correctly labelled", {
  
  # create the plot
  gg <- daVolcanoPlot(
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
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
        ggDataCoef$label
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "p")
        minVal <- .Machine$double.xmin
        topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
        topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
        topTableCoef[which(topTableCoef$ENTREZID %in% features), "ENTREZID"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataText, PANEL == 1)
      ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
      ggDataCoef$label
    },
    expected = {
      topTableCoef <- inputMixed$A[order(inputMixed$A$P.Value, decreasing = FALSE), ]
      minVal <- .Machine$double.xmin
      topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
      topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
      topTableCoef[which(topTableCoef$ENTREZID %in% features), "ENTREZID"]
    }
  )
  
})

testthat::test_that("The coefficient labels are specified as a list of functions", {
  
  # coefficients by 
  coefs <- c("L.LvsP", "L.PvsV", "B.LvsP", "B.PvsV")
  coefLabel <- list(
    function(x) sub("(.+)\\.(.+)", "\\1", x),
    function(x) sub("(.+)\\.(.+)", "\\2", x)
  )
  
  # create the plot
  gg <- daVolcanoPlot(
    input = res.limma,
    coef = coefs, coefLabel = coefLabel
  )
  
  ## check that data is sorted according to the coefficients
  
  # extract data behind the points
  idxGeomPoint <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
  ggDataPoint <- ggplot2::layer_data(plot = gg, i = idxGeomPoint)
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefs)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataPoint, PANEL == !!i)
        ggDataCoef <- ggDataCoef[order(ggDataCoef$y, decreasing = TRUE), ]
        list(
          pval = ggDataCoef$y,
          logfc = ggDataCoef$x
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefs[!!i], number = Inf, sort.by = "p")
        minVal <- .Machine$double.xmin
        topTableCoef$P.ValuePlot <- ifelse(topTableCoef[, "P.Value"] < minVal, minVal, topTableCoef[, "P.Value"])
        topTableCoef$P.ValuePlot <- -log10(topTableCoef$P.ValuePlot)
        list(
          pval = topTableCoef$P.ValuePlot,
          logfc = topTableCoef$logFC
        )
      }
    )
  }
  
  ## check that the correct facet labels are displayed
  
  # from plot (to improve)
  ggGrob <- ggplot2::ggplotGrob(gg)
  ggGrobFacets <- ggGrob$grobs[grep("^strip", ggGrob$layout$name)]
  ggFacetLabels <- sapply(ggGrobFacets, function(grob){
    child <- grob$grobs[[1]]$children
    titles <- child[sapply(child, inherits, "titleGrob")]
    titles[[1]]$children[[1]]$label
  })
  
  # from data
  dataFacetLabels <- lapply(coefLabel, function(fct) fct(coefs))
  # keep unique labels
  dataFacetLabels  <- unlist(lapply(dataFacetLabels , function(l) rle(l)$values))
  
  expect_identical(object = ggFacetLabels, expected = dataFacetLabels)
  
})
