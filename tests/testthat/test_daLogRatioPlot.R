context("Plot log ratios")

library(ggplot2)

testthat::test_that("The bars correctly visualize the logFC for a mixed input of limma and top table", {

  # create the plot
  gg <- daLogRatioPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features,
    errorBars = TRUE
  )
  
  # extract data behind the bar
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  ggDataBar <- ggplot2::layer_data(plot = gg, i = idxGeomBar)
  
  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataBar, PANEL == 1+!!i)
        ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "y"]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf)
        topTableCoef[match(features, topTableCoef$ENTREZID), "logFC"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataBar, PANEL == 1)
      ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "y"]
    },
    expected = {
      inputMixed$A[match(features, inputMixed$A$ENTREZID), "logFC"]
    }
  )
  
})

testthat::test_that("The bars correctly visualize the logFC with SE", {
  
  # create the plot
  gg <- suppressWarnings(daLogRatioPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features,
    errorBars = TRUE
  ))
  
  # extract data behind the error bar
  idxGeomErrorbar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomErrorbar"), logical(1)))
  ggDataErrorbar <- ggplot2::layer_data(plot = gg, i = idxGeomErrorbar)
  
  # check that logFC+SE correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataErrorbar, PANEL == 1+i)
        list(
          ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymin"],
          ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymax"]
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[i], number = Inf, sort.by = "none")
        seDataCoef <- c(res.limma$stdev.unscaled[, coefsModel[i]] * sqrt(res.limma$s2.post))
        topTableCoef$min <- topTableCoef[, "logFC"] - seDataCoef
        topTableCoef$max <- topTableCoef[, "logFC"] + seDataCoef
        list(
          topTableCoef[match(features, topTableCoef$ENTREZID), "min"],
          topTableCoef[match(features, topTableCoef$ENTREZID), "max"]
        )
      }
    )
  }
  
  # check that logFC+SE correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataErrorbar, PANEL == 1)
      list(
        ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymin"],
        ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymax"]
      )
    },
    expected = {
      topTable <- inputMixed[["A"]]
      topTable <- topTable[match(features, topTable$ENTREZID), ]
      list(
        with(topTable, logFC - se),
        with(topTable, logFC + se)
      )
    }
  )
  
})

testthat::test_that("The bars correctly visualize the logFC with SE for mixed input when SE column not present by default", {
  
  inputMixed$A <- inputMixed$A[, -grep("se", colnames(inputMixed$A), ignore.case = TRUE)]
    
  # create the plot
  gg <- suppressWarnings(daLogRatioPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features,
    errorBars = TRUE
  ))
  
  # extract data behind the error bar
  idxGeomErrorbar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomErrorbar"), logical(1)))
  ggDataErrorbar <- ggplot2::layer_data(plot = gg, i = idxGeomErrorbar)
  
  # check that logFC+SE correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataErrorbar, PANEL == 1+i)
        list(
          ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymin"],
          ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymax"]
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[i], number = Inf, sort.by = "none")
        seDataCoef <- c(res.limma$stdev.unscaled[, coefsModel[i]] * sqrt(res.limma$s2.post))
        topTableCoef$min <- topTableCoef[, "logFC"] - seDataCoef
        topTableCoef$max <- topTableCoef[, "logFC"] + seDataCoef
        list(
          topTableCoef[match(features, topTableCoef$ENTREZID), "min"],
          topTableCoef[match(features, topTableCoef$ENTREZID), "max"]
        )
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataErrorbar, PANEL == 1)
      all(is.na(c(ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymin"],
                  ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "ymax"])))
    },
    expected = {
      TRUE
    }
  )
  
})

testthat::test_that("The features are ordered by similarity", {
  
  # create the plot
  gg <- daLogRatioPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features,
    featuresOrder = "similarity",
    errorBars = TRUE
  )
  
  # extract data with the statistics
  topTableModel <- c(
    sapply(coefsModel, limma::topTable, fit = res.limma, number = Inf, simplify = FALSE, sort.by = "none"),
    list(A = inputMixed$A)
  )
  topTableStat <- sapply(topTableModel, function(topTable) 
    topTable[match(features, topTable$ENTREZID), "t"])
  featuresOrdered <- features[stats::hclust(stats::dist(topTableStat))$order]
  
  
  # TODO: update to extract the labels of the x-axis
  ggBuild <- ggplot2::ggplot_build(plot = gg)
  expect_identical(
    object = levels(ggBuild$plot$data$featuresVar),
    expected = rev(featuresOrdered)
  )
  
})

testthat::test_that("The features are ordered by significance", {
  
  # create the plot
  gg <- daLogRatioPlot(
    input = inputMixed,
    coef = coefMixed,
    features = features,
    featuresOrder = "significance",
    errorBars = TRUE
  )
  
  # extract data with the statistics
  topTableModel <- c(
    sapply(coefsModel, limma::topTable, fit = res.limma, number = Inf, simplify = FALSE),
    list(A = inputMixed$A)
  )
  topTablePValue <- sapply(topTableModel, function(topTable) 
    topTable[match(features, topTable$ENTREZID), "P.Value"])
  featuresOrdered <- features[order(rowMeans(-log10(topTablePValue)), decreasing = TRUE)]
  
  # TODO: update to extract the labels of the x-axis
  ggBuild <- ggplot2::ggplot_build(plot = gg)
  expect_identical(
    object = levels(ggBuild$plot$data$featuresVar),
    expected = rev(featuresOrdered)
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
  gg <- daLogRatioPlot(
    input = res.limma,
    coef = coefs, coefLabel = coefLabel,
    features = features,
    errorBars = TRUE
  )
  
  ## check that data is sorted according to the coefficients
  
  # extract data behind the bar
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  ggDataBar <- ggplot2::layer_data(plot = gg, i = idxGeomBar)
  
  # check that correct logFC displayed in each panel
  for(i in seq_along(coefs)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataBar, PANEL == !!i)
        ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "y"]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefs[!!i], number = Inf)
        topTableCoef[match(features, topTableCoef$ENTREZID), "logFC"]
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

testthat::test_that("An error is generated if the coefficient labels are duplicated", {
  
  # coefficients
  coefs <- c("L.LvsP", "L.PvsV")
  coefLabel <- list(
    c(`L.LvsP` = "L", L.PvsV = "L"),
    c(`L.LvsP` = "L1", L.PvsV = "L1")
  )
  
  # create the plot
  expect_error(
    daLogRatioPlot(
      input = res.limma,
      coef = coefs, coefLabel = coefLabel,
      features = features,
      errorBars = TRUE
    ),
    "labels are not unique"
  )
  
})

testthat::test_that("Log ratio plot correctly visualizes the logFC for a mixed input with different set of features", {
  
  # create the plot
  expect_warning(
    gg <- daLogRatioPlot(
      input = inputMixedRep,
      coef = coefMixed,
      features = features,
      errorBars = TRUE
    ),
    "66292, 68585) are not present for the coefficient A"
  )

  # extract data behind the bar
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  ggDataBar <- ggplot2::layer_data(plot = gg, i = idxGeomBar)

  # check that logFC correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataBar, PANEL == 1+!!i)
        ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "y"]
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limmaRep, coef = coefsModel[!!i], number = Inf)
        topTableCoef[match(features, topTableCoef$ENTREZID), "logFC"]
      }
    )
  }
  
  # check that logFC correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataBar, PANEL == 1)
      ggDataCoef[order(ggDataCoef$x, decreasing = TRUE), "y"]
    },
    expected = {
      inputMixedRep$A[match(features, inputMixedRep$A$ENTREZID), "logFC"]
    }
  )
  
})

testthat::test_that("The feature colors, when not named, are correctly set for multiple features", {
  
  featuresColor <- rep(c("black", "grey"), length.out = length(features))
  
  for(featuresOrder in list(NULL, "similarity", "significance")){

    gg <- daLogRatioPlot(
      input = res.limma,
      coef = "L.LvsP",
      features = features, 
      featuresOrder = featuresOrder,
      featuresColor = featuresColor
    )
    
    # extract labels of the y-axis (from bottom to top)
    idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
    yAxisLabel <- ggplot2::layer_scales(gg, idxGeomBar)$x$range$range

    # extract colors (from bottom to top)
    ggBuild <- ggplot2::ggplot_build(gg)
    yAxisColor <- ggBuild$plot$theme$axis.text.y$colour
    
    expect_identical(
      object = setNames(yAxisColor, yAxisLabel)[features],
      expected = setNames(featuresColor, features),
      info = paste0("Feature colors are not properly set if featuresOrder is: ",
        shQuote(featuresOrder))
    )
    
  }
  
})

testthat::test_that("A significance star is correctly displayed in a log ratio plot", {
  
  getSignifStar <- function(topTable){
    with(topTable, as.character(
      stats::symnum(
        x = adj.P.Val, 
        cutpoints = c(0, .001, .01, .05, .1, 1),
        symbols = c("***","**","*","."," "),
        corr = FALSE
      )
    ))
  }
  
  gg <- daLogRatioPlot(
    input = res.limma,
    coef = "L.LvsP",
    features = features,
    text = getSignifStar
  )
  
  # extract labels of the y-axis (from bottom to top)
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  yAxisLabel <- ggplot2::layer_scales(gg, idxGeomBar)$x$range$range
  
  # extract text
  idxGeomText <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomText"), logical(1)))
  ggDataStar <- ggplot2::layer_data(gg, idxGeomText)[, c("x", "label")]
  ggDataStar[, "featureID"] <- yAxisLabel 
  ggDataStar <- ggDataStar[match(features, ggDataStar$featureID), ]
  
  topTable <- limma::topTable(fit = res.limma, coef = "L.LvsP", n = Inf)
  topTable <- topTable[match(features, rownames(topTable)), ]
  topTable[, "star"] <- getSignifStar(topTable = topTable)
  
  expect_identical(
    object = ggDataStar[, "label"],
    expected = topTable[, "star"]
  )
  
})

testthat::test_that("The features are ordered by similarity if a text variable is specified", {
  
  coefs <- c("L.LvsP", "L.PvsV")
  input <- topTableList[coefs]
  
  # add extra column in the top table
  input <- sapply(input, function(x){
    x$text <- x$SYMBOL
    x
  }, simplify = FALSE)
  
  expect_silent(
    gg <- daLogRatioPlot(
      input = input,
      coef = coefs,
      features = features,
      featuresOrder = "similarity",
      text = "text"
    )
  )
  
  # extract labels of the y-axis (from bottom to top)
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  yAxisLabel <- ggplot2::layer_scales(gg, idxGeomBar)$x$range$range
  
  tStat <- sapply(coefs, function(coef){
    limma::topTable(fit = res.limma, n = Inf, coef = coef)[features, "t"]
  })
  featuresOrdered <- features[stats::hclust(stats::dist(tStat))$order]
  
  expect_identical(
    object = rev(yAxisLabel), # from top to the bottom
    expected = featuresOrdered
  )
  
})

testthat::test_that("The text variable is correctly displayed if error bars are not available for all coefficients", {
  
  coefs <- c("L.LvsP", "L.PvsV")
  input <- topTableList[coefs]
  
  # add extra column in the top table
  input <- sapply(input, function(x){
    x$text <- x$SYMBOL
    x
  }, simplify = FALSE)
  
  # remove standard errors for one coefficient
  input[[2]]$se <- NULL
  
  expect_warning(
    gg <- daLogRatioPlot(
      input = input,
      coef = coefs,
      features = features,
      featuresOrder = "similarity",
      text = "text"
    ),
    "Standard errors.+are not available"
  )
  
  # extract labels of the y-axis (from bottom to top)
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  yAxisLabel <- ggplot2::layer_scales(gg, idxGeomBar)$x$range$range
  
  # extract text
  idxGeomText <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomText"), logical(1)))
  ggText <- ggplot2::layer_data(gg, idxGeomText)
  ggText[, "featureID"] <- yAxisLabel[ggText$x]
  
  ## if error bars are available, text is positioned at logFC +- SE:
  
  # check if text is correctly extracted
  ggTextEB <- subset(ggText, PANEL == 1)
  ggTextEB <- ggTextEB[match(features, ggTextEB$featureID), ]
  inputEB <- input[[1]][features, ]
  expect_identical(
    object = ggTextEB[, "label"],
    expected = inputEB[, "text"]
  )   
  # check that the text position is available (not NA)
  expect_identical(
    object = ggTextEB[, "y"],
    expected = with(inputEB, ifelse(logFC < 0, logFC - se, logFC + se))
  )   
  
  ## if error bars are not available, text is positioned at logFC
  
  # check if text is correctly extracted
  ggTextNoEB <- subset(ggText, PANEL == 2)
  ggTextNoEB <- ggTextNoEB[match(features, ggTextNoEB$featureID), ]
  inputNoEB <- input[[2]][features, ]
  expect_identical(
    object = ggTextNoEB[, "label"],
    expected = inputNoEB[, "text"]
  )   
  # check that the text position is available (not NA)
  expect_identical(
    object = ggTextNoEB[, "y"],
    expected = with(inputNoEB, logFC)
  ) 
  
})

