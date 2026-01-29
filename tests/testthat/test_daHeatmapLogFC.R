context("Heatmap of logFC")

library(ggplot2)

testthat::test_that("Heatmap correctly visualizes the logFC for a mixed input of limma and top table, and features are correctly ordered", {

  # create the plot
  gg <- daHeatmapLogFC(
    input = inputMixed,
    coef = coefMixed,
    features = features
  ) + ggplot2::geom_text(aes_string(label = "logFC"))
  
  # extract data behind the tiles
  idxGeomTile <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomText"), logical(1)))
  ggDataTiles <- ggplot2::layer_data(plot = gg, i = idxGeomTile)
  
  # extract top tables
  topTableModel <- limma::topTable(fit = res.limma, coef = coefsModel, number = Inf, sort.by = "none")
  topTableModel <- topTableModel[match(features, topTableModel$ENTREZID), ]
  topTable <- inputMixed$A[match(features, inputMixed$A$ENTREZID), ]
  logfcTbl <- cbind(topTable$logFC, topTableModel[, coefsModel])
  
  logfc <- as.matrix(logfcTbl)
  rowm <- rowMeans(logfc, na.rm = FALSE)
  requireNamespace("stats")
  diste <- stats::dist(logfc) # distance matrix method = euclidean
  hcc <- stats::hclust(diste) # hierarchical clustering method = complete
  dend <- stats::as.dendrogram(hcc) 
  dend <- stats::reorder(dend, rowm)
  rowInd <- stats::order.dendrogram(dend)
  
  logfcTbl <- logfcTbl[rev(rowInd), ]
  
  # check that logFC correctly displayed for the mixed input
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataTiles, x == 1+!!i)
        ggDataCoef$label
      },
      expected = {
        logfcTbl[, 1+!!i, drop = TRUE]
      }
    )
  }
  
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataTiles, x == 1)
      ggDataCoef$label
    },
    expected = {
      logfcTbl[, 1, drop = TRUE]
    }
  )
  
  # check that the names of the features correctly displayed
  ggBuild <- ggplot2::ggplot_build(plot = gg)
  expect_identical(
    object = levels(ggBuild$plot$data$featuresVar),
    expected = rownames(logfcTbl)
  )
  
})

testthat::test_that("Heatmap correctly visualizes the logFC for a mixed input with different set of features", {
  
  # create the plot
  expect_warning(
    gg <- daHeatmapLogFC(
      input = inputMixedRep,
      coef = coefMixed,
      features = features, 
    ) + ggplot2::geom_text(aes_string(label = "logFC")),
    "66292, 68585) are not present for the coefficient A"
  )
  
  # extract data behind the tiles
  idxGeomTile <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomText"), logical(1)))
  ggDataTiles <- ggplot2::layer_data(plot = gg, i = idxGeomTile)

  # extract top tables
  topTableModel <- limma::topTable(fit = res.limmaRep, coef = coefsModel, number = Inf, sort.by = "none")
  topTableModel <- topTableModel[match(features, topTableModel$ENTREZID), ]
  topTable <- inputMixedRep$A[match(features, inputMixedRep$A$ENTREZID), ]
  rownames(topTable) <- features
  topTable$ENTREZID <- features
  logfcTbl <- cbind(topTable$logFC, topTableModel[, coefsModel])

  logfc <- as.matrix(logfcTbl)
  rowm <- rowMeans(logfc, na.rm = FALSE)
  requireNamespace("stats")
  diste <- stats::dist(logfc) # distance matrix method = euclidean
  hcc <- stats::hclust(diste) # hierarchical clustering method = complete
  dend <- stats::as.dendrogram(hcc)
  dend <- stats::reorder(dend, rowm)
  rowInd <- stats::order.dendrogram(dend)

  logfcTbl <- logfcTbl[rev(rowInd), ]

  # check that logFC correctly displayed for the mixed input
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataTiles, x == 1+!!i)
        ggDataCoef$label
      },
      expected = {
        logfcTbl[, 1+!!i, drop = TRUE]
      }
    )
  }
  
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataTiles, x == 1)
      ggDataCoef$label
    },
    expected = {
      logfcTbl[, 1, drop = TRUE]
    }
  )

  # check that the names of the features correctly displayed
  ggBuild <- ggplot2::ggplot_build(plot = gg)
  expect_identical(
    object = levels(ggBuild$plot$data$featuresVar),
    expected = rownames(logfcTbl)
  )
  
})

testthat::test_that("The feature colors, when not named, are correctly set for multiple features", {
  
  featuresColor <- rep(c("black", "grey"), length.out = length(features))
  
  gg <- daHeatmapLogFC(
    input = res.limma,
    features = features,
    featuresIdVar = "ENTREZID",
    coef = coefsModel,
    featuresColor = featuresColor
  )
    
    # extract labels of the y-axis (from bottom to top)
  idxGeomTile <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomTile"), logical(1)))
  yAxisLabel <- ggplot2::layer_scales(gg, idxGeomTile)$y$range$range
  
  # extract colors (from bottom to top)
  ggBuild <- ggplot2::ggplot_build(gg)
  yAxisColor <- ggBuild$plot$theme$axis.text.y$colour
  
  expect_identical(
    object = setNames(yAxisColor, yAxisLabel)[features],
    expected = setNames(featuresColor, features)
  )
    
})

