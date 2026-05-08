context("Upset plot")

testthat::test_that("Data for upset plot with up-regulated genes is correct", {

  # create data for upset plot
  dt <- createDataUpsetPlot(
    input = res.limma, 
    coef = colnames(res.limma), 
    coefLabel = colnames(res.limma), 
    featuresIdVar = "SYMBOL", 
    fdr = 0.05, 
    dir = "up"
  )
  
  # extract genes based on fdr and logfc
  geneList <- lapply(colnames(res.limma), function(x) {
    tmp <- limma::topTable(res.limma, coef = x, number = Inf)
    tmp[tmp$adj.P.Val <= 0.05 & tmp$logFC > 0, "SYMBOL"]
  })
  genes <- Reduce(union, geneList)
  genes[is.na(genes)] <- paste0("NA.", seq_len(length(genes[is.na(genes)])))
  
  expect_identical(
    object = nrow(dt),
    expected = length(genes)
  )
  
  expect_identical(
    object = sort(rownames(dt)),
    expected = sort(genes)
  )
  
})

testthat::test_that("Data for upset plot with down-regulated genes is correct", {
  
  # create data for upset plot
  dt <- createDataUpsetPlot(
    input = res.limma, 
    coef = colnames(res.limma), 
    coefLabel = colnames(res.limma), 
    featuresIdVar = "SYMBOL", 
    fdr = 0.05, 
    dir = "down"
  )
  
  # extract genes based on fdr and logfc
  geneList <- lapply(colnames(res.limma), function(x) {
    tmp <- limma::topTable(res.limma, coef = x, number = Inf)
    tmp[tmp$adj.P.Val <= 0.05 & tmp$logFC < 0, "SYMBOL"]
  })
  genes <- Reduce(union, geneList)
  genes[is.na(genes)] <- paste0("NA.", seq_len(length(genes[is.na(genes)])))
  
  expect_identical(
    object = nrow(dt),
    expected = length(genes)
  )
  
  expect_identical(
    object = sort(rownames(dt)),
    expected = sort(genes)
  )
  
})

testthat::test_that("Overlapping sets are correctly created for up-regulated genes", {

  dt <- createDataUpsetPlot(
    input = res.limma, 
    coef = colnames(res.limma), 
    coefLabel = colnames(res.limma), 
    featuresIdVar = "SYMBOL", 
    fdr = 0.05, 
    dir = "up"
  )
  
  sets <- extractFeatures(dt)
  
  for(i in seq_along(names(sets))){
    expect_identical(
      object = length(which(rowSums(dt) == !!i)),
      expected = sum(sapply(sets[[!!i]], length))
    )
  }
  
})

testthat::test_that("Overlapping sets are correctly created for down-regulated genes", {
  
  dt <- createDataUpsetPlot(
    input = res.limma, 
    coef = colnames(res.limma), 
    coefLabel = colnames(res.limma), 
    featuresIdVar = "SYMBOL", 
    fdr = 0.05, 
    dir = "down"
  )
  
  sets <- extractFeatures(dt)
  
  for(i in seq_along(names(sets))){
    expect_identical(
      object = length(which(rowSums(dt) == !!i)),
      expected = sum(sapply(sets[[!!i]], length))
    )
  }
  
})

testthat::test_that("Test for upset plot", {
  
  ttList <- topTableList
  ttList$L.LvsP$logFC <- ttList$L.PvsV$logFC <- ttList$B.LvsP$logFC <- -2
  
  expect_error(
    daUpset(
      input = ttList, 
      coef = names(ttList)[1], 
      featuresIdVar = "SYMBOL", 
      dir = "up"
    ),
    "At least 2 'coef' must be provided."
  )
  
  expect_error(
    daUpset(
      input = ttList, 
      coef = names(ttList), 
      featuresIdVar = "SYMBOL", 
      dir = "up"
    ),
    "There are up-regulated significant features only for one"
  )
  
  ttList$B.PvsV$logFC <- -2
  
  expect_error(
    daUpset(
      input = ttList, 
      coef = names(ttList), 
      featuresIdVar = "SYMBOL", 
      dir = "up"
    ),
    "No features are up-regulated"
  )
  
})

testthat::test_that("Test for upset plot", {
  
  # use suppressWarnings because UpSetR::upset uses aes_string internally 
  # this will throw a warning
  g <- suppressWarnings(daUpset(
    input = topTableList, 
    coef = names(topTableList), 
    featuresIdVar = "SYMBOL", 
    returnAnalysis = TRUE,
    dir = "up"
  ))
  
  expect_identical(
    object = names(g),
    expected = c("sets", "plot")
  )
  
  expect_identical(
    object = names(g$sets),
    expected = c("Intersect 1 set", paste0("Intersect ", 2:4, " sets"))
  )
  
})


