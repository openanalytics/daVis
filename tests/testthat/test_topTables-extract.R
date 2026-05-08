context("Extract top tables")

testthat::test_that("The standard errors are correctly extracted for a limma model", {
  
  coef <- colnames(coefficients(res.limma))[1]
  se <- getSEModel(input = res.limma, coef = coef)
  
  expect_equal(
    object = se,
    expected = setNames(
      with(res.limma, sqrt(s2.post) * stdev.unscaled[, coef]),
      rownames(res.limma)
    )
  )
  
})


testthat::test_that("The standard errors are missing for a edgeR model", {
  
  coef <- colnames(coefficients(res.edger))[1]
  expect_warning(
    se <- getSEModel(input = res.edger, coef = coef)
  )
  
  expect_equal(object = se, expected = NA_real_)
  
})


testthat::test_that("The top table is extracted correctly for DESeq2", {
  
  tt <- extractTopTables(
    input = res.deseq, 
    coef = coefsModel[1], 
    columns = c("ENTREZID", "logFC", "P.Value", "adj.P.Val"),
    hoverText = TRUE, 
    mean = TRUE, 
    stat = TRUE,
  )

  expect_true(all(colnames(tt) %in% c(
    "ENTREZID", "logFC", "P.Value", "adj.P.Val", "mean", "stat", "coef", 
    "featureID", "hoverText", "comparison"
  )))
  
})

testthat::test_that("Column names are reformatted to limma style", {
  
  tt <- extractTopTables(
    input = res.edger, 
    coef = coefsModel[1], 
    columns = c("ENTREZID", "logFC", "P.Value", "adj.P.Val"),
    mean = TRUE, 
    stat = TRUE,
  )
  
  expect_true(all(colnames(tt) %in% c(
    "ENTREZID", "logFC", "P.Value", "adj.P.Val", "mean", "F", "coef", 
    "featureID", "comparison"
  )))
  
})

testthat::test_that("The input is must be a model or top table", {
  
  expect_error(
    extractTopTables(
      input = topTableList$B.LvsP[, c("ENTREZID", "logFC")], 
      coef = coefsModel[1], 
      columns = c("ENTREZID", "logFC")
    ),
    "'input' must be an object of class"
  )
  
})


