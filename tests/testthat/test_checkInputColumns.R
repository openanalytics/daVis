context("Check input columns")

testthat::test_that("Check if correct messages are printed", {
  
  # test length of featureIdVar (should be always of length 1)
  expect_error(
    checkColumns(
      input = res.limma, 
      featuresIdVar = c("SYMBOL", "gene_id"), 
      coef = c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV"), 
      cols = c("logFC", "adj.P.Val"),
      error = TRUE
    ),
    "'featuresIdVar' must indicate only one column."
  )
  
  # test part isModel(input) with wrong columns
  expect_error(
    checkColumns(
      input = res.limma, 
      featuresIdVar = "SYMBOL", 
      coef = "B.LvsP", 
      cols = c("logFC", "adj.P.Val", "col1", "col2"),
      error = TRUE
    ),
    regexp = "The column(s): 'col1', 'col2' must be present in 'input'.",
    fixed = TRUE
  )
  
  # test part if(is.list(input))
  expect_error(
    checkColumns(
      input = topTableList, 
      featuresIdVar = "SYMBOL", 
      coef = "B.LvsP", 
      cols = c("logFC", "adj.P.Val", "col1", "col2"),
      error = TRUE
    ),
    regexp = "The column(s): 'col1', 'col2' must be present in 'input'.",
    fixed = TRUE
  )
  
})
