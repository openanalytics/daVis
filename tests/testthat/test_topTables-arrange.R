context("Arrange top tables")

testthat::test_that("Reshape long-form table", {
  
  dtList <- reshapeTable(ttLong, "coef")
  
  expect_identical(
    object = class(dtList),
    expected = "list"
  )
  
})


testthat::test_that("Check if top tables are correctly arranged", {
  
  dtList <- arrangeTopTables(
    ttLong, 
    featuresIdVar = NULL, 
    logFCrange = c(-2.5, 2.5),
    output = "list"
  )
  ttSubset <- lapply(topTableList, function(x) {
    x[x$logFC >= -2.5 & x$logFC <= 2.5, ]
  })
  
  expect_equal(
    object = sapply(dtList, nrow),
    expected = sapply(ttSubset, nrow)
  )
  
  expect_error(
    arrangeTopTables(
      topTableList, 
      featuresIdVar = NULL, 
      logFCrange = c(7, 10),
      output = "list"
    ),
    "There are no genes with logFC in range"
  )
  
  ttSubset$L.PvsV$ENTREZID <- "G"
  expect_error(
    arrangeTopTables(
      ttSubset, 
      featuresIdVar = "ENTREZID", 
      commonFeatures = TRUE, 
      output = "list"
    ),
    regex = "There are no common genes for B.LvsP, L.LvsP, B.PvsV, L.PvsV.",
    fixed = TRUE
  )
  
})

testthat::test_that("Chech if elements are made unique", {
  
  set.seed(123)
  vec <- sample(LETTERS[1:10], 15, replace = TRUE)
  
  expect_identical(
    object = length(unique(makeElementsUnique(vec))),
    expected = length(unique(paste0(vec, seq_len(length(vec)))))
  )
  
})
