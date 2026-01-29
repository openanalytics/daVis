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
