context("Get coefficient labels")

testthat::test_that("Test if labels are unique", {
  
  # coefficients
  coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
  # labels as vector
  coefsLabel1 <- LETTERS[seq_len(4)]
  # as a single function
  coefsLabel2 <- function(x) gsub("vs", " vs ", x)
  # as a list of vectors
  coefsLabel3 <- list(
    sub("(.+)\\.(.+)", "\\2", coefs),
    sub("(.+)\\.(.+)", "\\1", coefs)
  )
  # as a list of functions
  coefsLabel4 <- list(
    function(x) gsub("(.+)\\.(.+)", "\\2", x),
    function(x) gsub("(.+)\\.(.+)", "\\1", x)
  )
  coefsLabel <- list(coefsLabel1, coefsLabel2, coefsLabel3, coefsLabel4)
  
  expect_equal(
    object = {
      lapply(coefsLabel, getCoefLabel, coef = coefs)
    },
    expected = {
      list(
        LETTERS[seq_len(4)],
        setNames(c("B.L vs P", "L.L vs P", "B.P vs V", "L.P vs V"), coefs),
        list(
          c("LvsP", "LvsP", "PvsV", "PvsV"),
          c("B", "L", "B", "L")
        ),
        list(
          setNames(c("LvsP", "LvsP", "PvsV", "PvsV"), coefs),
          setNames(c("B", "L", "B", "L"), coefs)
        )
      )
    }
  )
  
})

testthat::test_that("An error is expected if the coefficient labels are duplicated", {
  
  # coefficients
  coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
  coefsLabel <- c("A", "A", "C", "D")
  
  # create the plot
  expect_error(
    getCoefLabel(coef = coefs, coefLabel = coefsLabel),
    "labels are not unique"
  )
  
})
