context("Get features")

testthat::test_that("The features of interest are correctly extracted for input", {
  
  set.seed(123)
  features <- sample(res.limma$genes$ENTREZID, 10)
  
  expect_equal(
    object = {
      obj <- subsetFeatures(res.limma, features, featuresIdVar = "ENTREZID")
      obj$genes$ENTREZID
    },
    expected = {
      features
    }
  )
  
  expect_equal(
    object = {
      obj <- subsetFeatures(res.edger, features, featuresIdVar = "ENTREZID")
      obj$genes$ENTREZID
    },
    expected = {
      features
    }
  )
  
  expect_equal(
    object = {
      obj <- subsetFeatures(res.deseq, features, featuresIdVar = "ENTREZID")
      obj$ENTREZID
    },
    expected = {
      features
    }
  )
  
  expect_equal(
    object = {
      obj <- subsetFeatures(topTableList, features, featuresIdVar = "ENTREZID")
      unique(unlist(lapply(obj, function(x) x$ENTREZID)))
    },
    expected = {
      features
    }
  )
  
  expect_error(
    subsetFeatures(res.limma, features, featuresIdVar = "SYMBOL"),
    "No 'features' matching 'featuresIdVar'."
  )
  
})

testthat::test_that("The features are correctly extracted from input", {
  
  features <- sort(res.limma$genes$ENTREZID)
  
  expect_equal(
    object = {
      sort(getFeatures(res.limma, featuresIdVar = "ENTREZID"))
    },
    expected = {
      features
    }
  )
  
  expect_equal(
    object = {
      sort(getFeatures(topTableList, featuresIdVar = "ENTREZID"))
    },
    expected = {
      features
    }
  )
  
  expect_error(
    sort(getFeatures(res.limma, featuresIdVar = "SYMBOL")),
    "'featuresIdVar' must indicate unique values."
  )
  
})
