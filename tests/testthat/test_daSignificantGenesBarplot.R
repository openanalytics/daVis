context("Significant genes barplot")

library(ggplot2)

testthat::test_that("The bars correctly visualize the number of significant genes for a mixed input of limma and top table", {

  # create the plot
  gg <- daSignificantGenesBarplot(
    input = inputMixed,
    coef = coefMixed
  )
  
  # extract data behind the points
  idxGeomBar <- which(vapply(gg$layers, function(x) inherits(x$geom, "GeomBar"), logical(1)))
  ggDataBar <- ggplot2::layer_data(plot = gg, i = idxGeomBar)
  
  # check that number of significant genes correctly displayed for the model
  for(i in seq_along(coefsModel)){
    expect_equal(
      object = {
        ggDataCoef <- subset(ggDataBar, x == 1+!!i)
        list(
          label = ggDataCoef$label,
          value = ggDataCoef$y
        )
      },
      expected = {
        topTableCoef <- limma::topTable(fit = res.limma, coef = coefsModel[!!i], number = Inf, sort.by = "none")
        up <- nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05 & topTableCoef$logFC > 0), ])
        down <- nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05 & topTableCoef$logFC < 0), ])
        list(
          label = c(up, down),
          value = c(up, nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05), ]))
        )
      }
    )
  }
  
  # check that number of significant genes correctly displayed for the top table
  expect_equal(
    object = {
      ggDataCoef <- subset(ggDataBar, x == 1)
      list(
        label = ggDataCoef$label,
        value = ggDataCoef$y
      )
    },
    expected = {
      topTableCoef <- inputMixed$A
      up <- nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05 & topTableCoef$logFC > 0), ])
      down <- nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05 & topTableCoef$logFC < 0), ])
      list(
        label = c(up, down),
        value = c(up, nrow(topTableCoef[which(topTableCoef$adj.P.Val <= 0.05), ]))
      )
    }
  )
  
})
