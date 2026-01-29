library(testthat)
library(daVis)

# example data for the test
tmpDir <- tempfile();dir.create(tmpDir)
exampleData <- createExampleData(path = tmpDir, output = c("limma", "topTable", "edgeR"))
res.limma <- exampleData$limma
res.edger <- exampleData$edgeR
topTableList <- exampleData$topTable

# consider a subset of the features and coefficients for the tests
set.seed(123)
features <- sample(rownames(res.limma), 10)
coefsModel <- sample(colnames(res.limma), 2)

# mixed input
inputMixed <- list(res.limma, A = topTableList[[1]])
coefMixed <- c("A", coefsModel)

# different set of features
res.limmaRep <- res.limma
res.limmaRep$genes[match(features, res.limmaRep$genes$ENTREZID), "SYMBOL"] <- "G"
inputMixedRep <- list(res.limma, A = topTableList[[1]][setdiff(features, c("66292", "68585")), ])
coefMixed <- c("A", coefsModel)
