#' Create example limma-, edgeR-object or a list of top tables.
#' 
#' @param path absolute path where data should be downloaded to
#' @param output which object should be created. Possible options: 
#' limma, edgeR, topTableList, DESeq2
#' @param quiet if TRUE, suppress status messages (if any), and the progress bar
#' (see \code{download.file} function)
#' @return list with specific objects (output from \code{limma}, \code{edgeR},
#' \code{topTable} and/or \code{DESeq2})
#' @examples 
#' tmpDir <- tempfile(); dir.create(tmpDir)
#' getData <- createExampleData(path = tmpDir, output = "limma")
#' 
#' @export
createExampleData <- function(
    path = ".", 
    output = c("limma", "topTable", "edgeR", "deseq2"),
    quiet = TRUE
) {
  
    output <- match.arg(output, several.ok = TRUE)
    requireNamespace("tools")
    path <- tools::file_path_as_absolute(path)
    
    downloadData(path = path, quiet = quiet)
    eset <- createExpressionSet(path = path)
    objectList <- modelExampleData(eset = eset)
    
    if ("edgeR" %in% output) res.edgeR <- runEdgeR(input = objectList)
    if ("limma" %in% output) res.limma <- runLimma(input = objectList)
    if ("deseq2" %in% output) res.deseq <- runDESeq2(input = objectList)
    if ("topTable" %in% output) topTableList <- runTopTable(input = objectList)
    
    # unlink("GSE60450_Lactation-samples.txt.gz")
    # unlink("GSE60450_Lactation-GenewiseCounts.txt.gz")
    
    c(
      if ("limma" %in% output) list(limma = res.limma),
      if ("topTable" %in% output) list(topTable = topTableList),
      if ("edgeR" %in% output) list(edgeR = res.edgeR),
      if ("deseq2" %in% output) list(DESeq2 = res.deseq)
    )
  
} 

#' Download example data from 
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450
#' @inheritParams createExampleData
#' @author Katarzyna Gorczak
#' @importFrom utils download.file
#' @return (invisibly) downloads files at specific location ('path')
downloadData <- function(path, quiet) {
  countsURL <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn",
    "/GSE60450/suppl/GSE60450%5FLactation%2DGenewiseCounts%2Etxt%2Egz"
  )
  download.file(
    url = countsURL, 
    destfile = paste0(path, "/GSE60450_Lactation-GenewiseCounts.txt.gz"),
    quiet = quiet
  )
  
  samplesURL <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn", 
    "/GSE60450/suppl/filelist.txt"
  )
  download.file(
    url = samplesURL, 
    destfile = paste0(path, "/GSE60450_Lactation-samples.txt.gz"),
    quiet = quiet
  )
}

#' Create sample annotation from example data
#' @inheritParams createExampleData
#' @importFrom utils read.table
#' @return data.frame with sample annotation
#' @author Katarzyna Gorczak
createSampleAnnotation <- function(path) {
  samples <- read.table(
    gzfile(paste0(path, "/GSE60450_Lactation-samples.txt.gz")),
    header = FALSE, skip = 2)
  samples$GEO <- unlist(lapply(samples$V2, function(x) 
    unlist(strsplit(x, "_"))[1]))
  samples$SampleName <- gsub(".txt.gz", "", samples$V2)
  samples$SampleName <- gsub("^[^_]*_(.*)", "\\1", 
                             gsub(".txt.gz", "", samples$SampleName))
  samples$SampleName <- gsub("-", ".", samples$SampleName)
  
  GSEannotation <- c(
    "GSM1480291", "Luminal virgin", "GSM1480292", "Luminal virgin", 
    "GSM1480293", "Luminal 18.5 dP", "GSM1480294", "Luminal 18.5 dP", 
    "GSM1480295", "Luminal 2 dL", "GSM1480296", "Luminal 2 dL", "GSM1480297", 
    "Basal virgin", "GSM1480298", "Basal virgin", "GSM1480299", "Basal 18.5 dP",
    "GSM1480300", "Basal 18.5 dP", "GSM1480301", "Basal 2 dL", "GSM1480302", 
    "Basal 2 dL")
  GEOs <- GSEannotation[seq(1, length(GSEannotation), 2)]
  type <- GSEannotation[seq(2, length(GSEannotation), 2)]
  
  if (all(samples$GEO == GEOs)) {
    samples$Type <- ifelse(grepl("Luminal", type), "L", "B")
    samples$Status <- ifelse(grepl("virgin", type), "virgin", 
                             ifelse(grepl("dP", type), "pregnant", "lactating"))
  }
  
  samples
}

#' Create ExpressionSet object from example data
#' @inheritParams createExampleData
#' @importFrom utils read.table
#' @return ExpressionSet object
#' @author Katarzyna Gorczak
createExpressionSet <- function(path) {
  counts <- read.table(
    gzfile(paste0(path, "/GSE60450_Lactation-GenewiseCounts.txt.gz")), 
    header = TRUE)
  samples <- createSampleAnnotation(path = path)
  samples <- samples[match(colnames(counts)[-c(1, 2)], samples$SampleName), ]
  rownames(samples) <- unlist(lapply(samples$SampleName, function(x) 
    unlist(strsplit(x, "_"))[1]))
  
  entrezid <- as.character(counts[, 1])
  requireNamespace("AnnotationDbi", quietly = TRUE)
  suppressPackageStartupMessages(requireNamespace("org.Mm.eg.db"))
  genes <- vapply(c("SYMBOL", "GENENAME", "ENSEMBL"), function(x) 
    suppressMessages( # don't show returned 1:1 or 1:many mapping between ...
      AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = entrezid, 
                            keytype = "ENTREZID", column = x)
    ), 
    character(nrow(counts)))
  genes <- as.data.frame(cbind(ENTREZID = entrezid, genes))
  counts <- as.matrix(counts[, -c(1, 2)])
  rownames(genes) <- rownames(counts) <- genes$ENTREZID
  colnames(counts) <- rownames(samples)
  
  requireNamespace("Biobase", quietly = TRUE)
  requireNamespace("methods")
  eset <- Biobase::ExpressionSet(
    assayData = counts, 
    phenoData = methods::new("AnnotatedDataFrame", data = samples),
    featureData = methods::new("AnnotatedDataFrame", data = genes)
  )
  
  eset
}

#' Normalize data, create model matrix and contrasts
#' @param eset ExpressionSet
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts
#' @return a list with data, model, design matrix and contrast matrix
#' @author Katarzyna Gorczak
modelExampleData <- function(eset) {
  
  if(!inherits(eset, what = 'ExpressionSet'))
    stop("'eset' must be of class 'ExpressionSet'.")
  
  group <- factor(paste(eset$Type, eset$Status, sep = "."))
  requireNamespace("Biobase", quietly = TRUE)
  y <- DGEList(counts = Biobase::exprs(eset), group = group, 
               genes = Biobase::fData(eset))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  contrasts <- c(B.LvsP = "B.lactating-B.pregnant",
                 L.LvsP = "L.lactating-L.pregnant",
                 B.PvsV = "B.pregnant-B.virgin",
                 L.PvsV = "L.pregnant-L.virgin",
                 LvsB = "(L.lactating-L.pregnant)-(B.lactating-B.pregnant)")
  contrastsOfInterest <- do.call(
    makeContrasts, c(contrasts, list(levels = design))
  )
  
  list(
    "contrasts" = contrastsOfInterest, 
    "data" = y,
    "design" = design,
    "eset" = eset
  )
}

#' Fit the model with edgeR
#' @param input list of objects (data, contrasts and model matrix)
#' @importFrom edgeR glmQLFit glmQLFTest
#' @return edgeR output from glmQLFTest() function
#' @author Katarzyna Gorczak
runEdgeR <- function(input) {
  
  y <- input$data
  design <- input$design
  contrast <- input$contrasts
  
  fit <- glmQLFit(y, design, robust = TRUE)
  fit$contrasts <- contrast
  res.edgeR <- glmQLFTest(fit, contrast = contrast[, seq_len(4)])
  # topTags(res.edgeR)
  
  res.edgeR
}

#' Fit the model with limma
#' @param input list of objects (data, contrasts and model matrix)
#' @importFrom limma voom lmFit contrasts.fit eBayes
#' @return limma outpur from eBayes() function
#' @author Katarzyna Gorczak
runLimma <- function(input) {
  
  y <- input$data
  design <- input$design
  contrast <- input$contrasts
  
  voom.data <- voom(y, design = design)
  voom.fit <- lmFit(voom.data, design)
  voom.fit <- contrasts.fit(fit = voom.fit, contrasts = contrast[, seq_len(4)])
  res.limma <- eBayes(voom.fit)
  # topTable(res.limma)
  
  res.limma
}

#' Fit the model with DESeq2
#' @param input list of objects (data, contrasts and model matrix)
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @return DESeq2 output from results() function
#' @author Katarzyna Gorczak
runDESeq2 <- function(input) {
  
  y <- input$data
  design <- input$design
  contrast <- input$contrasts
  eset <- input$eset
  
  requireNamespace("Biobase", quietly = TRUE)
  dds <- DESeqDataSetFromMatrix(
    countData = y$counts, colData = Biobase::pData(eset), design = design
  )
  dds <- suppressMessages(DESeq(dds))
  requireNamespace("S4Vectors", quietly = TRUE)
  S4Vectors::mcols(dds) <- cbind(S4Vectors::mcols(dds), y$genes)
  res.deseq <- results(dds, contrast = contrast[, 1])
  for (iCol in colnames(y$genes)) {
    res.deseq[[iCol]] <- y$genes[, iCol]
  }
  if (all(rownames(res.deseq) == res.deseq$ENTREZID)) {
    res.deseq 
  } else  stop("'rownames(res.deseq)' do not match 'ENTREZID'")
  
  # # all below give the same results
  # res1 <- results(dds) # results for the last column in 'design' ('L.virgin')
  # res2 <- results(dds, name = "L.virgin")
  # res3 <- results(dds, contrast = c(0,0,0,0,0,1)) # round to 9 digits
  
  res.deseq
}

#' Fit the model and return top table
#' @param input list of objects (data, contrasts and model matrix)
#' @importFrom limma voom lmFit contrasts.fit eBayes topTable
#' @return a list of top tables
#' @author Katarzyna Gorczak
runTopTable <- function(input) {
  
  y <- input$data
  design <- input$design
  contrast <- input$contrasts
  
  voom.data <- voom(y, design = design)
  voom.fit <- lmFit(voom.data, design)
  voom.fit <- contrasts.fit(fit = voom.fit, contrasts = contrast[, seq_len(4)])
  res.limma <- eBayes(voom.fit)
  
  topTableList <- lapply(colnames(contrast[, seq_len(4)]), function(x){
    topTable <- topTable(fit = res.limma, coef = x, number = Inf)
    topTable[, "se"] <- getSEModel(input = res.limma, coef = x)
    return(topTable)
  })
  names(topTableList) <- colnames(contrast[, seq_len(4)])
  
  topTableList
}