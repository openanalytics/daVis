#' Upset plot for up- or down-regulated genes
#'
#' This is a function to create a customized upset plot for up- or 
#' down-regulated genes. 
#' @param coef character, coefficient names.
#' @param coefLabel character vector of labels or a function to transform 
#' existing labels.
#' @param featuresIdVar column name with unique feature ids, empty by default.
#' @param fdr threshold considered for significance, 0.05 by default.
#' @param dir direction for feature regulation ('up' to select up-regulated 
#' features; 'down' to select down-regulated features)
#' @param ylab y-axis title 
#' @param xlab x-axis title 
#' @param axesCex cex for the axis text. If two numbers provided, the first 
#' one is used for the x-axis and the second one for y-axis.
#' @param axesTitleCex cex for the axis title. If two numbers provided, 
#' the first one is used for the x-axis and the second one for y-axis.
#' @param barsCex cex for the counts above the bars.
#' @param returnAnalysis logical, if TRUE (FALSE by default), return also the 
#' output of the analysis (list with all overlapping sets based on 
#' \code{featuresIdVar}), otherwise only the plot object
#' @inheritParams createDataUpsetPlot
#' @importFrom UpSetR upset
#' @return uspet plot; if \code{returnAnalysis} is TRUE, return a list with 
#' overlapping sets and plot
#' @examples 
#' exampleData <- createExampleData(path = ".", output = "limma")
#' model <- exampleData$limma
#' coefs <- c("B.LvsP", "L.LvsP", "B.PvsV", "L.PvsV")
#' 
#' # Significantly up-regulated genes
#' daUpset(input = model, coef = coefs, fdr = 0.05, dir = "up")
#'  
#' # see vignette for other examples
#' 
#' @author Kirsten Van Hoorde, Katarzyna Gorczak, Michela Pasetto
#' @export
daUpset <- function(
  input,
  coef = NULL, 
  coefLabel = coef,
  featuresIdVar = character(),
  fdr = 0.05, 
  dir = c("up", "down"),
  ylab = "Intersection Size",
  xlab = "Set Size",
  axesCex = c(1.2, 1.7),
  axesTitleCex = c(1.4, 1.9),
  barsCex = 1.5,
  returnAnalysis = FALSE
) {
  
  dir <- match.arg(dir)
  nCoefs <- length(coef)
  if (nCoefs < 2) stop("At least 2 'coef' must be provided.", call. = FALSE)
  
  # Create data for visualization
  dataFrameSG <- createDataUpsetPlot(
    input = input, coef = coef, coefLabel = coefLabel,
    featuresIdVar = featuresIdVar, fdr = fdr, dir = dir
  )
  
  if (returnAnalysis) sets <- extractFeatures(dataFrameSG)
  
  nBars <- nrow(unique(dataFrameSG))
  
  if(any(rowSums(dataFrameSG) > 1)){
    queryList   <- extractQueryList(dataFrameSG)		
    queriesPlot <- unlist(queryList, recursive = FALSE)
    queriesPlot <- queriesPlot[!unlist(lapply(queriesPlot, is.null))]
  }else  queriesPlot <- NULL
  
  pointSize <- ifelse(
    nBars < 5, 6, ifelse(nBars < 10, 3, ifelse(nBars < 20, 1.5, 0.75)))
  lineSize  <- ifelse(
    nBars < 5, 2, ifelse(nBars < 10, 1.3, ifelse(nBars < 20, 0.6, 0.1)))
  namesSize <- ifelse(
    nBars < 5, 2, ifelse(nBars < 10, 1.75, ifelse(nBars < 20, 1.5, 1)))
  
  plot <- upset(
    dataFrameSG,  
    nsets = ncol(dataFrameSG),
    nintersects = NA,
    main.bar.color = "black",
    queries = queriesPlot,
    order.by = c("degree", "freq"), 
    decreasing = c(TRUE, TRUE),
    mainbar.y.label = ylab,
    sets.x.label = xlab,
    line.size = lineSize, 
    text.scale = c(
      # Intersection size title
      ifelse(length(axesTitleCex) > 1, axesTitleCex[2], axesTitleCex), 
	  # intersection size tick labels
      ifelse(length(axesCex) > 1, axesCex[2], axesCex), 
	  # set size title
      ifelse(length(axesTitleCex) > 1, axesTitleCex[1], axesTitleCex),
	  # set size tick labels
      ifelse(length(axesCex) > 1, axesCex[1], axesCex),
	  # set names
      namesSize,
	  # numbers above bars
      barsCex),
    point.size = pointSize
  )
  
  if (returnAnalysis) {
    list(sets = sets,
         plot = plot)
  } else plot
  
}


#' Create data for upset plot
#' @inheritParams extractTopTables
#' @inheritParams daUpset
#' @return data.frame compatible with UpSetR
#' @author Kirsten Van Hoorde, Katarzyna Gorczak
createDataUpsetPlot <- function(
  input, 
  coef, coefLabel = NULL,
  featuresIdVar = character(), 
  fdr, 
  dir
) {
  
  dirPosNeg <- switch(dir,
                      "up" = "pos",
                      "down" = "neg")
  
  columns <- unique(c(featuresIdVar, "logFC", "adj.P.Val"))
  tblList <- extractTopTables(
    input = input, coef = coef, coefLabel = coefLabel, features = NULL, 
    n = Inf, columns = columns, hoverText = FALSE, output = "list"
  )
  if (length(featuresIdVar) == 0) featuresIdVar <- "featureID"
  tbl <- arrangeTopTables(
    input = tblList, fdr = fdr, commonFeatures = FALSE, dir = dirPosNeg, 
    featuresIdVar = featuresIdVar, output = "list"
  )
  
  # exclude contrasts with no significant features
  tbl <- tbl[which(lapply(tbl, nrow) > 0)]
  
  if(!sum(vapply(tbl, nrow, numeric(1))) > 0)
    stop("No features are ", 
         switch(dir, "up" = "up-regulated", "down" = "down-regulated"), ".",
         call. = FALSE)
  
  coefLabels <- unname(vapply(tbl, function(x) 
    unique(as.character(x[, "comparison"])), character(1)))
  tbl <- lapply(tbl, function(x) as.character(x[, featuresIdVar]))
  names(tbl) <- coefLabels
  
  if(sum(vapply(tbl, length, numeric(1)) != 0) == 1){
    idx <- which(vapply(tbl, length, numeric(1)) != 0)
    stop("There are up-regulated significant features only for one",
		" coefficient: ", names(tbl)[idx], ". The plot is not generated.",
		call. = FALSE)
  }
  
  elements <- unique(unlist(tbl))
  dt <- unlist(lapply(tbl, function(x) {
    as.vector(match(elements, x))
  }))
  dt[is.na(dt)] <- as.integer(0)
  dt[dt != 0] <- as.integer(1)
  dt <- data.frame(matrix(dt, ncol = length(tbl), byrow = FALSE))
  idx <- which(is.na(elements))
  if (length(idx) > 0) elements[idx] <- paste0("NA.", seq_len(length(idx)))
  rownames(dt) <- elements
  dt <- dt[which(rowSums(dt) != 0), ]
  names(dt) <- names(tbl)
  
  dt
}


#' Extract query list
#' 
#' @param dataFrameSG data.frame compatible with UpSetR
#' @importFrom grDevices colorRampPalette
#' @importFrom utils combn
#' @importFrom UpSetR intersects
#' @return a list with queries for upset plot
extractQueryList <- function(dataFrameSG) {
  
  lapply(seq(2, ncol(dataFrameSG)), function(n){
    
    cols <- colorRampPalette(
      c("darkblue", "blue", "steelblue3", "dodgerblue"))(ncol(dataFrameSG))
    dataFrameOI <- dataFrameSG[rowSums(dataFrameSG) == n,, drop = FALSE]
    sets <- combn(colnames(dataFrameSG), m = n)
    
    apply(sets, 2, function(s){
      if(any(rowSums(dataFrameOI[, s]) == n)){
        list(query = intersects,
             params = as.list(s),
             color = cols[n],
             active = TRUE)
      }else{
        NULL
      }
    })
  })
  
}

#' Extract feature ids for each overlapping set
#' 
#' @param data data.frame with 0s and 1s; 
#' the number of columns corresponds to the number of coefs and the number of 
#' rows corresponds to the number of unique feature ids
#' each column is filled with 0s or 1s depending on the overlapping set
#' @author Katarzyna Gorczak
#' @importFrom utils combn
#' @return list with features in overlapping sets
extractFeatures <- function(data){
  
  res <- lapply(seq_len(ncol(data)), function(n){
    
    dataFrameOI <- data[rowSums(data) == n,, drop = FALSE]
    sets <- combn(colnames(data), m = n)
    
    tmp <- lapply(seq_len(dim(sets)[2]), function(s){
      s <- sets[, s]
      if(any(rowSums(dataFrameOI[, s, drop = FALSE]) == n)){
        names(which(rowSums(dataFrameOI[, s, drop = FALSE]) == n) == TRUE)
      }else{
        NULL
      }
    })
    
    names(tmp) <- apply(sets, 2, paste, collapse = "; ")
    idx <- which(vapply(tmp, is.null, logical(1)) == FALSE)
    tmp[idx]
    
  })
  
  names(res) <- c(
    paste0("Intersect ", 1, " set"), 
    paste0("Intersect ", seq_len(ncol(data))[-1], " sets")
  )
  
  idx <- which(vapply(res, length, numeric(1)) > 0)
  res[idx]
  
}
