#' Common parameters for the functions of the \code{daVis} package.
#' @param input model or a list with top tables named with coefficients.\cr
#' For model: object of class \code{MArrayLM} 
#' (linear model, see \code{\link[limma]{eBayes}}), \code{DGELRT} 
#' (see \code{edgeR}) or \code{DESeqResults} (see \code{DESeq2}) are supported.
#' See more details on the top table format in section:
#' \link[=daVis-common-doc]{Top table format}.
#' @name daVis-common-args
#' @return No returned value
NULL

#' Common documentation for the \code{daVis} package.
#' @section Top table format:
#' The top table extracted from a specified model, or specified by the user
#' contains at least the columns:
#' \itemize{
#' \item{'logFC': log fold change}
#' \item{'AveExpr': average expression}
#' \item{'P.Value' and 'adj.PVal': raw and (multiplicity correction) adjusted 
#' p-values
#' \itemize{
#' \item{For \code{edgeR} model: this is extracted from the columns: 
#' 'PValue' and 'FDR' respectively.}
#' }}
#' \item{feature identifier as specified by: \code{featuresIdVar}}
#' }
#' @name daVis-common-doc
#' @return No returned value
NULL

