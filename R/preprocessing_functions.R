# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

#' Convert Matrix to SummarizedExperiment Object
#'
#'
#' @param matrix
#' @param row.anno
#' @param col.data
#'
#' @return
#'
#' @import S4Vectors
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#'
nmfExperimentFromMatrix <- function(matrix, row.anno = NULL, col.data = NULL) {
  # Convert to matrix
  if(!is.matrix(matrix)) matrix <- as.matrix(matrix)
  # If row or/and column annotation data is missing use data given from matrix
  if(is.null(col.data)) {
    col.data <- colnames(matrix)
    if(is.null(col.data)) col.data <- 1:ncol(matrix)
    col.data <- DataFrame("samples" = col.data)
  }
  if(is.null(row.anno)) {
    row.anno <- rownames(matrix)
    if(is.null(row.anno)) row.anno <- 1:nrow(matrix)
    row.anno <- DataFrame("row.ids" = row.anno)
  }
  # Convert matrix to summarized experiment object.
  sum.exp <- SummarizedExperiment(assays = list("raw" = matrix),
                                  rowData = row.anno, colData = col.data)
  sum.exp <- nmfExperiment(sum.exp)
  return(sum.exp)
}

#' Reads Matrix & Annotation Files to SummarizedExperiment Object
#'
#' @param matrix.file
#' @param rowAnno.file
#' @param colData.file
#'
#' @return
#' @importFrom rtracklayer import.bed
#' @export
#'
#' @examples
nmfExperimentFromFile <- function(matrix.file, rowAnno.file = NULL,
                                  colData.file = NULL) {
  # Read matrix file
  matrix <- read.table(matrix.file, header = T)
  # Read row/col annotation if given.
  if(!is.null(rowAnno.file)) {
    if(grepl(rowAnno.file, pattern = ".bed")) {
      row.anno <- import.bed(rowAnno.file)
    } else {
      row.anno <- read.table(rowAnno.file)
    }
  } else row.anno <- NULL
  if(!is.null(colData.file)){
    col.data <- read.table(colAnno.file)
  } else {
    col.data <- NULL
  }
  # Convert to summarized experiment
  sum.exp <- nmfExperimentFromMatrix(matrix = matrix,
                                     row.anno = row.anno,
                                     col.data = col.data)
  return(sum.exp)
}

#' Title
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
normalizeUpperQuartile <- function(matrix) {
  matrix.norm <- apply(matrix, 2, function(c) {
    nf <- c[c != 0]
    c <- c / quantile(nf, 0.75)
    return(c)
  })
  return(matrix.norm)
}

#' Title
#'
#' @param col.vector
#' @param q
#'
#' @return
#' @export
#'
#' @examples
sigmoidTransform <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 2 / (1 + exp((-2) * col.vector / q)) - 1
  return(x)
}

#' Title
#'
#' @param col.vector
#' @param q
#'
#' @return
#' @export
#'
#' @examples
sigmoidTransform2 <- function(col.vector, q = 0.95) {
  q <- as.numeric(quantile(col.vector, q))
  x <- 1 / (1 + exp((-2) * col.vector / q))
  return(x)
}

#' Title
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
rankTransform  <- function(matrix) {
  trans.matrix <- apply(matrix, 2, function(c) {
    rank(c) / length(c)
  })
  return(trans.matrix)
}

#' Title
#'
#' @param coords
#'
#' @return
#'
#' @import GenomicRanges
#' @export
#'
#' @examples
igvCoord2Grange <- function(coords) {
  # Functions converts str "chr20:1-1,345,363" to GRange
  chrs <- sub(":.*", "", coords)
  start <- as.numeric(gsub("(,|-.*|.*:)", "", coords))
  end <- as.numeric(gsub("(,|.*-|.*:)", "", coords))
  gr <- GRanges(chrs, IRanges(start, end))
  return(gr)
}

#' Function to order binary data matrix
#'
#' @param matrix
#'
#' @return
#' @export
#'
#' @examples
orderBinary <- function(matrix) {
  col.sum <- apply(matrix, 2, sum)
  unlist(sapply(unique(col.sum), function(s) which(col.sum == s)))
}
