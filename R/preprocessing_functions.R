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
    c <- c/quantile(nf, 0.75)
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
  x <- 2/(1 + exp(-2*col.vector/q)) - 1
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
  x <- 1/(1 + exp(-2*col.vector/q))
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
    rank(c)/length(c)
  })
  return(trans.matrix)
}

#' Title
#'
#' @param coords 
#'
#' @return
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