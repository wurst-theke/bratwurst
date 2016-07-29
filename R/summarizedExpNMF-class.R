#' Title
#'
#' @slot HMatrixList list. 
#' @slot HMatrix matrix. 
#' @slot WMatrixList list. 
#' @slot WMatrix matrix. 
#' @slot FrobError DataFrame. 
#'
#' @return
#' 
#' @import SummarizedExperiment
#' @export
#'
#' @examples
summarizedExpNMF <- setClass(Class = 'summarizedExpNMF', 
                             contains = 'SummarizedExperiment',
                             representation = representation(HMatrixList = 'list',
                                                             WMatrixList = 'list',
                                                             FrobError = 'DataFrame'))
#==============================================================================#
#                        Getter & Setter                                       #
#==============================================================================#
# FrobError
setGeneric('FrobError', function(x, ...) standardGeneric('FrobError'))

setMethod('FrobError', 'summarizedExpNMF',
          function(x, ...) x@FrobError
)

# H-Matrix
setGeneric('HMatrixList', function(x, ...) standardGeneric('HMatrixList'))

setMethod('HMatrixList', 'summarizedExpNMF',
          function(x, ...) x@HMatrixList
)

# W-Matrix
setGeneric('WMatrixList', function(x, ...) standardGeneric('WMatrixList'))

setMethod('WMatrixList', 'summarizedExpNMF',
          function(x, ...) x@WMatrixList
)


