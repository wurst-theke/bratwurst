#' NMF Experiment Class
#'
#' @slot HMatrixList list. 
#' @slot HMatrix list. 
#' @slot WMatrixList list. 
#' @slot WMatrix list. 
#' @slot FrobError DataFrame. 
#' @slot OptKStats DataFrame. 
#'
#' @return
#' 
#' @import SummarizedExperiment
#' @export
#'
#' @examples
nmfExperiment <- setClass(Class = 'nmfExperiment', 
                          contains = 'SummarizedExperiment',
                          representation = representation(HMatrixList = 'list',
                                                          WMatrixList = 'list',
                                                          FrobError = 'DataFrame',
                                                          OptKStats = 'DataFrame',
                                                          OptK = 'numeric', 
                                                          FeatureStats = 'DataFrame'))

#==============================================================================#
#                                 Getter & Setter                              #
#==============================================================================#
#### FrobError
# Getter
setGeneric('FrobError', function(x, ...) standardGeneric('FrobError'))

#' Frobenius Error getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('FrobError', 'nmfExperiment', function(x, ...) x@FrobError)

# Setter
setGeneric('setFrobError', function(nmfExperiment, FrobError) standardGeneric('setFrobError'))

#' Frobenius Error setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setFrobError', 'nmfExperiment', function(nmfExperiment, FrobError) {
  nmfExperiment@FrobError <- FrobError
  return(nmfExperiment)
})

#### H-Matrix List
# Getter
setGeneric('HMatrixList', function(x, k = NULL, ...) standardGeneric('HMatrixList'))

#' H-Matrix List getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('HMatrixList', 'nmfExperiment', function(x, k = NULL, ...){
  if(is.null(k)) {
    x@HMatrixList
  } else {
    x@HMatrixList[[as.character(k)]]
  }
})

# Setter
setGeneric('setHMatrixList', function(nmfExperiment, HMatrixList) standardGeneric('setHMatrixList'))

#' H-Matrix List setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setHMatrixList', 'nmfExperiment', function(nmfExperiment, HMatrixList){
  nmfExperiment@HMatrixList <- HMatrixList 
  return(nmfExperiment)
})

#### W-Matrix
# Getter
setGeneric('WMatrixList', function(x, k = NULL, ...) standardGeneric('WMatrixList'))

#' W-Matrix list getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('WMatrixList', 'nmfExperiment', function(x, k = NULL, ...) {
  if(is.null(k)) {
    x@WMatrixList
  } else {
    x@WMatrixList[[as.character(k)]]   
  }
})

# Setter
setGeneric('setWMatrixList', function(nmfExperiment, WMatrixList) standardGeneric('setWMatrixList'))

#' W-Matrix setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setWMatrixList', 'nmfExperiment', function(nmfExperiment, WMatrixList) {
  nmfExperiment@WMatrixList <- WMatrixList 
  return(nmfExperiment)
})

#### H-Matrix (H-Matrix with smallest frobError) 
# Getter
setGeneric('HMatrix', function(x, k = NULL, ...) standardGeneric('HMatrix'))

#' H-Matrix getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('HMatrix', 'nmfExperiment', function(x, k = NULL, ...) {
  i.min <- apply(x@FrobError, 2, which.min)
  if(is.null(k)) {
    H <- lapply(names(x@HMatrixList), function(k) {
      x@HMatrixList[[k]][[i.min[k]]]
    })
    names(H) <- names(x@HMatrixList)
  } else {
    k <- as.character(k)
    H <- x@HMatrixList[[k]][[i.min[k]]]
  }
  return(H)
})

#### W-Matrix (W-Matrix with smallest frobError) 
# Getter
setGeneric('WMatrix', function(x, k = NULL, ...) standardGeneric('WMatrix'))

#' W-Matrix getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('WMatrix', 'nmfExperiment', function(x, k = NULL, ...) {
  i.min <- apply(x@FrobError, 2, which.min)
  if (is.null(k)) {
    W <- lapply(names(x@WMatrixList), function(k) {
      x@WMatrixList[[k]][[i.min[k]]]
    })
    names(W) <- names(x@WMatrixList)
  } else {
    k <- as.character(k)
    W <- x@WMatrixList[[k]][[i.min[k]]]
  }
  return(W)
})

#### Optimal K Statistics
# Getter
setGeneric('OptKStats', function(x, ...) standardGeneric('OptKStats'))

#' Optimal K Statistics getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('OptKStats', 'nmfExperiment', function(x, ...) x@OptKStats)

# Setter
setGeneric('setOptKStats', function(nmfExperiment, OptKStats) standardGeneric('setOptKStats'))

#' Optimal K Statistics setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setOptKStats', 'nmfExperiment', function(nmfExperiment, OptKStats) {
  nmfExperiment@OptKStats <- OptKStats 
  return(nmfExperiment)
})

#### Optimal K 
# Getter
setGeneric('OptK', function(x, ...) standardGeneric('OptK'))

#' Optimal K 
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('OptK', 'nmfExperiment', function(x, ...) x@OptK)

# Setter
setGeneric('setOptK', function(nmfExperiment, OptK) standardGeneric('setOptK'))

#' Optimal K setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setOptK', 'nmfExperiment', function(nmfExperiment, OptK) {
  nmfExperiment@OptK<- OptK
  return(nmfExperiment)
})

#### Feature Statistics
# Getter
setGeneric('FeatureStats', function(x, ...) standardGeneric('FeatureStats'))

#' Feature Statistics getter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('FeatureStats', 'nmfExperiment', function(x, ...) x@FeatureStats)

# Setter
setGeneric('setFeatureStats', function(nmfExperiment, FeatureStats) standardGeneric('setFeatureStats'))

#' Feature Statistics setter
#'
#' @param nmfExperiment 
#'
#' @return
#' @export
#'
#' @examples
setMethod('setFeatureStats', 'nmfExperiment', function(nmfExperiment, FeatureStats) {
  if(nrow(nmfExperiment@FeatureStats) == 0) {
    nmfExperiment@FeatureStats <- FeatureStats
  } else {
    nmfExperiment@FeatureStats <- cbind(nmfExperiment@FeatureStats, FeatureStats)
  }
  return(nmfExperiment)
})