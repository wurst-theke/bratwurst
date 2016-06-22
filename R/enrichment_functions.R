#' Title
#'
#' @param peaks 
#' @param subregions 
#' @param universe 
#'
#' @return
#' @export
#'
#' @examples
computeRegionFC <- function(peaks, subregions, universe) {
  prop.universe <- sum(countOverlaps(universe, peaks))/length(universe)
  prop.subregions <- sum(countOverlaps(subregions, peaks))/length(subregions)
  fc <- prop.subregions/prop.universe
  return(fc)
}

# Compute FC.
#' Title
#'
#' @param peak 
#' @param subregion 
#' @param universe 
#'
#' @return
#' @export
#'
#' @examples
computeRegionFC4GRlist <- function(peak, subregion, universe) {
  # Compute proportion of peaks in universe.
  ov.universe <- findOverlaps(universe, peak)
  prop.universe <- table(subjectHits(ov.universe))/length(universe)
  # Compute proportion of peaks in subregions
  ov <- findOverlaps(subregion, peak)
  prop.subregions <- table(subjectHits(ov))/length(subregion)
  # Find Peaks with no-overlap to add them later.
  peak.ids <- c(names(prop.universe), names(prop.subregions))
  i.common <- peak.ids[duplicated(peak.ids)]
  zero.ids <- as.character(which(!as.character(1:length(peak))%in%i.common))
  # Compute Fold-change.
  fc <- prop.subregions[i.common]/prop.universe[i.common]
  fc[zero.ids] <- NA
  fc <- fc[order(as.numeric(names(fc)))]
  return(fc)
}

#' Title
#'
#' @param peaks 
#' @param subregions 
#' @param universe 
#'
#' @return
#' @export
#'
#' @examples
computeFisher4Regions <- function(peaks, subregions, universe) {
  x1 <- sum(countOverlaps(subregions, peaks) > 0)
  x2 <- length(subregions) - x1
  x3 <- sum(countOverlaps(universe, peaks) > 0) - x1
  x4 <- length(universe) - x3 - x2 - x1 
  m <- matrix(c(x1, x3, x2, x4), nrow = 2)
  p <- fisher.test(m)
  return(p$p.value)
}

### NES Functions
#' Title
#'
#' @param ranks 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
computeRecoveryAUC <- function(ranks, n) {
  ranks <- unique(c(1, sort(ranks[ranks < n]), n))
  rank.diff <- ranks[-1] - ranks[-length(ranks)] 
  #recov <- 0:c(length(rank.diff)-1)
  recov <- c(0, cumsum(rep(1/n, length(rank.diff)-1)))
  auc <- sum(recov*rank.diff)/n 
  return(auc)
}

#' Title
#'
#' @param ranks 
#' @param n 
#' @param iter 
#' @param threads 
#'
#' @return
#' @export
#'
#' @examples
computeRmdAUCs <- function(ranks, n, iter = 10^3, threads = 3) {
  n.hits <- sum(ranks <= n)
  rmd.aucs <- mclapply(1:iter, function(i) {
    rmd.ranks <- sort(sample(1:n, size = n.hits, replace = F))
    rmd.auc <- computeRecoveryAUC(rmd.ranks, n)
    return(rmd.auc)
  }, mc.cores = threads)
  return(unlist(rmd.aucs))
}

#' Title
#'
#' @param ranks 
#' @param n 
#' @param iter 
#' @param threads 
#'
#' @return
#' @export
#'
#' @examples
computeNES <- function(ranks, n, iter = 10^3, threads = 3) {
  auc <- computeRecoveryAUC(ranks, n)
  rmd.aucs <- computeRmdAUCs(ranks, n, iter = iter, threads = threads)
  nes <- (auc - mean(rmd.aucs))/sd(rmd.aucs)
  r <- sum(auc < rmd.aucs) + 1
  result <- data.frame('Rank' = r, 'NES' = nes, 'AUC' = auc)
  return(result)
}

#' Title
#'
#' @param peaks 
#' @param subregions 
#' @param universe 
#'
#' @return
#' @export
#'
#' @examples
computePeakStats4Regions <- function(peaks, subregions, universe) {
  n.regions <- lapply(subregions, length)
  n.ov <- lapply(subregions, function(r) sum(countOverlaps(r, peaks) > 0))
  # Compute Fold-enrichment
  fcs <- lapply(subregions, function(subregion){
    if(length(subregion) == 0) {return(NA)}
    computeRegionFC(peaks, subregion, universe)
  })
  # Compute Fisher-exact p-value
  pvalues <- lapply(subregions, function(subregion){
    if(length(subregion) == 0) {return(NA)}
    computeFisher4Regions(peaks, subregion, universe)
  })
  # Compute Normalized enrichment score
  nes <- lapply(subregions, function(subregion) {
    if(length(subregion) <= 1) {return(NA)}
    ov <- findOverlaps(subregion, peaks)
    promoter.rank <- unique(queryHits(ov))
    nes <- computeNES(promoter.rank, length(subregion))
    return(nes) 
  })
  nes <- do.call(rbind, nes)
  result  <- data.frame('NumberRegions' = unlist(n.regions), 'OV' = unlist(n.ov),
                        'FC' = unlist(fcs), 'FisherPvalue' = unlist(pvalues))
  result <- cbind(result, nes)
  return(result)
}

#' Title
#'
#' @param peaks 
#' @param subregions 
#' @param universe 
#'
#' @return
#' @export
#'
#' @examples
computeRunningPeakStats4Regions <- function(peaks, subregions, universe) {
  result.rankProp <- lapply(seq(0.1, 1, by = 0.1), function(region.prop) {
    # Get only top x% of region.
    filtered.subregions <- lapply(subregions, function(r) {
      n <- floor(length(r)*region.prop)
      return(r[0:n,])
    })
    result <- computePeakStats4Regions(peaks, filtered.subregions, universe)
    result$ClusterID <- 1:2
    result$RankProp <- region.prop
    return(result)
  })
  result.rankProp <- do.call(rbind, result.rankProp)
  return(result.rankProp)
}