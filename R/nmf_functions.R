#==============================================================================#
#                         NMF GPU Wrapper - FUNCTIONS                          #
#==============================================================================#
#' Title
#'
#' @param matrix 
#' @param tmp.path 
#' @param sep 
#'
#' @return
#' @export
#'
#' @examples
writeTmpMatrix <- function(matrix, tmp.path = '/tmp/nmf_tmp', sep = ' ') {
  # Write matrix to tmp file.
  matrix <- round(matrix, 7)
  dir.create(tmp.path)
  # Rmv final slash from path --> will lead to an error of nmfGpu wrapper
  tmp.path <- sub('/$', '', tmp.path)
  tmpMatrix.path <- tempfile(tmpdir = tmp.path, fileext = '.txt')
  write.table(matrix, file = tmpMatrix.path, quote = F, 
              col.names = F, row.names = F, sep = sep) 
  return(tmpMatrix.path)
}

#' Title
#'
#' @param tmpMatrix.path
#' @param k.min  
#' @param k.max 
#' @param outer.iter 
#' @param inner.iter  
#' @param conver.test.niter  
#' @param conver.test.stop.threshold 
#' @param out.dir 
#'
#' @return
#' @export
#'
#' @examples
runNmfGpu <- function(tmpMatrix.path, k.min= 2, k.max = 2, outer.iter = 10,
                      inner.iter = 10^4, conver.test.niter = 10, 
                      conver.test.stop.threshold = 40, out.dir = NULL) {
  
  # Define pattern to finde GPU_NMF output.
  tmp.dir <- dirname(tmpMatrix.path)
  h.pattern <- sprintf('%s_H.txt', basename(tmpMatrix.path))
  w.pattern <- sprintf('%s_W.txt', basename(tmpMatrix.path))
  
  # Check if deconv. matrix should be saved
  if(!is.null(out.dir)) dir.create(out.dir)
  
  # RUN NMF.
  dec.matrix <- lapply(k.min:k.max, function(k) {
    print(Sys.time())
    cat("Factorization rank: ", k, "\n")
    k.matrix <- lapply(1:outer.iter, function(i) {
      if(i%%10 == 0) { cat("\tIteration: ", i, "\n") }
      frob.error <- 1
      while(frob.error == 1 | is.na(frob.error)){
        # VERSION 1.0        
        # nmf.cmd <- sprintf('NMF_GPU %s -k %s -i %s', tmpMatrix.path, k, inner.iter)
        # nmf.stdout <- system(nmf.cmd, intern = T)
        nmf.cmd <- sprintf('%s -k %i -i %i -j %i -t %i', tmpMatrix.path, k,
                           inner.iter, conver.test.niter,
                           conver.test.stop.threshold)
        nmf.stdout <- system2('NMF_GPU',args = nmf.cmd, stdout = T, stderr = NULL)
        frob.error <- nmf.stdout[grep(nmf.stdout, pattern = 'Distance')]
        frob.error <- as.numeric(sub(".*: ", "", frob.error))      
        if(length(frob.error)==0) frob.error <- 1
      }
      h.file <- list.files(tmp.dir, pattern = h.pattern, full.names = T)
      h.matrix <- fread(h.file, header = F)
      w.file <- list.files(tmp.dir, pattern = w.pattern, full.names = T)
      w.matrix <- fread(w.file, header = F)
      if(!is.null(out.dir)) {
        file.copy(h.file, file.path(out.dir, sprintf('H_k%s_iter%s.txt', k, i)))     
        file.copy(w.file, file.path(out.dir, sprintf('W_k%s_iter%s.txt', k, i)))
      } 
      file.remove(c(h.file, w.file))
      return(list(H = h.matrix,
                  W = w.matrix,
                  Frob.error = frob.error))
    })
    names(k.matrix) <- 1:outer.iter
    return(k.matrix)
  })
  names(dec.matrix) <- k.min:k.max
  return(dec.matrix)
}

#==============================================================================#
#                               Getter FUNCTIONS                               #
#==============================================================================#
#' Title
#'
#' @param dec.matrix 
#'
#' @return
#' @export
#'
#' @examples
getFrobError <- function(dec.matrix) {
  frob.errorList <- lapply(dec.matrix, function(dec.m) lapply(dec.m, function(m) m$Frob.error))
  frob.errorMatrix <- melt(frob.errorList)
  frob.errorMatrix <- dcast(frob.errorMatrix, L2 ~ L1, value.var = 'value')
  frob.errorMatrix <- frob.errorMatrix[order(as.numeric(frob.errorMatrix$L2)),]
  frob.errorMatrix <- frob.errorMatrix[,-1]
  frob.errorMatrix <- frob.errorMatrix[,order(as.numeric(colnames(frob.errorMatrix)))]
  #colnames(frob.errorMatrix) <- 2:(ncol(frob.errorMatrix)+1)
  frob.errorMatrix[frob.errorMatrix == 1] <- NA
  return(frob.errorMatrix)
}

#' Title
#'
#' @param dec.matrix 
#'
#' @return
#' @export
#'
#' @examples
getHmatrixList <- function(dec.matrix) {
  frob.errorMatrix <- getFrobError(dec.matrix)
  frob.errorMatrix[is.na(frob.errorMatrix)] <- NA
  k.max <- ncol(frob.errorMatrix) + 1
  ## Choose optimal factorization rank by hand from plots.
  H.list <- lapply(2:k.max, function(k.opt) {
    i.minFrobError <- apply(frob.errorMatrix, 2, function(x) which(x == min(x, na.rm = TRUE))) 
    i.minFrobError <- i.minFrobError[as.character(k.opt)]
    i.minFrobError <- unlist(i.minFrobError)[1]
    H.opt <- as.data.frame(dec.matrix[[as.character(k.opt)]][[i.minFrobError]]$H)
    rownames(H.opt) <- paste('Faktor', 1:k.opt, sep = '') 
    return(H.opt)
  })
  names(H.list) <- 2:k.max
  return(H.list)
}

#' Title
#'
#' @param dec.matrix 
#'
#' @return
#' @export
#'
#' @examples
getWmatrixList <- function(dec.matrix) {
  frob.errorMatrix <- getFrobError(dec.matrix)
  W.list <- lapply(2:k.max, function(k.opt) {
    i.minFrobError <- apply(frob.errorMatrix, 2, function(x) which.min(x))
    i.minFrobError <- i.minFrobError[as.character(k.opt)]
    W.opt <- as.data.frame(dec.matrix[[as.character(k.opt)]][[i.minFrobError]]$W)
    colnames(W.opt) <- paste('Faktor', 1:k.opt, sep = '') 
    return(W.opt)
  })
  names(W.list) <- 2:k.max
  return(W.list)
}

#==============================================================================#
#             Criteria for optimal factorization rank - FUNCTIONS              #
#==============================================================================#
#' Title
#'
#' @param frob.errorMatrix 
#'
#' @return
#' @export
#'
#' @examples
computeFrobErrorStats <- function(frob.errorMatrix) {
  min.frobError <- apply(frob.errorMatrix, 2, function(x) min(x, na.rm = T))
  sd.frobError <- apply(frob.errorMatrix, 2, function(x) sd(x, na.rm = T))
  mean.frobError <- colMeans(frob.errorMatrix, na.rm = T)
  cv.frobError <- sd.frobError/mean.frobError
  frobError.data <- data.frame('k' = as.numeric(names(min.frobError)),
                               'min' = min.frobError, 
                               'mean' = mean.frobError,
                               'sd' = sd.frobError,
                               'cv' = cv.frobError)
  return(frobError.data)
}

# Compute p-value with t-test for different factorisation ranks.
#' Title
#'
#' @param frobError.matrix 
#'
#' @return
#' @export
#'
#' @examples
computePValue4FrobError <- function(frobError.matrix) {
  n.col <- ncol(frobError.matrix)
  p <- mapply(1:(n.col-1), 2:n.col, FUN = function(i,j) {
    p <- t.test(frobError.matrix[,i], frobError.matrix[,j])
    return(p$p.value)
  })
  #p <- p.adjust(unlist(p), method = 'BH') 
  return(-log10(unlist(p)))
}

# Cosine similarity
#' Title
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export
#'
#' @examples
cosineSim <- function(a,b){
  return(sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) )) 
}

# Cosine distance
#' Title
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export
#'
#' @examples
cosineDist <- function(a,b){
  return(1 - cosineSim(a,b)) 
}

# Create distance matrix with cosine similarity
#' Title
#'
#' @param in.matrix 
#' @param in.dimension 
#'
#' @return
#' @export
#'
#' @examples
cosineDiss <- function(in.matrix, in.dimension=2){
  if(in.dimension == 1) in.matrix <- t(in.matrix)
  cosineDist.list <- lapply(1:ncol(in.matrix), function(i.outCol) {
    cosine.dists <- unlist(lapply(1:ncol(in.matrix), function(i.innCol) {
      cosineDist(in.matrix[,i.outCol], in.matrix[,i.innCol])
    }))
  })
  diss.matrix <- do.call(rbind, cosineDist.list)  
  return(round(diss.matrix, digits=14))
}

# Create distance matrix with cosine similarity with matrix operations
#' Title
#'
#' @param in.matrix 
#' @param in.dimension 
#'
#' @return
#' @export
#'
#' @examples
cosineDissMat <- function(in.matrix, in.dimension=2){
  if(in.dimension == 1) in.matrix <- t(in.matrix)
  squaredVectorSum <- apply(in.matrix, 2, function(m) {sqrt(sum(m*m))})
  squaredVectorProduct <- squaredVectorSum %*% t(squaredVectorSum)
  squaredInputSum <- t(in.matrix) %*% in.matrix # sum(a*b) for any a,b in M
  diss.matrix <- 1 - squaredInputSum / squaredVectorProduct  # CosineDistance = 1 - CosineSimilarity
  return(round(diss.matrix, digits=14))
}

# Compute Alexandrov Criterion --> Silhoutte Width.
#' Title
#'
#' @param dec.matrix 
#'
#' @return
#' @export
#'
#' @examples
computeSilhoutteWidth <- function(dec.matrix) {
  sil.vec <- lapply(dec.matrix, function(m) {
    concat.matrix <- do.call(cbind, lapply(m, function(x) x$W))
    dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
    my.pam <- pam(dist.matrix, k = ncol(m[[1]]$W),  diss = T)
    sil.sum <- sum(my.pam$silinfo$widths[,'sil_width'])
    sil.mean <- mean(my.pam$silinfo$widths[,'sil_width'])
    return(list(sum = sil.sum,
                mean = sil.mean))
  })
  return(sil.vec)  
}

# Compute Cophenetic correlation coefficient, TO BE IMPROVED
#' Title
#'
#' @param dec.matrix 
#'
#' @return
#' @export
#'
#' @examples
computeCopheneticCoeff <- function(dec.matrix) {
  coph.coeff <-  lapply(dec.matrix, function(m) {
    concat.matrix <- do.call(cbind, lapply(m, function(x) x$W))
    dist.matrix <- cosineDissMat(as.matrix(concat.matrix))
    my.hclust <- hclust(as.dist(dist.matrix))
    dist.cophenetic <- as.matrix(cophenetic(my.hclust))
    
    # take distance matrices without diagonal elements
    diag(dist.matrix) = NA
    dist.matrix<-dist.matrix[which(!is.na(dist.matrix))]
    diag(dist.cophenetic) = NA
    dist.cophenetic<-dist.cophenetic[which(!is.na(dist.cophenetic))]
    
    return(cor(cbind(dist.cophenetic, dist.matrix))[1,2])
  })
}

#' Compute amari type distance between two matrices
#'
#' @param matrix.A,matrix.B of the same dimensionality
#'
#' @return The amari type distance of matrix.A & matrix.B according
#'        to [Wu et. al, PNAS 2016]
#' @export
#'
#' @examples
amariDistance <- function(matrix.A, matrix.B) {
  K <- dim(matrix.A)[2]
  C <- cor(matrix.A, matrix.B)
  return(1 - (sum(apply(C, 1, max)) + sum(apply(C,2,max))) / (2 * K))
}

#' Compute Amari Distances from [Wu et. al, PNAS 2016]
#'
#' @param dec.matrix list of NMF results for different k
#'
#' @return The average Amari-type error for each k
#' @export
#'
#' @examples
computeAmariDistances <- function(dec.matrix){
  distance.averages <- lapply(dec.matrix, function(m) {
    matrices <- lapply(m, function(x) as.matrix(x$W))
    B <- length(matrices)
    distances.list <- unlist(lapply(1:(B-1), function(b) {
      distances <- lapply((b+1):B, function(b.hat) {
        amariDistance(matrices[[b]], matrices[[b.hat]])
      })
    }))
    return(sum(distances.list) / (B*(B-1)/2))
  })
  return(distance.averages)
}

#==============================================================================#
#                         H-MATRIX ANALYSIS FUNCTIONS                          #
#==============================================================================#
# Compute unbaised signature names by applying a row k-means
# and classify signature according to cluster and highest cluster mean.
#' Title
#'
#' @param H 
#' @param meta.data 
#'
#' @return
#' @export
#'
#' @examples
getSignatureNames <- function(H, meta.data) {
  sig.names <- lapply(1:nrow(H), function(i) {
    k.cluster <- kmeans(as.numeric(H[i,]), 2)
    n.kCluster <- which(k.cluster$centers == max(k.cluster$centers))
    samples.kCluster <- colnames(H)[k.cluster$cluster == n.kCluster]
    m <- meta.data[meta.data$ID%in%samples.kCluster,]
    sig.name <- paste(sort(unique(m$Timepoint)), collapse = '/')
    sig.name <- paste(unique(m$Stage), sig.name, sep = '_')
    return(sig.name)
  })
  return(unlist(sig.names))
}

#==============================================================================#
#                         W-MATRIX ANALYSIS FUNCTIONS                          #
#==============================================================================#
# Compute 'shannon' entropy per region. 
# High Entropy means highly specific for one signature.
#' Title
#'
#' @param W 
#'
#' @return
#' @export
#'
#' @examples
computeEntropy <- function(W) {
  W.relativ <- t(apply(W, 1, function(x) x/sum(x)))
  W.entropy <- apply(W.relativ, 1, function(x) {
    p <- x*log2(length(x)*x)
    p[is.na(p)] <- 0
    h <- sum(p)
    return(h)
  })
  return(W.entropy)
}

# Compute delta between each signature per row
#' Title
#'
#' @param W 
#'
#' @return
#' @export
#'
#' @examples
computeSigAbsDelta <- function(W) {
  delta.regions <- lapply(1:ncol(W), function(k) {
    delta.vec <- W[,k]-(rowSums(W[,-k]))
    return(delta.vec)
  })
  delta.regions <- do.call(cbind, delta.regions)
  return(delta.regions) 
}

#' Title
#'
#' @param W 
#'
#' @return
#' @export
#'
#' @examples
computeCoefVar <- function(W) {
  apply(W, 1, function(r) sd(r)/mean(r))
}

# Perform Kmeans on W rows to extract all possible signature combinations.
#' Title
#'
#' @param W 
#'
#' @return
#' @export
#'
#' @examples
performRowKmeans <- function(W) {
  k.row <- apply(W, 1, function(x) {
    x.trans <- sigmoidTransform(x)
    k.cluster <- kmeans(x.trans, 2)
    d <- dist(as.data.frame(x.trans))
    sil.mean <- mean(silhouette(k.cluster$cluster, d)[,3])
    return(list(centers = t(k.cluster$centers),
                silmean = sil.mean,
                explainedVar = k.cluster$betweenss/k.cluster$totss,
                oddsVar = sum(k.cluster$withinss)/k.cluster$betweenss,
                attribution = k.cluster$cluster)) 
  })
  return(k.row)
}

#' Title
#'
#' @param W 
#' @param var.thres 
#'
#' @return
#' @export
#'
#' @examples
computeSignatureCombs <- function(W, var.thres = 0.25) {
  # Determine Region contributions
  k.row <- performRowKmeans(W)
  # Extract Signature combinations generated by k-means
  k.attribution <- lapply(k.row, function(r) abs(r$attribution))
  k.attribution <- do.call(rbind, k.attribution)
  k.ids <- apply(k.attribution, 1, function(r) paste(r, collapse = ''))
  # Determine regions which contribute to all signatures.
  k.varCoef <- computeCoefVar(W)
  all.signature <- which(k.varCoef < var.thres)
  k.ids[all.signature] <- gsub('2', '1', k.ids[all.signature])
  return(k.ids)
}

#' Title
#'
#' @param k.ids 
#'
#' @return
#' @export
#'
#' @examples
getIndex4SignatureCombs <- function(k.ids) {
  # Get for all signature combinations corresponding regions.
  ids <- sort(unique(k.ids))[-1]
  n.ids <- length(ids)+1
  i.regions <- lapply(1:(n.ids/2), function(i) {
    sub.id <- ids[c(i, n.ids - i)]
    i.region <- which(k.ids%in%sub.id)
    return(i.region)
  })
  names(i.regions) <- ids[1:(n.ids/2)]
  # Get regions for all signature id.
  allSig.id <- sort(unique(k.ids))[1]
  i.regions[[allSig.id]] <- which(k.ids%in%allSig.id)
  return(i.regions)
}

#' Title
#'
#' @param W 
#' @param i.regions 
#'
#' @return
#' @export
#'
#' @examples
computeMeanDiff4SignatureCombs <- function(W, i.regions) {
  # Compute Row mean diff between Cluster 1 and 2.
  w.diffs <- lapply(1:length(i.regions), function(i) {
    w <- W[i.regions[[i]],]
    if (!grepl(names(i.regions)[i], pattern = '2')) {
      w.diff <- rowMean(w)
    } else {
      signature.comb <- as.numeric(unlist(strsplit(names(i.regions)[[i]], split = '')))
      w.mean1 <- rowMean(w[,which(signature.comb == 1)])
      w.mean2 <- rowMean(w[,which(signature.comb == 2)])
      w.diff <- w.mean1 - w.mean2
    }
    return(w.diff)
  })
  names(w.diffs) <- names(i.regions)
}

### Get number of peaks per signature combination and create a barplot.
#' Title
#'
#' @param i.regions 
#' @param w.diffs 
#'
#' @return
#' @export
#'
#' @examples
getSignatureCombCounts <- function(i.regions, w.diffs) {
  n.peaks <- lapply(names(i.regions), function(region.n) {
    n.peak1 <- sum(w.diffs[[region.n]] > 0)
    n.peak2 <- sum(w.diffs[[region.n]] < 0)
    n.peak <- c(n.peak1, n.peak2)
    n.peak <- melt(n.peak)
    n.peak$clusterId <- c('1', '2')
    n.peak$sigCombId <- region.n
    return(n.peak) 
  })
  n.peaks <- do.call(rbind, n.peaks)
  return(n.peaks)
}

