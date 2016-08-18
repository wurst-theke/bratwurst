#==============================================================================#
# Author: Sebastian Steinhauser - s.steinhauser@gmx.net
# Date: 20.07.2016
# Comments: Bratwurst test script to run NMF-GPU on published ALL/AML RNA-seq data .
#==============================================================================#
# Load Libraries
library(Bratwurst)

#==============================================================================#
#                                   MAIN                                       #
#==============================================================================#
result.path <- '/home/thymin'
samples <- 'AML'
# Set datapath and find files to analyse.
data.path  <- file.path(getwd(), 'data')
matrix.file <- list.files(data.path, 'data.txt', full.names = T)
rowAnno.file <- list.files(data.path, 'micro.*anno.*txt', full.names = T)
rowAnno.bed <- list.files(data.path, '.bed', full.names = T)
colAnno.file <- list.files(data.path, 'sample.*anno.*txt', full.names = T)

# Read files to summarizedExperiment
nmf.exp <- nmfExperimentFromFile(matrix.file = matrix.file,
                                 rowAnno.file = rowAnno.file, 
                                 colData.file = colAnno.file)

### RUN NMF GPU
# RUN NMF.
k.max <- 4
outer.iter <- 10
inner.iter <- 10^4

nmf.exp<- runNmfGpu(nmf.exp = nmf.exp,
                    k.max = k.max,
                    outer.iter = outer.iter,
                    inner.iter = inner.iter)

## Check Getter Functions.
HMatrixList(nmf.exp, k = 2)
WMatrixList(nmf.exp)
FrobError(nmf.exp)
HMatrix(nmf.exp, k = 2)
WMatrix(nmf.exp)

#==============================================================================#
#                   Criteria for optimal factorization rank                    #
#==============================================================================#
# Get frob error from list of deconv. matrix and compute error statistics.
nmf.exp <- computeFrobErrorStats(nmf.exp)

# Generate Alexandrov Criterion plot
nmf.exp <- computeSilhoutteWidth(nmf.exp)

# Cophenetic correlation coefficient plot
nmf.exp <- computeCopheneticCoeff(nmf.exp)

# Compute amari type distance
nmf.exp <- computeAmariDistances(nmf.exp)

OptKStats(nmf.exp)

### Generate plots to estimate optimal K
gg.optK <- plotKStats(nmf.exp)

# Save plots for optimal k estimation.
optK.svg <- sprintf('%s_k%s_iter%s_optKPlot.%s', samples, k.max,
                    outer.iter, c('svg', 'png'))
sapply(1:length(optK.svg), function(i) {
  save_plot(gg.optK, filename = file.path(result.path, optK.svg[i]), base_aspect_ratio = 1.4)
})

# Generate ranked frob errors plot & save as svg/png.
gg.rankedFrobError <- plotRankedFrobErrors(nmf.exp)

rankedFrobError.svg <- sprintf('%s_k%s_iter%s_rankedFrobError.%s', 
                               samples, k.max, outer.iter, c('svg', 'png'))
sapply(1:length(rankedFrobError.svg), function(i) {
  save_plot(plot = gg.rankedFrobError, 
            filename = file.path(result.path, rankedFrobError.svg[i]))
})

#==============================================================================#
#                        H-MATRIX ANALYSIS/VISUALIZATION                       #
#==============================================================================#
# Plot Heatmaps for H over all k
# HMatrixHeat.svg <- sprintf('%s_k%s_iter%s_HmatrixHeatmap.svg', 
#                            samples, k.max, outer.iter) 
# svg(file.path(result.path, HMatrixHeat.svg), width = 15)
plotHeatmap4MatrixList(H.list, trans = T)
# dev.off()

heat.anno <- HeatmapAnnotation(df = meta.data[,c(1,3,4)],
                               col = list(Cell = cell.col, 
                                          Timepoint = time.col,
                                          Stage = stage.col))

# Create ouput dir structure.
hmatrix.resultPath <- file.path(result.path, 'H-matrix')
dir.create(hmatrix.resultPath)
sapply(c('png', 'svg'), function(suffix) {
    dir.create(file.path(hmatrix.resultPath, suffix))
})

# Print H-Matrix for all K's in png & svg.
sapply(1:length(H.list), function(i) {
  h.heatmap <- Heatmap(H.list[[i]], col = col.map, 
                       clustering_distance_columns = 'pearson',
                       heatmap_legend_param = list(color_bar = "continuous"),
                       show_column_names = F, cluster_rows = F,
                       top_annotation = heat.anno)
  png(file.path(hmatrix.resultPath, 'png', 
                sprintf('heart_Hmatrix_%sk.png', i+1)))
  draw(h.heatmap)
  dev.off() 
  svg(file.path(hmatrix.resultPath, 'svg', 
                sprintf('heart_Hmatrix_%sk.svg', i+1)))
  draw(h.heatmap)
  dev.off() 
})
#save(H.list, W.list, anno.peaks, regionDB, file = '/home/thymin/NMF-GPU_toolbox/heart_nmf.RData') #file.path(result.path, 'heart_nmf.RData'))

#==============================================================================#
#                        W-MATRIX ANALYSIS/VISUALIZATION                       #
#                       FIND REGIONS ASSOCIATED TO SIGNATURES                  #
#==============================================================================#
### Find representative regions.
# Get W for best K
k.opt <- 4
signature.names <- getSignatureNames(nmf.exp, k.opt)

# Compute 'Delta Signature contribution' and select regions for signature of interest.
W.delta <- computeSigAbsDelta(W)
# Compute Entropy per row
W.entropy <- computeEntropy(W)

# Perform Kmeans on W rows to extract all possible signature combinations.
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

# Compute different variance vars. to determine 11111 signature.
k.explainedVar <- unlist(lapply(k.row, function(r) r$explainedVar))
k.oddsVar <- unlist(lapply(k.row, function(r) r$oddsVar))
k.varCoef <- apply(W, 1, function(r) sd(r)/mean(r))

var.thres <- 0.25
all.signature <- which(k.varCoef < var.thres)

# Extract Signature combinations generated by k-means
k.attribution <- lapply(k.row, function(r) abs(r$attribution))
k.attribution <- do.call(rbind, k.attribution)
k.ids <- apply(k.attribution, 1, function(r) paste(r, collapse = ''))
k.ids[all.signature] <- gsub('2', '1', k.ids[all.signature])

# Extract Centroids for each row.
k.scores <- lapply(k.row, function(r) r$centers)
k.scores <- do.call(rbind, k.scores)

# Extract mean silhouette width for each row.
k.silmean <- lapply(k.row, function(r) r$silmean)
k.silmean <- unlist(k.silmean)

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

# Compute regions scores (delta centroid values) to distinguish between
# different signature combinations within a 'row cluster'
# region.scores <- lapply(1:(n.ids/2), function(i) {
#   k.score1 <- k.scores[which(k.ids == ids[i]),]
#   k.score2 <- k.scores[which(k.ids == ids[n.ids-i]), 2:1]
#   kscore.df <- rbind(k.score1, k.score2)
#   k.sil <- c(k.silmean[which(k.ids == ids[i])], 
#              k.silmean[which(k.ids == ids[n.ids-i])])
#   kscore.df <- cbind(kscore.df, k.sil, as.numeric(ids[i]))
#   return(kscore.df)
# })
# region.scores <- do.call(rbind, region.scores)
# colnames(region.scores) <- c('x', 'y', 'SilWidth', 'Sig')
# region.scores <- as.data.frame(region.scores)

# MA like score plotting useful?
# region.scores$d <- region.scores[,1] - region.scores[,2]
# region.scores$M <- log2(region.scores[,1]/region.scores[,2])
# region.scores$A <- 1/2*(region.scores[,1] + region.scores[,2])

# Compute Row mean diff between Cluster 1 and 2.
rowMean <- function(m) {apply(as.matrix(m), 1, mean)}

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

# gg.scatter <- ggplot(as.data.frame(region.scores), aes(x = x, y = y)) + geom_point()
# gg.scatter <- gg.scatter + stat_function(fun = function(x) {x}, col = 'red',
#                                          linetype = 'dashed', size = 1.25)
# gg.scatter <- gg.scatter + facet_wrap(~Sig, scales = 'free')
# gg.scatter <- gg.scatter + xlab('1') + ylab('2')
# gg.scatter <- gg.scatter + theme_bw() + science_theme 
# gg.scatter + theme(strip.background = element_rect(fill = 'white'))

# Plot W for given signature combination.
w <- W[i.regions[[1]],]
w <- w[region.scores[i.regions[[1]], 'd'] < 0,]
#w <- w[order(w[,5], decreasing = T),]
png(file.path(result.path, sprintf('heaeRrt_Wmatrix_%sk.png', 5)))
Heatmap(w, cluster_rows = F, cluster_columns = F, show_row_names = F,
        heatmap_legend_param = list(color_bar = "continuous")) 
dev.off() 

# Extract regions for a given signature combination
j <- 6# length(i.regions) - 1
p <- anno.peaks[i.regions[[j]],]
p[order(w.diffs[[j]], decreasing = T),]

#==============================================================================#
#      Analysis of Signature combination for regulatory element enrichment.    #
#==============================================================================#
### Get number of peaks per signature combination and create a barplot.
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

gg.bar <- ggplot(n.peaks, aes(x = sigCombId, y = value, fill = clusterId)) 
gg.bar <- gg.bar + geom_bar(colour='black', stat='identity', 
                            position=position_dodge())
gg.bar <- gg.bar + xlab('Signature combinations') + ylab('#Regions')
gg.bar <- gg.bar + scale_fill_manual(values = c('red', 'blue'),
                                     name = 'Cluster')
gg.bar <- gg.bar + theme_bw() + science_theme 
gg.bar <- gg.bar + theme(axis.text.x = element_text(angle = 90, hjust = 1))
save_plot(gg.bar, filename = file.path(result.path, 'clusterCombination_bar.svg'))

#### Promoter/Enhancer/Superenhancer enrichment in different signature combinations
## FUNCTIONS.
spliteRegionsByScore <- function(regions) {
  up.regions <- regions[mcols(regions)[,1] > 0,]
  up.regions <- up.regions[order(mcols(up.regions)[,1], decreasing = T),]
  down.regions <- regions[mcols(regions)[,1] < 0,]
  down.regions <- down.regions[order(mcols(down.regions)[,1]),]
  return(list('1' = up.regions, '2' = down.regions))
}

computeRegionFC <- function(peaks, subregions, universe) {
  prop.universe <- sum(countOverlaps(universe, peaks))/length(universe)
  prop.subregions <- sum(countOverlaps(subregions, peaks))/length(subregions)
  fc <- prop.subregions/prop.universe
  return(fc)
}

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
computeRecoveryAUC <- function(ranks, n) {
    ranks <- unique(c(1, sort(ranks[ranks < n]), n))
    rank.diff <- ranks[-1] - ranks[-length(ranks)] 
    #recov <- 0:c(length(rank.diff)-1)
    recov <- c(0, cumsum(rep(1/n, length(rank.diff)-1)))
    auc <- sum(recov*rank.diff)/n 
    return(auc)
}

computeRmdAUCs <- function(ranks, n, iter = 10^3, threads = 3) {
  n.hits <- sum(ranks <= n)
  rmd.aucs <- mclapply(1:iter, function(i) {
    rmd.ranks <- sort(sample(1:n, size = n.hits, replace = F))
    rmd.auc <- computeRecoveryAUC(rmd.ranks, n)
    return(rmd.auc)
  }, mc.cores = threads)
  return(unlist(rmd.aucs))
}

computeNES <- function(ranks, n, iter = 10^3, threads = 3) {
  auc <- computeRecoveryAUC(ranks, n)
  rmd.aucs <- computeRmdAUCs(ranks, n, iter = iter, threads = threads)
  nes <- (auc - mean(rmd.aucs))/sd(rmd.aucs)
  r <- sum(auc < rmd.aucs) + 1
  result <- data.frame('Rank' = r, 'NES' = nes, 'AUC' = auc)
  return(result)
}

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

## Analysis
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#gencode.gff <- '/home/thymin/genomes/mm10/'
promoter.regions <- promoters(TxDb.Mmusculus.UCSC.mm10.ensGene, 
                              upstream = 1500, downstream = 500)
universe <- GRanges(anno.peaks$V2, IRanges(anno.peaks$V3, anno.peaks$V4))

# Compute promoter overlap/enrichment.
promoter.enrichment <- mclapply(1:length(i.regions), function(i) {
  subregions <- universe[i.regions[[i]],]
  subregions$score <- w.diffs[[i]]
  subregions <- spliteRegionsByScore(subregions)
  peak.stats <- computeRunningPeakStats4Regions(peaks = promoter.regions, 
                                                subregions = subregions,
                                                universe =  universe)
  peak.stats$SignatureComb <- names(i.regions)[i]
  return(peak.stats)
}, mc.cores = 1)
promoter.enrichment <- do.call(rbind, promoter.enrichment)
promoter.enrichment <- promoter.enrichment[order(promoter.enrichment$ClusterID, promoter.enrichment$SignatureComb),]
promoter.enrichment <- promoter.enrichment[!is.na(promoter.enrichment$FC),]


promoter.data <- promoter.enrichment[promoter.enrichment$RankProp == 1,]
promoter.data$ClusterID <- factor(promoter.data$ClusterID)
promoter.data[order(promoter.data$NES)]

gg.enrichmentBar <- ggplot(promoter.data, aes(x = SignatureComb, y = FC, fill = ClusterID)) 
gg.enrichmentBar <- gg.enrichmentBar + geom_bar(colour='black', stat='identity', 
                                                position=position_dodge())
gg.enrichmentBar <- gg.enrichmentBar + scale_fill_manual(values = c('red', 'blue'),
                                                         name = 'Cluster')
gg.enrichmentBar <- gg.enrichmentBar+ theme_bw() + science_theme 
gg.enrichmentBar <- gg.enrichmentBar + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Compute superenhancer overlap/enrichment.
superenhancer.bed <- '/home/thymin/Projects/bloodEnhancer_project/heart_superenhancers.bed'
superenhancers <- import.bed(superenhancer.bed)

# Compute promoter overlap/enrichment.
superenhancer.enrichment <- mclapply(1:length(i.regions), function(i) {
  print(i)
  subregions <- universe[i.regions[[i]],]
  subregions$score <- w.diffs[[i]]
  subregions <- spliteRegionsByScore(subregions)
  peak.stats <- computeRunningPeakStats4Regions(peaks = superenhancers, 
                                                subregions = subregions,
                                                universe =  universe)
  peak.stats$SignatureComb <- names(i.regions)[i]
  peak.stats <- peak.stats[!is.na(peak.stats$FC),]
  return(peak.stats)
}, mc.cores = 1)
superenhancer.enrichment <- do.call(rbind, superenhancer.enrichment)

p <- promoter.enrichment[promoter.enrichment$RankProp == 1,]
p$anno <-
plot.matrix <- 
#barplot(p$FC)
signature.comb <- p[,c('ClusterID', 'SignatureComb')]
sortSignatureCombs <- function(signature.comb) {
  
}

# Compute enrichment of mouse encode chromHMM model.
chromstat.bed <- '/home/thymin/mouseEncode_chromHMM/heart_chromatinStates_mm10.bed'
chrom.stats <- import.bed(chromstat.bed)

chromStat.enrichment <- mclapply(1:length(i.regions), function(i) {
  subregions <- universe[i.regions[[i]],]
  subregions$score <- w.diffs[[i]]
  subregions <- spliteRegionsByScore(subregions)
  chromStat.peakStats <- lapply(sort(unique(chrom.stats$name)), function(n) {
    chromStat.n <- chrom.stats[chrom.stats$name == n,]
    peak.stats <- computePeakStats4Regions(peaks = chromStat.n, 
                                           subregions = subregions,
                                           universe =  universe)
    peak.stats$chromStat <- n 
    return(peak.stats)
  })
  chromStat.peakStats <- do.call(rbind, chromStat.peakStats)
  return(chromStat.peakStats)
}, mc.cores = 4)


##### 
# Compute promoter overlap/enrichment.
chip.peaks <- regionDB$regionGRL 


  fc <- computeRegionFC4GRlist(peak = chip.peaks, subregion = subregions[[1]], universe = promoter.universe) 


library('LOLA')
promoter.universe <- subsetByOverlaps(universe, promoter.regions)

promoter.chipEnrichment <- lapply(1:length(i.regions), function(i) {
  print(i)
  subregions <- universe[i.regions[[i]],]
  subregions$score <- w.diffs[[i]]
  subregions <- subsetByOverlaps(subregions, promoter.regions)
  subregions <- spliteRegionsByScore(subregions)
  lola.results <- lapply(1:length(subregions), function(i.subregion){
    lola.result <- runLOLA(userSets = subregions[[i.subregion]], 
                           userUniverse = promoter.universe,
                           regionDB, cores = 5)
    lola.result$adjPLog <- -log10(p.adjust(p = 10^(-lola.result$pValueLog),
                                           method = 'fdr'))
    lola.result$clusterId <- i.subregion
    return(lola.result)
  })
  lola.results <- do.call(rbind, lola.results)
  lola.results$SigCombId <- names(i.regions)[i]
  return(lola.results)
})

library('rGREAT')
test <- submitGreatJob(gr = subregions[[1]], species = 'mm10', bg = promoter.universe, request_interval = 30)  
availableCategories(test)
tb <- getEnrichmentTables(test, category = "Gene Expression")
==============================================================================#
#                         LOLA [PMID:26508757]-                                #
#       Perform enrichment analysis with signature associated regions          #
#                   with FISHER's EXACT TEST (--> Peak OV)                     #
#==============================================================================#

dbPath <- '/home/steinhau/NMF-GPU_toolbox/LOLACore/mm10/' 
dbPath <- '/home/thymin/Projects//NMF-GPU_toolbox/LOLACore/mm10/' 
regionDB <- loadRegionDB(dbPath, useCache = T)

# Define the universe
universe <- GRanges(anno.peaks$V2, 
                    IRanges(anno.peaks$V3, anno.peaks$V4))

lolaEnrichment.results <- lapply(names(i.regions), function(n.region) {
  print(n.region)
  regions <- universe[i.regions[[n.region]]]
  cluster.1 <- which(w.diffs[[n.region]] > 0)
  cluster.2 <- which(w.diffs[[n.region]] < 0)
  lola.results <- lapply(list(cluster.1, cluster.2), function(j.cluster) {
    lola.result <- runLOLA(regions[j.cluster,], universe, regionDB, cores=6)
    #lola.result$fdr <- p.adjust(10^-lola.result$pValueLog, method = 'fdr')
    return(lola.result)
  })
  return(lola.results)
}) 


### Prepare feature matrix for plotting.
x <- lapply(lolaEnrichment.results, function(lolaEnrichment.result) { 
  lapply(lolaEnrichment.result, function(l) {
    l$fdr <- p.adjust(10^-l$pValueLog, method = 'fdr')
    return(l)
  })
})
names(x) <- names(i.regions)

fdr.thres <- 0.01
term.dict <- lapply(x, function(sublist) {
  lapply(sublist, function(l) {
    l <- l[l$fdr <= fdr.thres,]
    l[1:5,]$filename
  })
})
term.dict <- unique(unlist(term.dict))
term.dict <- term.dict[!is.na(term.dict)]

x <- lapply(names(i.regions), function(n.regions) {
  n.cluster <- 1
  r <- lapply(x[[n.regions]], function(r) {
    r$clusterID <- n.cluster
    n.cluster <<- n.cluster + 1
    return(r)
  })
  r <- do.call(rbind, r)
  r$signatureID <- n.regions
  return(r)
})
x <- do.call(rbind, x)

filtered.x <- x[x$filename%in%term.dict,]
filtered.x$id <- paste(filtered.x$signatureID, filtered.x$clusterID)
signature.comb <- do.call(rbind, strsplit(filtered.x$signatureID, ''))
enrichment.matrix <- dcast(filtered.x, formula = filename ~ id, value.var = 'logOddsRatio') #fdr")
rownames(enrichment.matrix) <- enrichment.matrix$filename
enrichment.matrix <- enrichment.matrix[,-1]
#enrichment.matrix <- -log10(enrichment.matrix + 10^-120)

### Get row annotation data.
rowAnno.df <- unique(filtered.x[,c('cellType', 'antibody', 'filename'),with=F])
i.rmv <- grep(rowAnno.df$filename, pattern = 'Input')
rowAnno.df <- rowAnno.df[-i.rmv,]
enrichment.matrix <- enrichment.matrix[-i.rmv,]
i.na <- which(is.na(rowAnno.df$cellType))
rowAnno.df[i.na]$antibody <- gsub('.bed', '', rowAnno.df[i.na]$filename)

# Replace embryonic stem cell
i.stem <- grep(rowAnno.df$cellType, pattern = 'Embryonic', ignore.case = T)
rowAnno.df$cellType[i.stem] <- 'Embryonic Stem Cells'
rowAnno.df$anno <- paste(rowAnno.df$cellType,  rowAnno.df$antibody)
rowAnno.df$anno <- duplicated(rowAnno.df[,1:2,with=F]) + 1


###  Column Annotation.
sig.comb <- strsplit(colnames(enrichment.matrix), split = '')
sig.comb <- do.call(rbind, sig.comb)
sig.comb <- apply(sig.comb, 1, function(x) x[1:5] == x[7])
sig.comb <- apply(as.data.frame(sig.comb), 1, as.character)
rownames(sig.comb) <- 1:nrow(sig.comb)
colnames(sig.comb) <- paste(df = 'Factor', 1:5)

col <- c('TRUE' = 'black', 'FALSE' = 'white')
anno.colmap <- lapply(1:ncol(sig.comb), function(i) col)
names(anno.colmap) <- colnames(sig.comb)

# Col annotation order.
i.order <- apply(sig.comb, 1, function(c) sum(as.logical(c)))
i.order <- lapply(unique(sort(i.order)), function(n) {
  j <- which(i.order == n)
  j <- j[order(as.logical(sig.comb[j,5]), as.logical(sig.comb[j,4]), 
               as.logical(sig.comb[j,3]), as.logical(sig.comb[j,2]),
               decreasing = T)] 
  return(j)
})
i.order <- unlist(i.order)

# Prepare annotation heatmap.
sigComb.anno <- HeatmapAnnotation(as.data.frame(sig.comb)[i.order,],
                                  col = anno.colmap, show_legend = F,
                                  gp = gpar(col = 'black', lty = 1, lwd = 1))

Heatmap(enrichment.matrix[,i.order], cluster_rows = T, cluster_columns = F, 
        show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(color_bar = "continuous"),
        bottom_annotation = sigComb.anno)

filtered.x$fdr
table(x[x$fdr < 0.01,]$signatureID)

enriched.terms <- lapply(lola.results, function(r) r$filename[r$qValue < 0.01])

fdr.thres <- 0.1


#==============================================================================#
#                         Feature recovery analysis -                          #
#       Perform enrichment analysis with signature rankings by computing       #
#       normalized enrichment scores (NES) (--> Peak OV with x% top regions)   #
#==============================================================================#

### Compute NES
feature.peaks <- regionDB$regionGRL

j.regions  <- lapply(i.regions, function(i) {
  i[!is.promoter[i]]
})
j.regions <- i.regions

n.ranks <- 500 
featureNes.list <- lapply(1:length(i.regions), function(i) {
  print(names(i.regions)[i])
  # Get subregions and compute ranking score.
  j <- j.regions[[i]]
  subregions <- anno.peaks[j,]
  subregions.gr <- GRanges(subregions$V2, IRanges(subregions$V3, subregions$V4))
  i.keep <- names(w.diffs[[i]])%in%as.character(j)
  subregions.gr$score <- w.diffs[[i]][i.keep]
  feature.nes <- lapply(c(T, F), function(b) {
    print(b)
    if(b & sum(subregions.gr$score > 0) < n.ranks) { return(NA) }  
    if(!b & sum(subregions.gr$score < 0) < n.ranks) { return(NA) }
    # Compute feature overlap.
    sorted.regions <- subregions.gr[order(subregions.gr$score, decreasing = b)]
    ov <- findOverlaps(feature.peaks, sorted.regions)
    # Compute Normalized Enrichment Scores.
    feature.nes <- mclapply(1:length(feature.peaks), FUN = function(j) {
      ranks <- subjectHits(ov)[queryHits(ov) == j]
      ranks <- sort(unique(ranks))
      computeNES(ranks, n.ranks, threads = 1)
    }, mc.cores = 6)
    feature.nes <- unlist(feature.nes)
    return(feature.nes)
  })
  names(feature.nes) <- c('1', '2')
  return(feature.nes)
})
names(featureNes.list) <- names(i.regions)

# Extract features with nes thres.
nes.thres <- 2.5
feature.dict <- lapply(featureNes.list, function(feature.nes) {
  lapply(feature.nes, function(nes) {
    i.feature <- which(nes >= nes.thres)
    i.feature <- i.feature[!i.feature%in%is.na(nes)]
    i.feature <- i.feature[order(nes[i.feature], decreasing = T)]
    return(i.feature[1:3])
  })
})
feature.dict <- sort(unique(unlist(feature.dict)))
enriched.features <- regionDB$regionAnno[feature.dict,]
enriched.features <- enriched.features[enriched.features$collection != 'ucsc_features',]

# Create large data.table containing NES scores.
features.dt <- lapply(names(featureNes.list), function(n) {
  features <- enriched.features
  features.dt <- lapply(names(featureNes.list[[n]]), function(id) {
    features$clusterId <- id
    features$combClusterId <- n
    features$nes <- featureNes.list[[n]][[id]][feature.dict]
    return(features) 
  })
  features.dt <- do.call(rbind, features.dt)
  return(features.dt)
})
features.dt <- do.call(rbind, features.dt)
features.dt$id <- paste(features.dt$combClusterId, features.dt$clusterId)
features.dt <- features.dt[!is.na(features.dt$nes)]
features.dt <- features.dt[features.dt$collection != 'ucsc_features',]

### Compute enrichment matrix.
enrichment.matrix <- dcast(features.dt, formula = filename ~ id, value.var = 'nes')
rownames(enrichment.matrix) <- enrichment.matrix$filename
enrichment.matrix <- enrichment.matrix[,-1]

### Get row annotation data.
rowAnno.df <- unique(enriched.features[,c('cellType', 'antibody', 'filename'),with=F])
i.rmv <- grep(rowAnno.df$filename, pattern = 'Input')
rowAnno.df <- rowAnno.df[-i.rmv,]
enrichment.matrix <- enrichment.matrix[-i.rmv,]
i.na <- which(is.na(rowAnno.df$cellType))
rowAnno.df[i.na]$antibody <- gsub('.bed', '', rowAnno.df[i.na]$filename)

# Replace embryonic stem cell
i.stem <- grep(rowAnno.df$cellType, pattern = 'Embryonic', ignore.case = T)
rowAnno.df$cellType[i.stem] <- 'ESC'
rowAnno.df$anno <- paste(rowAnno.df$cellType,  rowAnno.df$antibody)
i.rep <- lapply(unique(rowAnno.df$anno), function(a) {
  j <- which(rowAnno.df$anno%in%a) 
  names(j) <- 1:length(j)
  return(j)
})
i.rep <- sort(unlist(i.rep))
rowAnno.df$anno <- names(i.rep)
rowAnno.df$anno <- paste(rowAnno.df$antibody, rowAnno.df$cellType,
                         rowAnno.df$anno)
rownames(enrichment.matrix) <- rowAnno.df$anno

###  Column Annotation.
sig.comb <- strsplit(colnames(enrichment.matrix), split = '')
sig.comb <- do.call(rbind, sig.comb)
sig.comb <- apply(sig.comb, 1, function(x) x[1:5] == x[7])
sig.comb <- apply(as.data.frame(sig.comb), 1, as.character)
rownames(sig.comb) <- 1:nrow(sig.comb)
colnames(sig.comb) <- paste(df = 'Factor', 1:5)

col <- c('TRUE' = 'black', 'FALSE' = 'white')
anno.colmap <- lapply(1:ncol(sig.comb), function(i) col)
names(anno.colmap) <- colnames(sig.comb)

# Col annotation order.
i.order <- apply(sig.comb, 1, function(c) sum(as.logical(c)))
i.order <- lapply(unique(sort(i.order)), function(n) {
  j <- which(i.order == n)
  j <- j[order(as.logical(sig.comb[j,5]), as.logical(sig.comb[j,4]), 
               as.logical(sig.comb[j,3]), as.logical(sig.comb[j,2]),
               decreasing = T)] 
  return(j)
})
i.order <- unlist(i.order)

# Prepare annotation heatmap.
sigComb.anno <- HeatmapAnnotation(as.data.frame(sig.comb)[i.order,],
                                  col = anno.colmap, show_legend = F,
                                  gp = gpar(col = 'black', lty = 1, lwd = 1))

enrichment.matrix[is.na(enrichment.matrix)] <- 0
# enrichment.matrix[enrichment.matrix > 5] <- 5

h <- Heatmap(enrichment.matrix[,i.order], cluster_rows = T, cluster_columns = T, 
             show_row_names = T, show_column_names = F,
             heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal",
                                         legend_width = unit(5, "cm"), title_position = "lefttop"),
             bottom_annotation = sigComb.anno)


png(file.path(result.path, 'heart_enhancer_lolaPeaksNES_1000.png'))
draw(h, heatmap_legend_side = "bottom")
dev.off() 
svg(file.path(result.path, 'heart_enhancer_lolaPeaksNES_1000.svg'), width = 12, height = 10)
draw(h, heatmap_legend_side = "bottom")
dev.off() 


################################################################################
motif.bed <- '/home/steinhau/NMF-GPU_toolbox/data/mENCODE/heart/H3K27ac_fimoMotifs.bed'
motifs <- fread(motif.bed, sep = '\t')
motifs.gr <- GRanges(motifs$V2, IRanges(motifs$V3, motifs$V4), 'motif' = motifs$V10)
motif.names <- unique(motifs$V10)

motif.nes <- lapply(motif.names, function(motif.n) {
  print(motif.n)
    motif.gr <- motifs.gr[motifs.gr$motif == motif.n]
    nes <- sapply(1:ncol(W.zscore), function(j) {
      j.sorted <- order(W.zscore[,j], decreasing = T)[1:n.top]
      nes <- computeNES(motif.gr, universe[j.sorted,])
      return(nes)
    })
  return(nes)
})
motif.nes <- do.call(rbind, motif.nes)
rownames(motif.nes) <- motif.names

max.nes <- apply(motif.nes, 1, function(x) max(x, na.rm = T))
motif.nes[max.nes > 3,]


library(JASPAR2014)
library(TFBSTools)
PFMatrixList <- getMatrixSet(JASPAR2014, opts = list())
pwmList <- toPWM(PFMatrixList)

peak.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, universe)

sitesetList <- searchSeq(pwmList, peak.seq, min.score="60%", strand="*")

#==============================================================================#
#          Feature enrichment analysis inspired by GREAT [20436461] -          #
#               Perform enrichment analysis with binominal test                #
#           (--> Peak Annotation OV against 'genome background')               #
#==============================================================================#
library(BSgenome.Mmusculus.UCSC.mm10)
chrs <- unique(as.character(seqnames(regionDB$regionGRL[[1]])))
genome.size <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[chrs]
genome.size <- sum(as.numeric(seqlengths(genome.size)))


binom.pvalues <- lapply(1:length(lola.regions), function(i.region) {
  print(i.region)
  regions <- lola.regions[[i.region]]
  anno.size <- sum(width(regions))/genome.size
  pvalues <- lapply(i.specific, function(i.spec){
    n.regions <- length(i.spec)
    ov <- findOverlaps(regions, universe[i.spec,])
    anno.hits <- length(unique(subjectHits(ov))) 
    pvalue <- binom.test(x = anno.hits, n = n.regions, p = anno.size)
    return(pvalue$p.value)
  })
  pvalues <- unlist(pvalues)
  return(pvalues)
})
