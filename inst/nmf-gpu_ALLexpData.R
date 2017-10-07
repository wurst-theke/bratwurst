# Copyright © 2015-2017  The Bratwurst package contributors
# This file is part of the Bratwurst package. The Bratwurst package is licenced
# under GPL-3

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
lapply(seq(2, k.max), function(k) {
  plotHMatrix(nmf.exp, k)
})

#==============================================================================#
#                        W-MATRIX ANALYSIS/VISUALIZATION                       #
#                       FIND REGIONS ASSOCIATED TO SIGNATURES                  #
#==============================================================================#
### Find representative regions.
# Get W for best K
nmf.exp <- setOptK(nmf.exp, 4)
OptK(nmf.exp)

signature.names <- getSignatureNames(nmf.exp, OptK(nmf.exp))

FeatureStats(nmf.exp)

nmf.exp <- computeFeatureStats(nmf.exp)

#nmf.exp <- computeEntropy4OptK(nmf.exp)

#nmf.exp <- computeAbsDelta4OptK(nmf.exp)

FeatureStats(nmf.exp)

#==============================================================================#
#      Analysis of Signature combination for regulatory element enrichment.    #
#==============================================================================#
# Plot all possible signature combinations
plotSignatureFeatures(nmf.exp)

# Plot only signature combinations
plotSignatureFeatures(nmf.exp, sig.combs = F)

sig.id <- '1000'
m <- WMatrix(nmf.exp, k = OptK(nmf.exp))[FeatureStats(nmf.exp)[,1] == sig.id,]
m <- m[order(m[,1]),]
c <- getColorMap(m)
#m <- t(apply(m, 1, function(r) (r - mean(r))/sd(r)))

Heatmap(m, col = c, cluster_rows = T, cluster_columns = F)

# genes <- rowData(nmf.exp)$V3[FeatureStats(nmf.exp)$cluster == '0001']
# write.table(genes, file = 'test.txt', quote = F, sep = '\n', row.names = F, col.names = F)
# write.table(rowData(nmf.exp)$V3, file = 'test_all.txt', quote = F, sep = '\n', row.names = F, col.names = F)
