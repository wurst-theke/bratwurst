#==============================================================================#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#==============================================================================#
# Define global plot theme.
require('ComplexHeatmap')
require('ggplot2')
require('cowplot')
#' Title
#'
#' @return
#' @export
#'
#' @examples
science_theme <- function() {
  science_theme <- theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                         axis.line = element_line(size = 0.7, color = "black"),
                         text = element_text(size = 14)) 
  science_theme <- science_theme + background_grid(major = "xy", minor = "none")
  return(science_theme)
}


# Define empty x-axis theme
#' Title
#'
#' @return
#' @export
#'
#' @examples
emptyXAxis_theme  <- function() {
  emptyXAxis_theme <- theme(axis.ticks.x = element_blank(),
                            axis.line.x  = element_blank(),
                            axis.title.x  = element_blank(),
                            axis.text.x = element_blank())
  return(emptyXAxis_theme)
}

# Get Colormap
#' Title
#'
#' @param matrix.list 
#'
#' @return
#' @export
#'
#' @examples
getColorMap <- function(matrix.list) {
  matrix.max <- quantile(unlist(matrix.list), 0.95)
  matrix.min <- quantile(unlist(matrix.list), 0.05)
  if(matrix.min > 0) {
    col.map <- colorRamp2(c(0, matrix.max), c("white", "red"))
  } else {
    col.map <- colorRamp2(c(matrix.min, 0, matrix.max), c("blue", "white", "red"))
  }
  return(col.map)
}

#' Title
#'
#' @param frob.errorMatrix 
#'
#' @return
#' @export
#'
#' @examples
plotRankedFrobErrors <- function(frob.errorMatrix) {
  # Prepare frob error data for plotting.
  frobErrorRanks.df <- apply(frobError.matrix, 2, function(matrix.col) {
    ranks <- order(matrix.col, decreasing = T)
    rank.df <- data.frame(rank = 1:length(ranks), 
                          froberror = matrix.col[ranks])
    return(rank.df)
  })
  frobErrorRanks.df <- melt(frobErrorRanks.df, value.name = 'frobError',
                            id.vars = 'rank')
  grid.labels <- paste('K =', frobErrorRanks.df$L1)
  frobErrorRanks.df$L1 <- factor(grid.labels, levels = unique(grid.labels))
  # Plot ranked frob error
  gg.rankedFrobError <- ggplot(frobErrorRanks.df, aes(x = rank, y = frobError))
  gg.rankedFrobError <- gg.rankedFrobError + geom_line(col = 'black', size = 0.75)
  gg.rankedFrobError <- gg.rankedFrobError + facet_wrap(~L1, scales = 'free')
  gg.rankedFrobError <- gg.rankedFrobError + xlab('Error rank') + ylab('Frobenius error')
  gg.rankedFrobError <- gg.rankedFrobError + theme_bw() + science_theme 
  gg.rankedFrobError <- gg.rankedFrobError + theme(strip.background = element_rect(fill = 'white'))
  return(gg.rankedFrobError)
}

#' Title
#'
#' @param frob.errorMatrix 
#' @param keep.x 
#'
#' @return
#' @export
#'
#' @examples
plotFrobError <- function(frob.errorMatrix, keep.x = F) {
  frobError.melt <- melt(frob.errorMatrix)
  frobError.melt$variable <- as.numeric(as.character(frobError.melt$variable))
  frobError.data <- computeFrobErrorStats(frob.errorMatrix)  
  
  # Frob error Plot
  gg.frobErrorScatter <- ggplot() + geom_point(data = frobError.melt,
                                               aes(x=variable, y=value), 
                                               col = 'black', size = 0.5)
  gg.frobErrorScatter <- gg.frobErrorScatter + geom_point(data = frobError.data,
                                                          aes(x=k, y=mean),
                                                          col = 'red')
  gg.frobErrorScatter <- gg.frobErrorScatter + ylab('Frobenius Error') 
  gg.frobErrorScatter <- gg.frobErrorScatter + theme_bw() + science_theme 
  gg.frobErrorScatter <- gg.frobErrorScatter + emptyXAxis_theme
  
  # coef of variation = sd/mean
  gg.frobErrorCv <- ggplot() + geom_point(data = frobError.data, 
                                          aes(x=k, y=cv), col = 'darkgreen')
  gg.frobErrorCv <- gg.frobErrorCv + xlab('Factorization rank') 
  gg.frobErrorCv <- gg.frobErrorCv + ylab('Coef. of variation')
  gg.frobErrorCv <- gg.frobErrorCv + theme_bw() + science_theme
  if(!keep.x) {
    gg.frobErrorCv <- gg.frobErrorCv + emptyXAxis_theme
  }
  
  
  gg.frobError <- plot_grid(gg.frobErrorScatter, gg.frobErrorCv,
                            nrow = 2, align = 'v')
  
  return(gg.frobError) 
}

#' Title
#'
#' @param sil.vec 
#' @param keep.x 
#'
#' @return
#' @export
#'
#' @examples
plotSilhoutteStas <- function(sil.vec, keep.x = F) {
  sil.melt.df <- melt(sil.vec)
  sil.melt.df$L1 <- as.numeric(as.character(sil.melt.df$L1))
  Sum.df <- sil.melt.df[which(sil.melt.df$L2 == "sum"),]
  Mean.df <- sil.melt.df[which(sil.melt.df$L2 == "mean"),]
  
  # Scatterplot for Sum Silhoutte values.
  gg.silScatterSum <- ggplot() + geom_point(data = Sum.df,  aes(x=L1, y=value),
                                            col = 'blue')
  gg.silScatterSum <- gg.silScatterSum + xlab('Factorization rank') 
  gg.silScatterSum <- gg.silScatterSum + ylab('Sum over\nSilhouttes')
  gg.silScatterSum <- gg.silScatterSum + theme_bw() + science_theme
  gg.silScatterSum <- gg.silScatterSum + emptyXAxis_theme
  
  # Scatterplot for mean Silhoutte values.
  gg.silScatterMean <- ggplot() + geom_point(data = Mean.df, aes(x=L1, y=value),
                                             col = 'blue')
  gg.silScatterMean <- gg.silScatterMean + xlab('Factorization rank') 
  gg.silScatterMean <- gg.silScatterMean + ylab('Mean over\nSilhouttes')
  gg.silScatterMean <- gg.silScatterMean + theme_bw() + science_theme
  if(!keep.x) {
    gg.silScatterMean <- gg.silScatterMean + emptyXAxis_theme 
  }
  
  gg.silScatter <- plot_grid(gg.silScatterSum, gg.silScatterMean,
                             nrow = 2, align = 'v')
  return(gg.silScatter)
}

#' Title
#'
#' @param coph.coeff 
#' @param keep.x 
#'
#' @return
#' @export
#'
#' @examples
plotCopheneticCoeff <- function(coph.coeff, keep.x = F){
  coph.melt.df <- melt(coph.coeff)
  coph.melt.df$L1 <- as.numeric(as.character(coph.melt.df$L1))
  gg.CophScatterMean <- ggplot() + geom_point(data = coph.melt.df, 
                                              aes(x=L1, y=value), col = 'blue')
  gg.CophScatterMean <- gg.CophScatterMean + xlab('Factorization rank')
  gg.CophScatterMean <- gg.CophScatterMean + ylab('Cophenetic\ncorr. coeff.')
  gg.CophScatterMean <- gg.CophScatterMean + theme_bw() + science_theme
  if(!keep.x) {
    gg.CophScatterMean <- gg.CophScatterMean + emptyXAxis_theme 
  }
  return(gg.CophScatterMean)
}

#' Title
#'
#' @param matrix.list 
#' @param titles 
#' @param trans 
#' @param cluster_columns 
#' @param show.row 
#'
#' @return
#' @export
#'
#' @examples
plotHeatmap4MatrixList <- function(matrix.list, titles = NULL, trans = F,
                                   cluster_columns = F, show.row = F) {
  if(trans) {
    matrix.list <- lapply(matrix.list, function(m) t(m))
  }
  matrix.max <- quantile(unlist(matrix.list), 0.95)
  matrix.min <- quantile(unlist(matrix.list), 0.05)
  if(matrix.min > 0) {
    col.map <- colorRamp2(c(0, matrix.max), c("white", "red"))
  } else {
    col.map <- colorRamp2(c(matrix.min, 0, matrix.max), c("blue", "white", "red"))
  }
  heatmap.list <- Heatmap(matrix.list[[1]], 
                          col = col.map,
                          cluster_columns = F,
                          show_row_names = show.row,
                          clustering_distance_rows = 'pearson',
                          heatmap_legend_param = list(color_bar = "continuous"),
                          column_title = titles[1])
  for (i in 2:length(matrix.list)){
    heatmap.list <- heatmap.list + Heatmap(matrix.list[[i]],
                                           col = col.map,
                                           cluster_columns = F,
                                           show_row_names = show.row,
                                           clustering_distance_rows = 'pearson',
                                           show_heatmap_legend = FALSE, 
                                           column_title = titles[i])
  }
  return(heatmap.list)
}
