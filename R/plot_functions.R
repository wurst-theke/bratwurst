#==============================================================================#
#                      NMF-GPU plot generation - FUNCTIONS                     #
#==============================================================================#
#' Science like figures - ggplot2 theme
#'
#' @return
#' 
#' @import ggplot2
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

#' Universal ploting function for Optimal K Statistics
#'
#' @param nmf.exp 
#'
#' @return
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
plotKStats <- function(nmf.exp) {
  frobError.df <- melt(as.data.frame(FrobError(nmf.exp)))
  frobError.df[,1] <- factor(gsub('X', '', frobError.df[,1]))
  frobError.df <- data.frame('k' = frobError.df[,1], 
                             'variable' = 'FrobError',
                             'value' = frobError.df[,2])
  optKStats.df <- melt(as.data.frame(OptKStats(nmf.exp)), id.vars = 'k')
  
  meanError.df <- optKStats.df[optKStats.df$variable == 'mean',]
  meanError.df$variable <- unique(frobError.df$variable)
  plot.vars <- as.character(unique(optKStats.df$variable)[-1:-3])
  optKStats.df <- optKStats.df[optKStats.df$variable%in%plot.vars,]
  optKStats.df <- rbind(frobError.df, optKStats.df)
  optKStats.df$k <- as.numeric(as.character(optKStats.df$k))
  
  
  
  gg.optK <- ggplot() + geom_point(data = optKStats.df, aes(x=k, y=value), 
                                   col = 'black', size = 0.75)
  gg.optK <- gg.optK + geom_point(data = meanError.df,
                                  aes(x=k, y=value), col = 'red')
  gg.optK <- gg.optK + facet_wrap(~variable, scales = 'free_y')
  gg.optK <- gg.optK + xlab('K') + ylab('') + theme_bw() + science_theme()
  gg.optK <- gg.optK + theme(strip.background	= element_rect(fill = 'white'))
  return(gg.optK)
}

#' Ploting function for Rank vs. Frobenius Error
#'
#' @param nmf.exp 
#'
#' @return
#' 
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
plotRankedFrobErrors <- function(nmf.exp) {
  # Prepare frob error data for plotting.
  frobErrorRanks.df <- apply(FrobError(nmf.exp), 2, function(matrix.col) {
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
  gg.rankedFrobError <- gg.rankedFrobError + theme_bw() + theme_cowplot()
  gg.rankedFrobError <- gg.rankedFrobError + theme(strip.background = element_rect(fill = 'white'))
  return(gg.rankedFrobError)
}

#' Generate ColorMap for ComplexHeatmaps
#'
#' @param matrix.list 
#'
#' @return
#' 
#' @importFrom circlize colorRamp2
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

plotHeatmap4KMatrix <- function(nmf.exp, k){
  # Get matrix
  if (H) m <- HMatrix(nmf.exp, k = k)
  
  # Get signature names
  signature.names <- getSignatureNames(nmf.exp, k)
  
  # Define colormap for heatmaps
  col.map <- getColorMap(m)
  
  # Create annotations
  n.cols <- ncol(colData(nmf.exp))
  colnames(m) <- colData(nmf.exp)[,1]
  rownames(m) <- signature.names
  heat.anno <- HeatmapAnnotation(colData(nmf.exp)[,2:n.cols])
  
  heatmap.list <- Heatmap(matrix = m, 
                          col = col.map,
                          cluster_rows = F,
                          top_annotation = heat.anno,
                          clustering_distance_rows = 'pearson',
                          heatmap_legend_param = list(color_bar = 'continuous'),
                          column_title = sprintf('H-Matrix for K = %s', k))
  
  
}


#' Plot H or W-Matrix as heatmaps for a given range of K's
#'
#' @param matrix.list 
#' @param titles 
#' @param trans 
#' @param cluster_columns 
#' @param show.row 
#'
#' @return
#' 
#' @import ComplexHeatmap
#' @export
#'
#' @examples
plotHeatmap4MatrixList <- function(nmf.exp, H = T, W = F, titles = NULL, trans = F,
                                   cluster_columns = F, show.row = F) {
  # Get Matrix List
  matrix.list <- HMatrix(nmf.exp)
  
  # Define colormap for heatmaps
  col.map <- getColorMap(matrix.list)
  
  # Create annotations
  heat.anno <- HeatmapAnnotation(colData(nmf.exp))
  
  # Plot heatmaps
  heatmap.list <- Heatmap(matrix.list[[1]], 
                          col = col.map,
                          cluster_columns = F,
                          show_row_names = show.row,
                          clustering_distance_rows = 'pearson',
                          heatmap_legend_param = list(color_bar = 'continuous'),
                          
                          column_title = colData(nmf.exp))
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
