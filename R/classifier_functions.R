#' Generate classifier
#'
#' Details
#'
#' @param matrix_W
#'
#' @return k.row, a list of lists (one list per row of the input matrix)
#' @export
#'
#' @examples
generate_classifier <- function(matrix_W){
  k.row <- apply(matrix_W, 1, function(x) {
    x.trans <- x
    k.cluster <- kmeans(x.trans, 2)
    d <- dist(as.data.frame(x.trans))
    sil.mean <- mean(silhouette(k.cluster$cluster, d)[, 3])
    size <- k.cluster$size
    return(list(data = x,
                centers = t(k.cluster$centers),
                silmean = sil.mean,
                explainedVar = k.cluster$betweenss / k.cluster$totss,
                oddsVar = sum(k.cluster$withinss) / k.cluster$betweenss,
                attribution = k.cluster$cluster,
                size = size))
  })
  return(k.row)
}

#' Filter features which are specific for a single signature
#'
#' @param k.row
#' @param specific
#'
#' @return
#' @export
#'
#' @examples
filter_specific_features <- function(k.row, specific=TRUE){
  if(isTRUE(specific)){
    specific_list <- lapply(k.row, FUN = function(t){
      if(min(t$size) == 1){
        spec_temp <- which(t$size == 1)
        unspec_temp <- which(t$size != 1)
        if(t$centers[spec_temp] > t$centers[unspec_temp]){
          x <- t$data
        } else {
          x <- as.vector(matrix(0, 1, length(t$data)))
        }
      } else {
        x <- as.vector(matrix(0, 1, length(t$data)))
      }
      return(x)
    })
  }else{
    specific_list <- lapply(k.row, FUN = function(t){
      if(min(t$size) == 1){
        x <- t$data
      }else{
        x <- as.vector(matrix(0, 1, length(t$data)))
      }
      return(x)
    })
  }
  specific_df <- do.call(rbind, specific_list)
  return(specific_df)
}

#' Filter features
#'
#' Filter features which are specific for a single signature if min(size) == 1
#' is already set
#'
#' @param k.row
#' @param out_int
#' @param specific_signatures
#' @param genes_df
#'
#' @return
#' @export
#'
#' @examples
filter_features <- function(k.row,
                            out_int,
                            specific_signatures,
                            genes_df){
  k.row <- k.row[out_int]
  k.row <- k.row[specific_signatures]
  specific_list <- lapply(k.row, FUN = function(t){
    spec_temp <- which(t$size == 1)
    unspec_temp <- which(t$size != 1)
    if(t$centers[spec_temp] > t$centers[unspec_temp]){
      x <- t$data
    }else{
      x <- as.vector(matrix(0, 1, length(t$data)))
    }
    return(x)
  })
  temp_df <- do.call(rbind, specific_list)
  rowsums <- rowSums(temp_df)
  temp_int <- which(rowsums != 0)
  reduced_signatures <- genes_df[temp_int, ]
  return(reduced_signatures)
}

#' Match and count
#'
#' @param LCD_result_df
#'
#' @return
#' @export
#'
#' @examples
match_and_count <- function(LCD_result_df){
  counter <- 0
  match_value <- 0
  for(x in 1:dim(LCD_result_df)[2]){
    max_col <- which.max(LCD_result_df[, x])
    temp_col <- gsub("\\.[0-9]*$", "", colnames(LCD_result_df)[x])
    temp_row <- gsub("\\.[0-9]*$", "", rownames(LCD_result_df)[max_col])
    if(temp_col == temp_row){
      counter <- counter + LCD_result_df[max_col, x]
      match_value <- match_value + 1
    }else{
      counter <- counter + 0
      match_value <- match_value + 0
    }
  }
return(c(counter = counter, match = match_value))
}

#' Calculate F1 score
#'
#' @param LCD_result_df
#'
#' @return
#'
#' @importFrom caret confusionMatrix
#' @export
#'
#' @examples
calculate_FI <- function(LCD_result_df){
  conf_data <- vector()
  loop_counter <- 0
  signatures <- gsub("\\.[0-9]*$", "", rownames(LCD_result_df))
  number_of_signatures <- length(unique(signatures))
  for(x in 1:dim(LCD_result_df)[2]){
    loop_counter <- loop_counter + 1
    max_col <- which.max(LCD_result_df[, x])
    temp_row <- gsub("\\.[0-9]*$", "", rownames(LCD_result_df)[max_col])
    conf_data[loop_counter] <- temp_row
  }

  reference_data <- gsub("\\.[0-9]*$", "", colnames(LCD_result_df))

  if(length(unique(conf_data)) == number_of_signatures){
    conf_matrix <- confusionMatrix(conf_data, reference_data)
    class_list <- as.list(c(1:number_of_signatures))
    F1_list <- lapply(class_list, FUN = function(x){
      recall <- conf_matrix$byClass[x, 1]
      precision <- conf_matrix$byClass[x, 3]
      FI_temp <- 2 * ((recall * precision) / (recall + precision))
      return(FI_temp)
    })
    names_FI <- gsub("Class: ", "", rownames(conf_matrix$byClass))
    names(F1_list) <- names_FI
    FI_vector <- unlist(F1_list)
    FI <- mean(unlist(F1_list))
  } else {
    FI_vector <- 0
    FI <- 0
    warning(paste0("ConfusionMatrix could not been used: levels of ",
                   "reference exceed levels of matches"))
  }
  return(list(FI_vector = FI_vector, FI = FI))
}


#==============================================================================#
#                   Caluclate F1 etc for a range of classfiers                 #
#==============================================================================#

#' Title
#'
#' @param signature_matrix
#' @param matrix_V
#' @param test_odds
#' @param test_dist
#' @param oddsVar_min
#' @param k.dist_max
#'
#' @return
#' @export
#'
#' @examples
calculate_for_classifier_range <- function(signature_matrix,
                                           matrix_V,
                                           test_odds,
                                           test_dist,
                                           oddsVar_min = 0,
                                           k.dist_max = 100000){
  all_results_list <- lapply(test_odds, FUN = function(thres.outer){
    results_list <- lapply(test_dist, FUN = function(thres.inner){
      all.signature <- which(k.dist > thres.inner[1] &
                               k.dist < k.dist_max &
                               oddsVar < thres.outer[1] &
                               oddsVar > oddsVar_min)
      if(length(all.signature) >= 2){
        #####SAMPLE MATRIX
        sample_matrix <- matrix_V[all.signature, ]
        sample_matrix <- data.frame(apply(sample_matrix, 2,
                                          function(x) x / sum(x)))
        #####SIGNATURE MATRIX
        signature_matrix <- signature_matrix[all.signature, ]
        signature_matrix <- data.frame(apply(signature_matrix, 2,
                                             function(x) x / sum(x)))
        current_cutoff_vector <- rep(0, dim(signature_matrix)[2])
        LCD_abs <- try(LCD_complex_cutoff(
          sample_matrix,
          signature_matrix,
          in_cutoff_vector = current_cutoff_vector,
          in_filename = NULL,
          in_method = "abs"))
        if(class(LCD_abs) == "try-error"){
          counter <- 0
          match_value <- 0
          sample_size <- 1
          FI <- 0
          FI_vector <- 0
        } else {
          LCD_result_df <- LCD_abs$exposures
          #####Evaluation of run
          ## Match and Count
          m_and_c_vector <- match_and_count(LCD_result_df)
          counter <- m_and_c_vector[1]
          match_value <- m_and_c_vector[2]
          ##Confusion matrix
          F1_result_list <- calculate_FI(LCD_result_df)
          FI_vector <- F1_result_list$FI_vector
          FI <- F1_result_list$FI
          sample_size <- dim(LCD_result_df)[2]
        }
      } else {
        counter <- 0
        match_value <- 0
        sample_size <- 1
        FI <- 0
        FI_vector <- 0
      }
      results <- list(Counter = (counter / match_value),
                      Match = (match_value / sample_size),
                      Number_of_genes = (length(all.signature)),
                      F1_for_each_class = FI_vector,
                      F1_mean = FI)
      return(results)
    })
    names <- unlist(test_dist)
    names <- paste0("dist > ", round(names, 4))
    counter <- unlist(lapply(results_list, FUN = function(x) x[1]))
    names(counter) <- names
    matches <- unlist(lapply(results_list, FUN = function(x) x[2]))
    names(matches) <- names
    number <- unlist(lapply(results_list, FUN = function(x) x[3]))
    names(number) <- names
    FI_vector <- lapply(results_list, FUN = function(x) x[4])
    names(FI_vector) <- names
    FI <- unlist(lapply(results_list, FUN = function(x) x[5]))
    names(FI) <- names
    all_results <- list(Counter = counter, Matches = matches, Numbers = number,
                        FI_vector = FI_vector, FI_mean = FI)
    return(all_results)
  })
  names <- unlist(test_odds)
  names(all_results_list) <- paste0("oddsVar < ", round(names, 4))
  return(all_results_list)
}


#==============================================================================#
#                   EXTRACT GENES AND BUILD GRANGES OBJECT                     #
#==============================================================================#

#' Extract genes specfic for each signature
#'
#' @param genes_df
#' @param annotation_df
#' @param match_by
#'
#' @return
#' @export
#'
#' @examples
extract_signature_specific_genes <- function(genes_df, annotation_df,
                                             match_by = "external_gene_name"){
  if(dim(genes_df)[1] == dim(annotation_df)[1]){
    genes_df <- genes_df
  }else{
    rownames(genes_df) <- gsub("\\.[0-9]*$", "", rownames(genes_df))
    col_int <- grep(paste0(match_by), colnames(annotation_df))
    match_int <- match(tolower(annotation_df[, col_int]),
                       tolower(rownames(genes_df)))
    genes_df <- genes_df[match_int, ]
  }
  temp_list <- as.list(colnames(genes_df))
  signature_specfic <- lapply(temp_list, FUN = function(signature){
    temp_int <- c(0)
    for(i in 1:dim(genes_df)[1]){
      max_row <- which.max(genes_df[i, ])
      if(names(max_row) == signature){
        temp_int <- cbind(temp_int, i)
      }
    }
    temp_int <- as.vector(temp_int)
    temp_int <- temp_int[-1]
    specific_genes <- annotation_df[temp_int, ]
    return(specific_genes)
  })
  names(signature_specfic) <- c(unlist(temp_list))
  return(signature_specfic)
}

#' Build GRanges object from signature list
#'
#' @param specific_signatures_list
#'
#' @return
#' @export
#'
#' @examples
build_granges_object <- function(specific_signatures_list){
  gr_list <- lapply(specific_signatures_list, FUN = function(x){
    x$strand <- gsub("^1$", "+", x$strand)
    x$strand <- gsub("^-1$", "-", x$strand)
    x$chromosome_name <- paste0("chr", x$chromosome_name)
    gr <- GRanges(seqnames = Rle(c(x$chromosome_name)),
                  ranges = IRanges(x$start_position, x$end_position),
                  strand = Rle(strand(c(x$strand))),
                  genes = c(x$external_gene_name))
    return(gr)
  })
  names(gr_list) <- names(specific_signatures_list)
  return(gr_list)
}
