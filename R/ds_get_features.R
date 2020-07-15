#' Feature extraction via an autoencoder model
#'
#' This function extracts features from the reference dataset (scRNA-seq)
#' @param enc The encoder model output of the ds_encoder function.
#' @param data A scaled matrix of gene expressions like in the `scale.data`
#' of the Seurat object. Rows are genes and columns are cells from the reference
#' dataset.
#' @param genes The set of genes you want to use before the dimensionality reduction. Usually, the genes are those in common with the query dataset.
#' @return The features extracted by the encoder.
#' @export
#' @examples
#' #not run
#' enc <- ds_encoder(data = scaled,genes = gg,dims = 2000,verbose = T,hnodes = c(5000))
#' features <- ds_get_features(enc = enc,data = scaled,genes=gg)


ds_get_features <- function(enc,data,genes){

  library(keras)
  data_x <- t(as.matrix(data[genes,]))
  features <- enc %>% predict(data_x)
  rownames(features) <- colnames(data)


  return(features)
  }
