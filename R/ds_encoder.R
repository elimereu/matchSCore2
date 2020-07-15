#' Creation of an autoencoder model
#'
#' This function creates an autoencoder from the reference dataset (scRNA-seq)
#'
#' @param data A scaled matrix of gene expressions like in the `scale.data`
#' of the Seurat object. Rows are genes and columns are cells from the reference
#' dataset.
#' @param genes The set of genes you want to use before the dimensionality reduction. Usually, the genes are those in common with the query dataset.
#' @param dims The dimensions of the space you want to reduce the original data.
#' @param hnodes The number of the nodes of the intermediate layers of the encoder. The first and last layers are not included.
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param activation Keras parameter for the layer_dense function (Default= "tanh").
#' @param ... Keras parameters for the layer_dense function.
#' @return The encoder model to reduce the dimensions of the input data.
#'
#' @export
#'
#' @examples
#' #not run
#' #enc <- encoder(data = scaled,genes = gg,dims = 2000,verbose = T,hnodes = c(5000))


ds_encoder <- function(data,genes,dims,hnodes,verbose=T,activation="tanh",...){

  library(keras)

  enc1 <- keras_model_sequential()  %>%
    layer_dense(units = hnodes[1], activation = activation, input_shape = nrow(data))

  nleyers <- length(hnodes)
  N <- nleyers -1
  if(nleyers>1){
    for(i in c(2:N)){
      enc1 <- enc1 %>%
        layer_dense(units = hnodes[i], activation = activation)
    }}

  enc <- enc1 %>%
    layer_dense(units = dims, activation = "tanh")




return(enc)
}
