#' This function trains, fits and evaluates a deep neural network (DNN) model.
#'
#' This function creates a DNN model from the reference dataset (scRNA-seq) by using the keras package.
#' @param out The output of the ds_split_data_encoder function. The input data used to train the DNN model.
#' @param hnodes The number of the nodes of the intermediate layers of the encoder. The first and last layers are not included.
#' @param epochs Keras parameter for the layer_dense function (Default= 30).
#' @param batch_size Keras parameter for the layer_dense function (Default= 32).
#' @param activation Keras parameter for the layer_dense function (Default= "relu").
#' @param add_dropout Keras parameter for the layer_dense function (Default= TRUE).
#' @param pct_dropout Keras parameter for the layer_dense function (Default= 0.2).
#' @param name_mod name of the model. This is mandatory if you want to include the model in a stacked model.
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param ... Other keras parameters.
#' @return The DNN model.
#' @export
#' @examples
#' #not run
#' enc <- ds_encoder(data = scaled,genes = gg,dims = 2000,verbose = T,hnodes = c(5000))
#' features <- ds_get_features(enc = enc,data = scaled,genes=gg)
#' out <- ds_split_data_encoder(features = features,clus = cluster,prop = 0.8,verbose = T)
#' model <- ds_dnn_model(out = out,hnodes = c(1000),verbose = T,epochs = 10,batch_size = 32)




ds_dnn_model <- function(out,hnodes,epochs = 30,batch_size = 32,activation="relu",add_dropout=T,
         pct_dropout=0.2,name_mod="mod",verbose=T, ...){

  library(keras)

  train_x <- out$train_x
  train_y <- out$train_y

  test_x <- out$test_x
  test_y <- out$test_y
  Nclass <- ncol(train_y)

  model <- keras_model_sequential() %>%
    layer_dense(units = hnodes[1], activation = activation, input_shape = ncol(train_x),name = paste(name_mod,"dense_1",sep="_"))

  if(add_dropout==TRUE){
    model <- model %>%
      layer_dropout(rate = pct_dropout,name = paste(name_mod,"dropout_1",sep="_"))
  }


  N <- length(hnodes)
  for(i in c(1:N)){
    model <- model %>%
      layer_dense(units = hnodes[i], activation = activation,name = paste(name_mod,"dense",i+1,sep="_"))
    if(add_dropout==TRUE){
      model <- model %>%
        layer_dropout(rate = pct_dropout,name = paste(name_mod,"dropout",i+1,sep="_"))
    }
  }

  model <- model %>% layer_dense(units = Nclass,activation = "softmax",name = paste(name_mod,"out",sep="_")) %>%

    # backpropagation
    compile(
      optimizer = "adam",
      loss = "categorical_crossentropy",
      metrics = c("accuracy")
    )

  learn <- model %>% fit(
    x = train_x,
    y = train_y,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = .2,
    verbose = verbose
  )

  if(verbose==TRUE){
    eval <- model %>% evaluate(test_x, test_y)
    print("The accuracy of the model is: ")
    print(eval$acc)
  }


  return(model)

}
