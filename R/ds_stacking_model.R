#' This function creates a stacking ensemble of 2 or more DNN models.
#'
#' This function creates a DNN model from the reference dataset (scRNA-seq) by using the keras package.
#' @param list.models A named list with the different models. The names of the models must be different and correspond to the name given to the name_mod in the ds_dnn_model function.
#' @param hnodes The number of the nodes of the intermediate layers of the encoder. The first and last layers are not included.
#' @param activation Keras parameter for the layer_dense function (Default= "relu").
#' @param add_dropout Keras parameter for the layer_dense function (Default= FALSE).
#' @param pct_dropout Keras parameter for the layer_dense function (Default= 0.2).
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param ... Other keras parameters.
#' @return The stacked model.
#' @export
#' @examples
#' #not run
#' stack.mod <- ds_stacking_model(list.models,hnodes = c(1000),activation = "relu",verbose = T)




ds_stacking_model <- function(list.models,hnodes,activation="relu",add_dropout=F,
                           pct_dropout=0.2,verbose=T, ...){

  for(i in c(1:length(list.models))){
    freeze_weights(list.models[[i]])
  }

  ensemble_input <- lapply(list.models,function(x) x$input)
  ensemble_out <- lapply(list.models, function(x) x$output)
  merge <- layer_concatenate(ensemble_out,name = "stack_mod")


  N <- length(hnodes)
    for(i in c(1:N)){
      hidden <- merge %>%
        layer_dense(units = hnodes[i], activation = activation,name = paste("stack_mod","hidden",i,sep="_"))
      if(add_dropout==TRUE){
        hidden <- hidden %>%
          layer_dropout(rate = pct_dropout,name = paste("stack_mod","dropout",i,sep="_"))
      }
    }

  out.shape <- list.models[[1]]$output$shape[[2]]
  output <- layer_dense(object = hidden,units = out.shape,activation = "softmax",name = "stacked_ensemble_out")

  model <- keras_model(inputs = ensemble_input,outputs = output)

  compile(model,optimizer = "adam",loss = "categorical_crossentropy",metrics = "accuracy")

  return(model)
}
