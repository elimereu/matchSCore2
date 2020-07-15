#' This function fits the stacking ensemble of 2 or more DNN models.
#'
#' @param model The stacking model output of the function ds_stacking_model
#' @param list.models A named list with the different models. The names of the models must be different and correspond to the name given to the name_mod in the ds_dnn_model function.
#' @param train_x A new_train dataset to fit the parameter of the model. Usually part of the test_set can be used.
#' @param activation Keras parameter for the layer_dense function (Default= "relu").
#' @param add_dropout Keras parameter for the layer_dense function (Default= FALSE).
#' @param pct_dropout Keras parameter for the layer_dense function (Default= 0.2).
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @param ... Other keras parameters.
#' @export
#' @examples
#' #not run
#' enc <- ds_encoder(data = scaled,genes = gg,dims = l,verbose = T,hnodes = c(5000))
#' features <- ds_get_features(enc = enc,data = scaled,genes=gg)
#' out <- ds_split_data_encoder(features = features,clus = cluster,prop = 0.8,verbose = T)
#' mod1 <- ds_dnn_model(out = out,hnodes = c(100),verbose = T,epochs = 10,batch_size = 32,name_mod = "mod1")
#' mod2 <- ds_dnn_model(out = out,hnodes = c(200),verbose = T,epochs = 10,batch_size = 32,name_mod = "mod2")
#' mod3 <- ds_dnn_model(out = out,hnodes = c(800),verbose = T,epochs = 10,batch_size = 32,name_mod = "mod3")
#' list.models <- list(mod1,mod2,mod3)
#' new_train_x <- out$test_x
#' new_train_y <- out$test_y
#' stack.mod <- ds_stacking_model(list.models,hnodes = c(1000),activation = "relu",verbose = T)
#' ds_fit_stacking_model(model = stack.mod,list.models = list.models,train_x = new_train_x,train_y = new_train_y,epochs = 10,
#'                   batch_size = 32,validation_split = 0.5)


ds_fit_stacking_model <- function(model,list.models,train_x,train_y,epochs = 10,batch_size = 32,validation_split=0.2,
                               verbose=T,...){

  l <- length(list.models)
  inputX <- lapply(c(1:l), function(x) train_x)

  learn <- model %>% fit(
    x = inputX,
    y = train_y,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = validation_split,
    verbose = verbose
    )

}
