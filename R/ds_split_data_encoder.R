#' This function split the encoded data in training and test set.
#'
#' This function extracts features from the reference dataset (scRNA-seq)
#' @param features The encoder model output of the ds_encoder function.
#' @param clus A named factor with the cell type annotations from the reference data.
#' @param prop The proportion of cells for the training dataset (default=0.5).
#' @return A list with the following data: train_x= training_data,train_y= categorical cell type variables,
#' test_x =test_data,test_y=categorical celltype from the test data, classes=original cell types).
#' @export
#' @examples
#' #not run
#' enc <- ds_encoder(data = scaled,genes = gg,dims = 2000,verbose = T,hnodes = c(5000))
#' features <- ds_get_features(enc = enc,data = scaled,genes=gg)
#' out <- ds_split_data_encoder(features = features,clus = cluster,prop = 0.8,verbose = T)


ds_split_data_encoder <- function (features,clus,prop = NULL,
                        verbose = TRUE)
{
  classes <- levels(clus)
  levels(clus) <- seq(1:length(levels(clus)))
  if (verbose)
    message("Splitting the refence into train and test datasets...")
  total <- 6
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  Sys.sleep(0.1)
  progress <- 1
  setTxtProgressBar(pb, progress)
  start.time <- Sys.time()
  sizes <- unlist(lapply(levels(clus), function(x) length(clus[which(clus %in%
                                                                       x)])))
  names(sizes) <- levels(clus)
  progress <- 2
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  if (is.null(prop)) {
    prop <- 0.5
  }
  train.sample <- unlist(sapply(levels(clus), function(x) sample(x = which(clus %in%
                                                                             x), size = sizes[names(sizes) == x] * prop)))
  s <- rownames(features)[train.sample]
  test.sample <- which(!rownames(features) %in% s)
  train.data <- features[train.sample,]
  progress <- 3
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  test.data <- features[test.sample,]
  progress <- 4
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  out.train <- clus[train.sample]
  out.test <- clus[test.sample]

  #var.train <- apply(train.data, 1, function(x) (x - min(x))/(max(x) - min(x)))
  progress <- 5
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)

  #var.test <- apply(test.data, 1, function(x) (x - min(x))/(max(x) - min(x)))
  progress <- 6
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  model.train <- data.frame(out.train, train.data,check.names = F)
  n <- nrow(model.train)
  model.train <- model.train[sample(x=seq(1:n),size = n),]

  #model.test <- data.frame(out.test, var.test,check.names = F)

  library(keras)
  train_x <- as.matrix(model.train[,-1])
  train_y <- to_categorical(model.train[,1])
  test_x <- test.data
  test_y <- to_categorical(out.test)




  return(list(train_x=train_x,train_y=train_y,test_x =test_x,test_y=test_y,classes=classes))
}
