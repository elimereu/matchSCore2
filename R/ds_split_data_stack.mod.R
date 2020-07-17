#' This function split the reference data in validation and test set to use it as input for the ds_stacking_model function.
#' In this way the joint model will be fitted with data that are different from the original training set.
#' @param out The original output of the ds_split_data_dnn function for one of the model
#' @param clus A named factor with the cell type annotations from the reference data.
#' @param prop The proportion of cells for the training dataset (default=0.5).
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#' @return A list with the following data: val_x= validation_data,val_y= categorical cell type variables,
#' test_x =test_data,test_y=categorical celltype from the test data, classes=original cell types).
#' @export
#' @examples
#' #not run
#  out1 <- ds_split_data(scale.data = scaled,clus = cluster,genes = sel_gg,prop = 0.6)
#  outj <- ds_split_stack.mod(out = out1,clus = cluster,prop = 0.8)


ds_split_stack.mod <- function (out,clus,prop = NULL,
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
  test_x <- out$test_x
  clus <- clus[rownames(test_x)]
  sizes <- unlist(lapply(levels(clus), function(x) length(clus[which(clus %in%
                                                                       x)])))
  names(sizes) <- levels(clus)
  progress <- 2
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  if (is.null(prop)) {
    prop <- 0.5
  }
  val.sample <- unlist(sapply(levels(clus), function(x) sample(x = which(clus %in%
                                                                             x), size = sizes[names(sizes) == x] * prop)))
  s <- names(clus)[val.sample]
  test.sample <- names(clus)[which(!names(clus) %in% s)]

  progress <- 3
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)

  progress <- 4
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  out.val <- clus[val.sample]
  out.test <- clus[test.sample]


  progress <- 5
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)


  progress <- 6
  Sys.sleep(0.1)
  setTxtProgressBar(pb, progress)
  model.val <- data.frame(out.val,test_x[names(out.val),],check.names = F)
  n <- nrow(model.val)
  model.val <- model.val[sample(x=seq(1:n),size = n),]


  library(keras)
  val_x <- as.matrix(model.val[,-1])
  val_y <- to_categorical(model.val[,1])
  test_x <- test_x[out.test,]
  test_y <- to_categorical(out.test)



  return(list(val_x=val_x,val_y=val_y,test_x =test_x,test_y=test_y,classes=classes))
}
