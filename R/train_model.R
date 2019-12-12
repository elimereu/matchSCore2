#' Model training
#'
#' This function train a model from the reference dataset
#' @param scale.data A scaled matrix of gene expressions like in the `scale.data`
#' of the Seurat object. Rows are genes and columns are cells from the reference
#' dataset.
#' @param clus A factor with identities from the reference dataset.
#' @param gene_cl.ref A list of cluster-specific markers. Each element of the list
#' contains markers of a cell type. Usually only top100 ranked markers are used.
#' @param prop The proportion of the reference data used to train the model.
#' Default=0.5.
#' @param p.threshold Probability threshold to consider a cell classified.
#' Default=0.65.
#'
#' @return A multinomial fitted model as in the `nnet` package.
#'
#' @export
#'
#' @examples
#' # TODO
train_model <- function(scale.data,clus,gene_cl.ref,prop=NULL,p.threshold=NULL,...){

  print("Splitting the refence into train and test datasets...")

  total <- 10
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)


  Sys.sleep(0.1)
  # update progress bar
  progress <- 1
  setTxtProgressBar(pb,progress)

  start.time <- Sys.time()
  sizes <- unlist(lapply(levels(clus),function(x) length(clus[which(clus %in% x)])))
  names(sizes) <- levels(clus)

  progress <- 2
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)

  if(is.null(prop)){prop <- 0.5}

  train.sample <- unlist(sapply(levels(clus),function(x) sample(x = which(clus %in% x),size = sizes[names(sizes)==x]*prop)))
  s <- colnames(scale.data)[train.sample]
  test.sample <- which(!colnames(scale.data) %in% s)
  train.data <- scale.data[,train.sample]


  progress <- 3
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  test.data <- scale.data[,test.sample]

  progress <- 4
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  out.train <- clus[train.sample]
  out.test <- clus[test.sample]

  var.train <- sapply(gene_cl.ref,function(x) Matrix::colSums(train.data[x,]))
  var.train <- t(apply(var.train,1, function(x) (x-min(x))/(max(x)-min(x))))

  progress <- 5
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  var.test <- sapply(gene_cl.ref,function(x) Matrix::colSums(test.data[x,]))
  var.test <- t(apply(var.test,1, function(x) (x-min(x))/(max(x)-min(x))))

  progress <- 6
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  model.train <- data.frame(out.train,var.train)
  model.test <- data.frame(out.test,var.test)

  progress <- 7
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  require(Matrix)
  library(nnet)

  cat("\n Learning the model from the training dataset...\n")

  mod <- multinom(out.train ~ ., data = model.train,decay=0.0001,maxit = 500)
  print(summary(mod))

  progress <- 8
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  fitted.results <- predict(mod, newdata = model.test, "probs")

  cat("\n Fitting identities in the set data.. \n")

  fit <- apply(fitted.results,1,function(x) colnames(fitted.results)[which(x==max(x))])
  if(is.null(p.threshold)){p.threshold <- 0.65}
  fit <- apply(fitted.results,1,function(x) ifelse(max(x)>p.threshold,colnames(fitted.results)[which(x==max(x))],"unclassified"))
  fit_res <- data.frame(out.test,fit,fitted.results)

  progress <- 9
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  cat("\n Measuring the accuracy of the fitting..\n")

  acc <- length(which(as.character(fit_res$out.test)==as.character(fit_res$fit)))/length(fit_res$out.test)

  cat("\n",paste("The accuracy of the model is:",round(acc,digits = 2),sep=" "),"\n")

  print(table(id_test=fit_res$out.test,class=fit_res$fit))

  progress <- 10
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)
  close(pb)
  end.time <- Sys.time()
  time <- difftime(end.time,start.time,units="mins")
  print(paste("The runtime is:",time,"min",sep=" "))

return(mod)

  }
