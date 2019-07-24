#' Cell projections onto the reference dataset
#'
#' This function infers cell identities by using the model learned by the reference dataset. 
#' @param scale.data A scaled matrix of gene expressions like in the scale.data of the Seurat object. Rows are genes and columns are cells from the reference dataset.
#' @param model The model as in output of the train_model function.
#' @param gene_cl.ref A list of cluster-specific markers. Each element of the list contains markers of a cell type. Usually only top100 ranked markers are used.
#' @param p.threshold Probability threshold to consider a cell classified. Default=0.5.
#' @return A list with a vector with cell identities predicted by the model (ids) and a dataframe with probabilities at each identity class.
#' @export
#' @examples
#' 


identity_map <- function(scale.data,model,gene_cl.ref,p.threshold=NULL){
  
  require(Matrix)
  p<- lapply(gene_cl.ref,function(x) rownames(scale.data)[which(rownames(scale.data) %in% x)])
  
  var.test <- sapply(p,function(x) Matrix::colSums(scale.data[x,]))
  var.test <- t(apply(var.test,1, function(x) (x-min(x))/(max(x)-min(x))))
  
  model.test <- data.frame(var.test)
  fitted.results <- predict(model, newdata = model.test, "probs")
  
  fit <- apply(fitted.results,1,function(x) colnames(fitted.results)[which(x==max(x))])
  if(is.null(p.threshold)){p.threshold <- 0.5}
  fit <- apply(fitted.results,1,function(x) ifelse(max(x)>p.threshold,colnames(fitted.results)[which(x==max(x))],"unclassified"))
  
  
  return(list(ids=fit,fit.prob=data.frame(fitted.results)))
}