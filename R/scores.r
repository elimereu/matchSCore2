#' Compute matchSCore for different thresholds of specificity and ntop ranked genes
#'
#' This function computes the matchSCore for different cutoffs of marker specificity and top cluster markers
#' @param sim A Splatter simulation object.
#' @param tool_out The output of the used tool to cluster your data (It should be as the output of the seurat_run function).
#' @param ntop A vector with different proportions of top-ranked markers for each cluster. Usually seq(250,2000,250).
#' @param tool_run The name of the function used to run the clusters without quote (Default= seurat_run).
#' @param labels Cluster labels as in the output of the compute_labels function.
#' @return A matrix. Each row contains the matchSCore at the specified level of specificity and ntop.  
#' @export
#' @examples
#' 


scores <- function(sim,tool_out,ntop,tool_run=seurat_run,labels){
  
  
  for(y in seq(0.1,1,0.1)){
    
    score=sapply(ntop,function(x) tool_scores_by_specificity(sim,specificity = y,tool_out,x,seurat_run,labels))
    if(y==0.1){
      s=cbind(specificity=rep(y,length(ntop)),ntop,score)
    }else{
      s=rbind(s,cbind(specificity=rep(y,length(ntop)),ntop,score))}
    
    
  }
  
return(s)
}