#' Compute matchSCore for a specified value of specificity and ntop ranked genes
#'
#' This function computes the matchSCore for different cutoffs of marker specificity and top cluster markers
#' @param sim A Splatter simulation object.
#' @param specificity The proportion of top-ranked markers for each simulated group.
#' @param tool_out The output of the used tool to cluster your data (It should be as the output of the seurat_run function).
#' @param ntop The proportion of top-ranked markers for each cluster.
#' @param tool_run The name of the function used to run the clusters without quote (Default= seurat_run).
#' @param labels Cluster labels as in the output of the compute_labels function.
#' @return A matchSCore value
#' @export
#' @examples
#' 


tool_scores_by_specificity = function(sim,specificity,tool_out,ntop,tool_run=seurat_run,labels){
  
  fd=rowData(sim)
  pd=colData(sim)
  groups=factor(pd$Group)
  lev=levels(groups)
  k=length(lev)
    
 
  gene_cl=do.call(tool_run,args=list(sim,ntop,tool_out))
  
  rank_df =rank_sim(sim)
  markers=markers_by_specificity(rank_df,specificity,k)
  
  
  score=matchSCore(markers,gene_cl,labels)
  
  
  return(score)
  
}