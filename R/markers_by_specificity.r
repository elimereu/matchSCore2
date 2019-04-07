#' Simulated group markers 
#'
#' This function gives you the ranked list of group markers at the specified proportion of top markers (specificity)
#' @param rank_df The data frame with the ranked group genes as returned by the rank_sim function. 
#' @param spec The proportion of top ranked genes. It has to be a number between 0 and 1. 
#' @param n_clust The number of simulated cell groups. 
#' @return A list with n_clust elements representing their corresponding group markers. 
#' Each element of the list contains the relative set of marker indexes as ordered in the original rowData(sim). 
#' The list doesn't keep the order of the ranking. 
#' @export
#' @examples
#' 
markers_by_specificity = function(rank_df,spec,n_clust){
  
  k=5+n_clust
  info=rank_df[,c(k:(k+n_clust-1))] ## order
  p=grep("^Rank",names(rank_df))
  R=rank_df[,p] ## rank
  
  sorted=sapply(1:n_clust,function(x) order(R[,x]))
  k=10*n_clust
  length=sapply(1:n_clust,function(x) length(which(info[,x]<k)))
  
  l=sapply(1:n_clust,function(x) round(spec*length[x]))
  markers=lapply(1:n_clust,function(x) sorted[1:l[x],x])
  return(markers)
}
