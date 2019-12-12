#' Compute the cellular phenotype matching across simulated cell groups and computational clusters
#'
#' This function computes the matchSCore across gene markers provided from the simulation and 
#' cluster-specific genes computationally obtained from the tool. 
#' @param gene_markers A list of gene markers. Each element of the list represents the set of markers for each biological group.
#' @param gene_cl A list of cluster-specific gene sets to match with gene_markers.
#' @param labels A vector of gene_cl indexes correspondant to the gene_markers group level. 
#' @return A value from 0 to 1 representing the matchSCore between gene_cl and gene_markers. A matchSCore of 1 is the total match. 
#' @export
#' @examples
#' 

matchSCore <- function(gene_markers,gene_cl,labels){
  
  score=0
  
  for(i in 1:length(labels)){
    ind=labels[i]
    len1=length(gene_markers[[i]])
    len2=length(gene_cl[[ind]])
    I=length(intersect(gene_markers[[i]],gene_cl[[ind]]))
    J=I/(len1+len2-I)
    score=sum(score,J)
  }
  score=score/length(gene_markers)
  return(score)
}

