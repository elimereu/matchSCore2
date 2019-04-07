#' Markers and differentially expressed genes of the simulated data set
#'
#' This function gives you information about differentially expressed genes, markers and and markers specific of only
#' one group. The definition is based according to the DEFac value provided by the Splatter simulation. 
#' Here a gene is defined marker of a group if its DEFac>1 in that group.
#' @param sim A Splatter sim object with simulated cellular groups or paths.
#' @return A list with four outputs: a data frame (de) indicating differentially expressed (DE) genes (1) and not (0). 
#' Columns represent groups and rows are genes.
#' The list with sets of de genes per group (gr_de).
#' The list with sets of group markers (markers). 
#' The list with sets of markers specific only in one group (sp_markers).
#' 
#' @export
#' @examples
#' 

defac_matrix = function(sim){
  
  pd=colData(sim)
  fd=rowData(sim)
  group=factor(pd$Group)
  k=length(levels(group))
  
  fc=fd[,grep("^DEFac",names(fd))]
  
  l=length(levels(group))
  de=sapply(1:l,function(x) ifelse(fc[,x]!=1,1,0))
  
  clust_set = apply(de,2,function(x) which(x==1))
  
  de=sapply(1:l,function(x) ifelse(fc[,x]>1,1,0))
  markers = apply(de,2,function(x) which(x==1))
  
  dd=apply(fc,1,function(x) ifelse(length(which(x!=1))==1,which(x>1),0))
  sp_markers=lapply(1:k,function(x) which(dd==x))
  
  return(list(de=de,gr_de=clust_set,markers=markers,sp_markers=sp_markers))
  
}
