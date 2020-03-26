#' Rank simulated group markers
#'
#' This function rank group markers according their relative specificity based on fold changes. 
#' @param sim A Splatter sim object with simulated cellular groups or paths.
#' @return A data frame containing the group gene ranking and other additional information for each gene from the simulated data set.
#' @export
#' @examples
#' 


rank_sim = function(sim){
  
  pd=colData(sim)
  fd=rowData(sim)
  group=factor(pd$Group)
  l=length(levels(group))
  ng=length(fd$Gene)
  
  fc=fd[,grep("^DEFac",names(fd))]
  
  de=sapply(1:l,function(x) ifelse(fc[,x]>1,1,0))
  rank_info=sapply(1:l,function(x) ifelse(de[,x]*rowSums(de)>=1,rowSums(de),l*10))
  rank_info=data.frame(fc,rank_info)
  
  ord=sapply(1:l,function(x) order(rank_info[,(x+l)],-rank_info[,x]))
  rank=apply(ord,2,function(x) order(x))
  
  rank_df=data.frame(fd[,c(1:4)],r_info=rank_info,rank_info,Rank=rank)
  

  
  return(rank_df)
}