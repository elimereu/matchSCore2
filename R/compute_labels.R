#' Clustering matching 
#'
#' This function identifies true label groups between reference groups and clusters
#' @param idc A character vector of reference group labels.
#' @param idc1 Vector of labels for clusters.
#' @return Vector of corresponding reference labels for clusters. 
#' @export
#' @examples
#' 



compute_labels=function(idc,idc1){

  idc=factor(idc)
  idc1=factor(idc1)
  lev=levels(idc)
  lev1=levels(idc1)

for(i in 1:length(lev)){
  p1=which(idc==i)
  cl=as.numeric(idc1[p1])
  x=as.numeric(names(sort(table(cl),decreasing=TRUE)[1])) 
  if(i==1){lab=x} else{lab=append(lab,values = x,after = length(lab))}
}
  
  return(lab)
}


