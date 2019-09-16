#markers1: data.frame1 (from protocol1) of cluster-specific markers, as computed by the seurat function FindAllMarkers  
#markers2: data.frame2 (from protocol2) of cluster-specific markers, as computed by the seurat function FindAllMarkers  
jaccard_index <- function(markers1,markers2,ntop){
  
  library(matchSCore2)
  clus <- levels(factor(markers1$cluster))
  gene_cl1 <- cut_markers(clusters = clus,markers = markers1,ntop = ntop)
  gene_cl2 <- cut_markers(clusters = clus,markers = markers2,ntop = ntop)
  
  gene_cl1 <- lapply(gene_cl1, function(x) x[!is.na(x)])
  gene_cl2 <- lapply(gene_cl2, function(x) x[!is.na(x)])
  
  ji <- sapply(clus,function(x) length(intersect(gene_cl1[[x]],gene_cl2[[x]]))/(length(union(gene_cl1[[x]],gene_cl2[[x]]))))
  
  return(ji)
}
