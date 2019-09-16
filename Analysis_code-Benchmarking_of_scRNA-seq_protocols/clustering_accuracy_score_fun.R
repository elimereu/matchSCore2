## data: seurat object with stored the labels of clusters (in data@active.ident) and identities from the classification (in data@meta.data$nnet2)
## tech: character vector with names of protocols. Lenght of tech corresponds to length of protocols.
## cell_types: character vector with names of cell_types to be consider
cell_type_separation <-  function(data,tech,cell_types){
  
  
  
  clus <- data@active.ident
  nnet <- factor(data@meta.data$nnet2)
  df <- table(clus,nnet)
  clus.map <- apply(df,1,function(x) colnames(df)[which(x==max(x))])
  
  
  cells.p <- which(data@meta.data$batch==tech)
  clus <- data@active.ident[cells.p]
  levels(clus) <- clus.map
  nnet <- factor(data@meta.data$nnet2[cells.p])
  
  
  df <- matrix(table(clus,nnet)[,cell_types],nrow = length(levels(clus)))
  colnames(df) <- cell_types
  rownames(df) <- rownames(table(clus,nnet)[,cell_types])
  
  acc <- apply(df,2,function(x) max(x)/sum(x))
  mapping <- apply(df,2,function(x) levels(clus)[which(x==max(x))])
  names(mapping) <- cell_types
  
  acc <- sapply(seq(1:length(cell_types)),function(x) ifelse(mapping[x]==cell_types[x],acc[x],sum(df[which(cell_types[x]==levels(clus)),x])/sum(df[,x])))
  ave_acc <- sum(acc)/length(levels(clus)) ## 
  
  return(ave_acc)  
}

