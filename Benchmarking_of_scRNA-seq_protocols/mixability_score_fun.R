## data: seurat object with stored the labels of clusters (in data@active.ident) and identities from the classification (in data@meta.data$nnet2)
## cell_types: character vector with names of cell_types to be considered

mixing_score <- function(data,cell_type){
  
  print(cell_type)
  nnet <- factor(data@meta.data$nnet2)
  names(nnet) <- rownames(data@meta.data)
  batch <- factor(data@meta.data$batch)
  names(batch) <- names(nnet)
  cells <- names(nnet)[which(nnet %in% cell_type)]
  
  nnet <- factor(nnet[cells])
  batch <- factor(batch[cells])
  
  
  t <- table(nnet,batch)
  t <- matrix(t,nrow =length(levels(nnet)))
  sizes <- apply(t,1,function(x) min(x))
  names(sizes) <- levels(nnet)
  
  sub.tech <- sapply(levels(batch),function(x) cells[which(batch==x)])
  sub.nnet <- sapply(levels(nnet),function(x) sapply(sub.tech,function(y) sample(x = y[which(nnet[y]==x)],size = sizes[which(names(sizes)==x)])))
  samples <- unlist(sub.nnet)
  
  
  e <- as.matrix(data@reductions$umap@cell.embeddings)
  exp <- e[samples,]
  
  
  batch <- batch[samples]
  nnet <- nnet[samples]
  design <- data.frame(replicates=batch,condition=nnet)
  
  
  library(kBET)
  
  batch.estimate <- kBET(exp, batch= batch,verbose = T,heuristic = F,testSize = 50)
  index.out <- batch.estimate$outsider$index 
  x <- seq(1:length(batch))
  outsiders <- sapply(x,function(y) ifelse(y %in% index.out,TRUE,FALSE))
  info_test <- data.frame(design,res=batch.estimate$results,outsiders)
  p <- lapply(levels(batch),function(x) info_test$res.kBET.pvalue.test[which(batch==x)])
  names(p) <- levels(batch)
  
  acceptance_rate <- sapply(p,function(x) length(which(x>0.05))/length(x))
  return(acceptance_rate)
}

