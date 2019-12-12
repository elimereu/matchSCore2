#' Integration of different datasets
#'
#' @param dataset_list A named list of annotated Seurat objects. The cell
#' annotations are in the `active.ident` slots.
#' @param marker_list A list of cell-type specific markers from the input dataset.
#'
#' @return A list containing:
#' - `integrated`: the integrated scaled matrix.
#' - `counts`: the combined count matrix.
#' - `annotation_label`: the cell annotations.
#' - `dataset_label`: a vector indicating the dataset of origin for each cell.
#' - `genes`: the set of common genes used.
#'
#' @export
#'
#' @examples
#' # TODO
align_run <- function(dataset_list,marker_list){
  # require(Seurat)
  # require(Matrix)
  print("Defining the set of common genes")

  total <- 10
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  Sys.sleep(0.1)
  # update progress bar
  progress <- 1
  setTxtProgressBar(pb,progress)

  start.time <- Sys.time()


    genes <-unique(unlist(marker_list))

    bg <- lapply(dataset_list,function(x) intersect(rownames(x),genes))
    l <- length(bg)
    t <- table(unlist(bg))
    genes <- names(t)[which(t==l)]

  progress <- 2
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)



  dataset_list <- lapply(dataset_list,function(x) ScaleData(x,features = genes,verbose = FALSE))
  d_list <- lapply(dataset_list,function(x) x@assays$RNA@scale.data[genes,])
  cl <- unique(unlist(lapply(dataset_list,function(x) factor(x@active.ident))))

  progress <- 3
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  cc <- lapply(dataset_list,function(x) sapply(levels(cl),function(y) rowMeans(x@assays$RNA@scale.data[genes,which(x@active.ident==y)]) ))
  cc <- lapply(cc,function(x) apply(x,1,function(y) median(y,na.rm = TRUE)))

  progress <- 4
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  H <- lapply(seq(1,l),function(x) d_list[[x]]-cc[[x]])
  cov <- lapply(seq(2,l), function(x) cov(H[[1]],H[[x]]))
  svd_out <- lapply(cov, function(x) svd(x))

  progress <- 5
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  Ropt <- lapply(svd_out, function(x) x$v %*% t(x$u))
  new_coord <- lapply(Ropt,function(x) x %*% t(H[[1]]))


  progress <- 6
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)

  for(i in seq(1,l)){
    if(i>1){
      integrated <- cbind(integrated,t(new_coord[[i-1]]))
    }else{
      integrated <- H[[1]]
    }

  }

  cells <- unlist(lapply(dataset_list,function(x) names(x@active.ident)))
  colnames(integrated) <- cells

  counts_list <- lapply(dataset_list,function(x) x@assays$RNA@counts[genes,names(x@active.ident)])
  #md <- lapply(dataset_list,function(x) x@meta.data[names(x@active.ident),])

  progress <- 7
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  for(i in seq(1,l)){
     if(i>1){
       counts <- cbind(counts,counts_list[[i]])
     }else{
       counts <- counts_list[[1]]
     }}

  progress <- 8
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)


  #require(plyr)

  # for(i in seq(1,l)){
  #
  #   if(i>1){
  #     meta.data <- rbind.fill(meta.data,md[[i]])
  #   }else{
  #     meta.data <- md[[1]]
  #   }}


    # meta.data$clus_annotation <- unlist(lapply(dataset_list,function(x) factor(x@active.ident)))
    # rownames(meta.data) <- colnames(integrated)
  annotation <- unlist(lapply(dataset_list,function(x) as.character(x@active.ident)))

  progress <- 9
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)

  dataset <- unlist(lapply(seq(1:length(dataset_list)),function(x) rep(names(dataset_list)[x],ncol(dataset_list[[x]]))))

  progress <- 10
  Sys.sleep(0.1)
  setTxtProgressBar(pb,progress)
  close(pb)
  end.time <- Sys.time()
  time <- difftime(end.time,start.time,units="mins")
  print(paste("The runtime is:",time,"min",sep=" "))

  return(list(counts=counts,integrated=integrated,annotation_label=annotation,dataset_label=dataset,genes=genes))

}
