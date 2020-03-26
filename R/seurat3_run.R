#' This function runs Seurat3 by using the output of the `align_run` function.
#'
#' @param out_align The output of the function align_run. The combined count matrix is used
#' to create the Seurat object and the integrated is provided to the slot @data.
#' @param dims Seurat parameter. It is the dimension of the PCA space.
#' @param res Seurat resolution parameter.
#' @param col_anno TODO
#' @param col_data TODO
#'
#' @return The integrated Seurat object. Two UMAP plots with cells coloured by
#' annotation and dataset will be generated.
#'
#' @export
#'
#' @examples
#' # TODO


seurat3_run <- function(out_align,dims=c(1:10),res=0.2,col_anno=NULL,col_data=NULL){

  # require(Seurat)
  # require(cowplot)
  # require(ggplot2)

  counts <- out_align$counts
  integrated <- out_align$integrated
  annotation <- out_align$annotation_label
  dataset <- out_align$dataset_label

  message("Running Seurat by using as normalized data the integrated matrix..")

  counts <- counts[,colnames(integrated)]
  data <- CreateSeuratObject(counts = counts, min.features = 0, min.cells = 0,project = "integrated")
  data@meta.data$cluster <- annotation
  data@meta.data$dataset <- dataset
  data <- NormalizeData(object = data)
  VariableFeatures(data) <- rownames(integrated)

  data@assays$RNA@scale.data <- integrated


  data <- RunPCA(data, features = VariableFeatures(object = data))
  plot(ElbowPlot(data))
  data <- FindNeighbors(data, dims = dims)
  data <- FindClusters(data, resolution = res)
  data <- RunUMAP(data, dims = dims)

  return(data)
}

