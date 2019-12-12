#' This function runs Seurat3 by using the output of the `align_run` function.
#'
#' @param counts The output of the function `align_run`. The combined count
#' matrix used to create the new integrated Seurat object.
#' @param integrated The output of the function `align_run`. It will be used as
#' the `scale.data` of the new integrated Seurat object.
#' @param annotation A vector with the cell annotation like the `annotation_label`
#' output of the align_run function.
#' @param dataset a vector indicating the dataset of origin for each cell like
#' the `dataset_label` output of the `align_run` function.
#' @param dims Seurat parameter. It is the dimension of the PCA space.
#' @param res Seurat resolution parameter.
#'
#' @return The integrated Seurat object. Two UMAP plots with cells coloured by
#' annotation and dataset will be generated.
#'
#' @export
#'
#' @examples
#' # TODO
seurat3_run <- function(counts,integrated,annotation,dataset,dims=c(1:10),res=0.2,col_anno=NULL,col_data=NULL){

  require(Seurat)
  require(cowplot)
  require(ggplot2)

  print("Running Seurat by using as scale.data the integrated matrix..")

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



  DimPlot(object = data, reduction.use = 'umap',group.by = "cluster",no.axes = TRUE,cols = col_anno)+theme_void()+
         theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

  ggsave("plot1.png")

  DimPlot(object = data, reduction.use = 'umap',group.by = "dataset",no.axes = TRUE,cols = col_data)+theme_void()+
    theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

  ggsave("plot2.png")
  return(data)
}

