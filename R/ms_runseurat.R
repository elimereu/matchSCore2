#' This function runs Seurat3 by using the output of the `ms_integrate` function.
#'
#' @param out_align The output of the function `ms_integrate`. The combined count
#' matrix is used to create the Seurat object and the integrated is provided
#' to the slot `@data`. # TODOELI: we would need to say a little more on how this
#' list is expected to look like
#' @param dims Seurat parameter. It is the dimension of the PCA space.
#' @param res Seurat resolution parameter.
#' @param col_anno TODO
#' @param col_data TODO
#' @param verbose Logical, controls the displaying of additional messages while
#' running the function. Defaults to `TRUE`.
#'
#' @return The integrated Seurat object. Two UMAP plots with cells coloured by
#' annotation and dataset will be generated.
#'
#' @export
#'
#' @examples
#' # TODO
ms_runseurat <- function(out_align,
                        dims = c(1:10),
                        res = 0.2,
                        col_anno = NULL,
                        col_data = NULL,
                        verbose = TRUE) {
  counts <- out_align$counts
  integrated <- out_align$integrated
  annotation <- out_align$annotation_label
  dataset <- out_align$dataset_label

  if (verbose) message("Running Seurat by using the integrated matrix as normalized data...")

  counts <- counts[, colnames(integrated)]
  data <- CreateSeuratObject(counts = counts, min.features = 0, min.cells = 0, project = "integrated")

  data <- AddMetaData(object = data, metadata = annotation, col.name = "cluster")
  data <- AddMetaData(object = data, metadata = dataset, col.name = "dataset")

  data <- NormalizeData(object = data)
  VariableFeatures(data) <- rownames(integrated)

  data@assays$RNA@scale.data <- integrated  # TODOELI: is there a Seurat specific function for this?

  data <- RunPCA(data, features = VariableFeatures(object = data))
  plot(ElbowPlot(data))
  data <- FindNeighbors(data, dims = dims)
  data <- FindClusters(data, resolution = res)
  data <- RunUMAP(data, dims = dims)

  return(data)
}
