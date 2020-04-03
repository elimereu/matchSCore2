#' Select the top n cluster markers
#'
#' This function identifies true label groups between reference groups and clusters.
#'
#' @param clusters A vector of cluster labels.
#' @param markers A `data.frame` as in the output of the `FindAllMarkers` Seurat
#' function.
#' @param ntop The number of top markers you want.
#'
#' @return A list of ntop markers per cluster.
#'
#' @export
#'
#' @examples
#' # TODO
cut_markers <- function(clusters,
                        markers,
                        ntop) {
  gene_cl <- lapply(
    clusters,
    function(x) markers$gene[markers$cluster == x][1:ntop]
  )
  gene_cl <- lapply(
    gene_cl,
    function(x) x[!is.na(x)]
  )
  names(gene_cl) <- clusters

  return(gene_cl)
}
