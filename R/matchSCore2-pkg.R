#' matchSCore2
#'
#' `matchSCore2` is a Bioconductor package that provides ... # TODOELI, we need a nice short para for this
#'
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom methods new
#' @importFrom stats cor median na.omit predict
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom  Matrix colSums
#' @importFrom topGO GenTable getSigGroups GOKSTest printGraph annFUN.org
#' annFUN.GO2genes score
#' @importFrom ggplot2 aes aes_string element_blank element_line element_text
#' geom_bar geom_line geom_point geom_smooth geom_text geom_tile ggplot
#' labs scale_colour_gradientn scale_fill_gradient scale_fill_gradientn
#' scale_x_discrete theme theme_bw
#' @importFrom org.Hs.eg.db org.Hs.egALIAS2EG org.Hs.egGO2ALLEGS
#' @importFrom org.Mm.eg.db org.Mm.egALIAS2EG org.Mm.egGO2ALLEGS
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @importFrom corrplot corrplot
#' @importFrom nnet multinom
#' @importFrom cowplot plot_grid
#' @importFrom S4Vectors DataFrame
#' @importFrom Seurat CreateSeuratObject ElbowPlot FindClusters FindNeighbors
#' NormalizeData RunPCA RunUMAP VariableFeatures `VariableFeatures<-`
#' AddMetaData
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' logcounts
#' @importFrom SummarizedExperiment `assay<-` `colData<-`
#' @importFrom corpcor fast.svd
#'
#' @name matchSCore2-pkg
#' @docType package
NULL

globalVariables(c("value"))
