#' matchSCore2
#'
#' `matchSCore2` is a Bioconductor package that provides ...
#'
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom methods new
#' @importFrom stats cor median na.omit predict
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import Matrix
#' @importFrom topGO GenTable getSigGroups GOKSTest printGraph annFUN.org
#' annFUN.GO2genes score
#' @import ggplot2
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @importFrom corrplot corrplot
#' @importFrom nnet multinom
#' @import igraph
#' @import cowplot
#' @import Seurat
#' @importFrom SingleCellExperiment colData rowData
#'
#' @name matchSCore2-pkg
#' @docType package
NULL

# cowplot although here we use ggarrange, TODO? (is it from ggpubr)
