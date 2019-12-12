#' This function computes the correlation between cells of the same cellular type
#' that are from different batches (e.g. protocols).
#'
#' @param raw A combined matrix of counts with gene expressions from all batches.
#' Rows are genes and columns are cells.
#' @param nnet A named factor with the annotation per cell.
#' @param batch A named factor with the batch label per cell.
#' @param cell_types A character with names (keywords for nnet levels) for cell
#' types you want to compute the correlation (e.g. "B|HEK|Monocytes").
#' @param n Number of cells sampled for the computation of the expression average
#' for cells that are from the same type. Default=minimum number of cells across
#' all batches.
#' @param genes A set of genes you want to use to compute the correlation.
#'
#' @return A matrix with correlations between batches.
#'
#' @export
#'
#' @examples
#' # TODO
cor_batch <- function(raw,nnet,cell_types,batch,n=NULL,genes=NULL){


  id <- factor(nnet[grep(cell_types,nnet)])
  cells <- names(id)

  if(is.null(genes)){
    raw <- raw[,cells]
  }else{

    raw <- raw[which(rownames(raw) %in% genes),cells]

  }


  batch <- factor(batch[cells])
  t <- table(id,batch)
  t
  sizes <- apply(t,1,function(x) min(x))


  sub.tech <- sapply(levels(batch),function(x) cells[which(batch==x)])
  sub.id <- sapply(levels(id),function(x) sapply(sub.tech,function(y) sample(x=y[which(id[y]==x)],size = sizes[which(names(sizes)==x)])))

  merge <- lapply(sub.id,function(x) apply(x,2,function(y) rowSums(raw[,y])))
  corr <- lapply(merge,function(x) cor(x,use="pairwise.complete.obs",method="pearson"))

  library(corrplot)
  col <- colorRampPalette(c(rep("white",3),"#FFECB3","#E85285","#6A1B9A"))

  lapply(corr,function(x) corrplot(x, method = "square",order = "hclust",hclust.method = "ward.D2",type="upper",tl.col="black",tl.cex=2,number.cex = 20))


  return(corr)
}
