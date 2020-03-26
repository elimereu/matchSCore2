#' Cumulative barplot of cell mapped compositions
#'
#' This function identifies true label groups between reference groups and clusters.
#' @param class.fac A named vector of classified cells (e.g. the `output$ids`
#' from the [identity_map()] function).
#' @param obs.fac A named vector of clusters.
#'
#' @return A cumulative barplot describing the cell identity compositions of the
#' observed clusters onto the reference
#'
#' @export
#'
#' @examples
#' # TODO
summary_barplot <- function(class.fac,obs.fac){

  # library(ggplot2)
  # library(reshape2) # TODO: really needed?
  # library(grid) # TODO: really needed?

  t <- table(class.fac,obs.fac)
  df <- data.frame(t)
  n <- colSums(t)
  percentage.cell_type <- sapply(names(n),function(x) round(t[,x]/n[x]*100,digits = 2))
  df2 <- data.frame(percentage.cell_type)

  y <- melt(df2)
  y$class.fac <- df$class.fac
  y$obs.fac <- df$obs.fac
  y$Freq <- df$Freq
  names(y)[2] <- "Cell.type.percent"


  gg <- ggplot(y,aes(x=obs.fac,y=Cell.type.percent)) + geom_bar(aes(y = Cell.type.percent, x = obs.fac, fill = class.fac),stat="identity")+theme_bw()+
    theme(axis.title = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
          axis.text.x = element_text(colour = "black",size = 14,angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(colour = "black",size = 16),
          axis.line = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.border = element_blank())

          return(gg)
}
