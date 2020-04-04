#' Heatmaps of class probabilities
#'
#' This function identifies true label groups between reference groups and clusters.
#'
#' @param out output of the function `ms_identity_map`. # TODOELI: here as well, some param name more informative of what needs to be provided?
#'
#' @return A ggplot object, with the heatmap displaying the probability of each
#' cell to belong to the identity classes according the model.
#'
#' @export
#' @examples
#' # TODO
ms_identity_heatmap <- function(out) {
  ids.ord <- order(out$ids)
  ord.p <- out$fit.prob[ids.ord, ]

  my_df <- data.frame(ord.p, check.names = FALSE, check.rows = FALSE)
  my_df.melt <- melt(cbind(x = 1:nrow(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x)

  gg <- ggplot(my_df.melt, aes_string(x = "x", y = "variable", fill = "value")) +
    labs(x = "Cells", y = "Cell identity") +
    geom_tile(aes_string(fill = "value")) +
    scale_fill_gradientn(
      colours = c("gray100", "blue", "yellow", "deeppink"),
      name = "Probability\n"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "black", size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )

  return(gg)
}
