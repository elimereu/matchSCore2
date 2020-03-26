# TODO: to be exported - or anyway documented?
summary_ggplot <- function(data, ylab, xlab) {

  my_df <- data.frame(t(data), check.names = FALSE, check.rows = FALSE)
  my_df.melt <- melt(cbind(x = 1:nrow(my_df), my_df), id = "x")
  my_df.melt$x <- factor(my_df.melt$x)

  gg <- ggplot(my_df.melt, aes_string(x = "x", y = "variable", fill = "value")) +
    labs(x = xlab, y = ylab) +
    geom_tile(aes_string(fill = "value")) +
    scale_x_discrete(lab = rownames(my_df)) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16)) +
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient(low = "white", high = "red", name = "Jaccard Index")

  return(gg)
}
