#' Plots matchSCore trends for a specific tool
#'
#' This function identifies true label groups between reference groups and clusters
#' @param scores A character vector of reference group labels.
#' @export
#' @examples

tscores_plots = function(scores_df){

  
  names(scores_df)[c(2,3)]=c("top_genes","matchSCore")
  
  p1 <- ggplot(scores_df, aes_string(x="specificity", y="matchSCore", colour="top_genes",group="top_genes")) +
    geom_line() + scale_colour_gradientn(colours = rainbow(15)) +theme_bw() 
  
  p2 <- ggplot(scores_df, aes(x=top_genes, y=matchSCore, colour=specificity)) +
    geom_point(alpha=0.7,size=2.5) + scale_colour_gradientn(colours = c("gold","gray","blue")) +
    geom_smooth(span=1.3,se=F) + theme_bw()
  
  p3 <- ggplot(scores_df, aes(specificity, matchSCore, fill = top_genes))+ geom_bar( stat="identity",position = "dodge",width=0.03,col="white") + 
    scale_fill_gradientn(colours = rainbow(15)) +  theme_bw()
  
  
  ggarrange(p1,p2,p3,ncol = 1,nrow = 3)
  
}