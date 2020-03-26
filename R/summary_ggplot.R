# TODO: to be exported - or anyway documented?
summary_ggplot <- function(data,ylab,xlab){

  # library(ggplot2)
  # library(reshape2) # TODO: really needed?
  # library(grid) # TODO: really needed?

  my_df <- data.frame(t(data),check.names = F,check.rows = F)
  my_df.melt <-  melt(cbind(x=1:nrow(my_df),my_df),id ="x")

  gg <- ggplot(my_df.melt, aes(x=factor(x),y=variable,fill=value)) + labs(x=xlab,y=ylab)+
    geom_tile(aes(fill = value)) + scale_x_discrete(lab=rownames(my_df))+
    theme(axis.text.x=element_text(angle=30,hjust = 1,size=16),axis.text.y = element_text(size=16),axis.title = element_text(size=16))+
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient(low = "white", high = "red",name="Jaccard Index")

  return(gg)
}

