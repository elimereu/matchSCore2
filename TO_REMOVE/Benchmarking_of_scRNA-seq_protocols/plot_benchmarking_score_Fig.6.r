
load(file="ave.gene.detection.Rev.RData") 
load(file="clus.accuracy.RData") ## ave. clus. accuracy DS20K 
load(file="ave.mappability.RData") ## classification DS20K 
load(file="ave.expression.markers.RData") ## average expression for markers in HEK,B, Monocytes 
load(file="~/Dropbox/HCA_benchmarking_sample/Revision/Integration/AVEtools_integrability.RData") ### integratability scores (clus accuracy and mixability)


bench_scores <-data.frame(gene_detection=gene_detection,
                          marker_exp=ave.markers[names(gene_detection)],
                          clus_accuracy=clus_accuracy[names(gene_detection)],
                          mappability=ave.mappability[names(gene_detection)],
                          cell_type_merging=acc_by_ct[names(gene_detection)],
                          mixedness_merging=scores[names(gene_detection)])

library(psych)
bscores <- apply(bench_scores,1,function(x) harmonic.mean(x))
bench_scores <- data.frame(method=rownames(bench_scores),bench_scores,score=bscores)


bench_scores <- bench_scores[order(bench_scores$score,decreasing = T),]
library(reshape2)
bench_scores.df.m <- melt(bench_scores) 


library(ggplot2)

df <- bench_scores.df.m
gg <-ggplot(df, aes(x=variable,y=method)) + geom_point(aes(colour=value,size=value))+scale_color_viridis_c(option="plasma")+
  theme(axis.text.y= element_text(colour = "black",size=16),axis.title = element_blank(),axis.line = element_blank(),
        axis.ticks = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(panel.background = element_rect(fill = "white",colour = "white", linetype = "solid"))

gg$data$method <- factor(gg$data$method,levels=names(sort(bscores)))


png(file="Benchmarking_score.png",width = 5, height = 5, units = 'in', res = 600)
ggplot(gg$data[grep("score",gg$data$variable,invert = T),], aes(x=variable,y=method)) + geom_point(aes(colour=value,size=value))+scale_color_viridis_c(option="plasma")+
  theme(axis.text.y= element_text(colour = "black",size=18),axis.title = element_blank(),axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(panel.background = element_rect(fill = "white",colour = "white", linetype = "solid"))
dev.off()                                                                            
            