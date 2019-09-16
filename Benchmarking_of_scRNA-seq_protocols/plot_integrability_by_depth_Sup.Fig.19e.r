load(file = "Integration/Seurat/cell_type_sep_vs_mixedness.RData") ### df: data.frame with two columns. Column 1: clustering accuracy after integration across protocols. Column 2: mixability. Data were integrated with Seurat after downsampling to 20K 
ds20k <- data.frame(df,method=rownames(df),depth=rep("20K",dim(df)[1]))

load(file="Integration/Seurat.DS10K/cell_type_sep_vs_mixedness.RData") ### df: data.frame with two columns. Column 1. clustering accuracy after integration. Column 2: mixability. Data were integrated with Seurat after downsampling to 10K 
ds10k <- data.frame(df,method=rownames(df),depth=rep("10K",dim(df)[1]))

load(file="col.batch.RData")

df <- rbind(ds10k,ds20k)

library(ggplot2)


png(file="cell_type_sep_vs_mixedness_10K_and_20K_depth.png",width = 6, height = 6,units = 'in', res = 600)
ggplot(df,aes(x=acc,y=mix))+geom_point(aes(colour=method,size=depth))+labs(colour="Method")+theme_bw()+
  xlim(c(0.5,1))+ylim(c(0.5,0.95))+ 
  geom_line(aes(colour=method),linetype = 2)+
  xlab("Clustering Accuracy")+ylab("Mixability")+ scale_color_manual(values = col.batch)+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1),
        axis.text = element_text(size=20),text = element_text(size = 20)) 

dev.off()

png(file="cell_type_sep_vs_mixedness_10K_and_20K_depth_nolegend.png",width = 6, height = 6,units = 'in', res = 600)
ggplot(df,aes(x=acc,y=mix))+geom_point(aes(colour=method,size=depth))+labs(colour="Method")+theme_bw()+
  xlim(c(0.5,1))+ylim(c(0.5,0.9))+ geom_line(aes(colour=method),linetype = 2)+
  xlab("Clustering Accuracy")+ylab("Mixability")+ scale_color_manual(values = col.batch)+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1),
        axis.text = element_text(size=20),text = element_text(size = 20),legend.position = "none") 

dev.off()


png(file="Mixability_scores_10K_and_20K_depth_nolegend.png",width = 6, height = 6,units = 'in', res = 600)
ggplot(df,aes(x=method,y=mix,fill=depth))+geom_bar(aes(fill=depth),stat="identity",position=position_dodge())+theme_bw()+
  xlab("")+ylab("Mixability")+ scale_fill_brewer(palette = "red")+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.line = element_line(size=1),
        axis.text.x = element_text(size=20,angle = 65,hjust = 1,vjust = 1),axis.text.y = element_text(size=20),text = element_text(size = 20),legend.position = "top") 

dev.off()
