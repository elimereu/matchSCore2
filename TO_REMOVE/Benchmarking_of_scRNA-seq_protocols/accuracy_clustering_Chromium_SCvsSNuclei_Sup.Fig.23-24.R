load(file="acc_by_ct.harmony.RData") ### clustering accuracy by harmony
load(file="acc_by_ct.seurat.RData") ### clustering accuracy by seurat
load(file="acc_by_ct.scmerge.RData") ### clustering accuracy by scMerge

acc_by_ct.seurat$CellType <- rownames(acc_by_ct.seurat)
acc_by_ct.harmony$CellType <- rownames(acc_by_ct.harmony)
acc_by_ct.scmerge$CellType <- rownames(acc_by_ct.scmerge)

acc <- rbind(acc_by_ct.scmerge,acc_by_ct.seurat,acc_by_ct.harmony)
acc
library(reshape2)
acc.m <- melt(acc) 

library(ggplot2)
load(file="col.batch.RData")
df <- acc.m


df$CellType <- factor(df$CellType)


png(file="Accuracy_clustering_scVSsn2.png",width = 6, height = 4,units = 'in', res = 600)
ggplot(df, aes(x=Tool,y=`Clustering Accuracy`))+geom_boxplot(aes(colour=Protocol),outlier.shape = NA)+ 
  geom_jitter(aes(colour=Protocol,shape=CellType),size=3)+
  scale_color_manual(values = col.batch)+theme_bw()+
  theme(axis.text.y= element_text(colour = "black",size=16),
        axis.title.y = element_text(size=16),axis.title.x = element_blank(),axis.line = element_line(size=1),
        axis.text.x = element_text(size=16,angle = 45,hjust = 1,vjust = 1),
        legend.text = element_text(size=14),
        legend.title= element_text(size = 14))
dev.off()


