setwd("Integration/Seurat/")
load(file="data_seu.obj.RData")
load(file="col_PBMC_human.RData")
load(file="col.batch.RData")
load(file="gene_cl.ref.RData") ### gene signatures

tech <- unique(data@meta.data$batch) ## protocols
data@meta.data$batch <- as.factor(data@meta.data$batch) ### cell protocol

nnet <- factor(data@meta.data$nnet2) ## matchSCore identity
cell_types <- levels(nnet)[c(1:7,9)] ## cell types of interest
cell_types




exp.data <- as.matrix(data@assays$RNA@data) ## expression
exp.markers <- sapply(gene_cl.ref,function(x) Matrix::colMeans(exp.data[which(rownames(exp.data) %in% x),])) ### average gene signature x cell 

exp.markers <- data.frame(exp.markers)
colnames(exp.markers)

exp.markers <- exp.markers[,c(1:7,9)] ## average expression for cell type of interest

cells <- rownames(data@meta.data)[which(nnet %in% cell_types)]
cells.p <-which(nnet %in% cell_types)
df <- data.frame(exp.markers[cells,],id=factor(data@meta.data$nnet2[cells.p]),batch=factor(data@meta.data$batch[cells.p]))
head(df)
table(df$id)

table(df$batch)
df$cell <- rownames(df)
head(df)

df.x <-lapply(levels(df$id),function(x) df[which(df$id==x),])
names(df.x) <- levels(df$id)



library(ggplot2)

gg1 <- ggplot(df.x$`B cell`,aes(x=cell,y=B.cells))+geom_bar(aes(fill=batch),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=20),axis.line = element_line(size=2))+
  xlab("B cells")+ylab("B cell signature")


p1 <- which(df.x$`CD14+ Monocyte`$CD14..Monocytes>0)
gg2 <- ggplot(df.x$`CD14+ Monocyte`[p1,],aes(x=cell[p1],y=CD14..Monocytes[p1]))+geom_bar(aes(fill=batch[p1]),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=20),axis.line = element_line(size=2))+
  xlab("CD14+ Monocytes")+ylab("CD14+ Monocyte signature") 


gg3 <- ggplot(df.x$`HEK cell`,aes(x=cell,y=HEK.cells))+geom_bar(aes(fill=batch),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=20),axis.line = element_line(size=2))+
  xlab("HEK293T cells")+ylab("HEK293T signature")



gg4 <- ggplot(df.x$`CD4 T cell`,aes(x=cell,y=CD4.T.cells))+geom_bar(aes(fill=batch),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=18),axis.line = element_line(size=2))+
  xlab("CD4+ T-cells")+ylab("CD4+ T-cell signature") 

gg5 <- ggplot(df.x$`CD8 T cells`,aes(x=cell,y=CD8.T.cells))+geom_bar(aes(fill=batch),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=18),axis.line = element_line(size=2))+
  xlab("CD8+ T-cells")+ylab("CD8+ T-cell signature") 

gg6 <- ggplot(df.x$`NK cells`,aes(x=cell,y=NK.cells))+geom_bar(aes(fill=batch),stat = "identity",position = "dodge")+
  scale_fill_manual(values=col.batch)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=24),
                                            axis.title = element_text(size=18),axis.line = element_line(size=2))+
  xlab("NK cells")+ylab("NK cell signature") 


library(ggpubr)

png(file="barplot_markers_signatures_by_protocol_HEK_B_Monocytes.png",width = 18, height = 6, units = 'in', res = 600)
ggarrange(gg3,gg1,gg2,ncol = 3,nrow = 1,common.legend = T,legend ="none")+
  theme(legend.position = "none")
dev.off()
