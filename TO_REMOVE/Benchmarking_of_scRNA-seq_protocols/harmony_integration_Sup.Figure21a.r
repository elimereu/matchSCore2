load(file="Human/DS_20K/merged/sce.all_classified.technologies.RData")

library(Seurat)
library(harmony)




pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
counts <- sce@assays$data$counts


library(ggplot2)
library(dplyr)
library(harmony)
data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data <- NormalizeData(data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
data <- RunHarmony(data, group.by.vars = "batch")
data <- RunUMAP(data, reduction = "harmony", dims = 1:16)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:16) %>% FindClusters()
table(data@active.ident)
0    1    2    3    4    5    6    7    8    9   10   11   12 
2822 2643 2439 2335 1308 1239 1195 1194 1142  878  425  291  123 

DimPlot(data, group.by = c("batch", "nnet2"), ncol = 2)



load(file="~/Dropbox/HCA_benchmarking_sample/col.batch.RData")
load(file="~/Dropbox/HCA_benchmarking_sample/col_PBMC_human.RData")

DimPlot(data, reduction = "umap", group.by = "batch",cols = col.batch)
DimPlot(data, reduction = "umap", group.by = "nnet2",cols = col) 


save(data,file="data_seu.obj.RData")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_nnet2.png",width = 3, height = 3,units = 'in', res = 600)
umap+theme(legend.position = "none")
dev.off()

png(file="UMAP_merged_by_nnet2_leg.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()

umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 3, height = 3,units = 'in', res = 600)
umap+theme(legend.position = "none")
dev.off()

png(file="UMAP_merged_by_protocols_leg.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data, reduction.use = 'umap',no.axes = T)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap



