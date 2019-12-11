
library(harmony)
library(Seurat)


load(file="../scMerge/sce.all_classified.technologies.RData") 

pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
counts <- sce@assays$data$counts


library(ggplot2)
data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data <- NormalizeData(data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
data <- RunHarmony(data, group.by.vars = "batch")
data <- RunUMAP(data, reduction = "harmony", dims = 1:10)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution=0.5)


DimPlot(data, group.by = c("batch", "nnet2"), ncol = 2)



load(file="col.batch.RData")
load(file="col_PBMC_human.RData")



save(data,file="data_seu.obj.RData")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_nnet2.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data, reduction.use = 'umap',no.axes = T)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap

