
load(file="sce.all_classified.technologies.RData") 

pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
levels(pd$batch)[2] <- "Chromium(sn)"
counts <- sce@assays$data$counts

library(Seurat)
library(ggplot2)

data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data.list <- SplitObject(data, split.by = "batch")


for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


reference.list <- data.list[c("Chromium","Chromium(sn)")]
data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:10)

data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:10)


library(cowplot)

DefaultAssay(data.integrated) <- "integrated"

data.integrated <- ScaleData(data.integrated, verbose = T)
data.integrated <- RunPCA(data.integrated, npcs = 10, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:10)
data <- FindNeighbors(data.integrated, dims = 1:10)
data <- FindClusters(data,resolution = 0.5)


load(file="col.batch.RData")
load(file="col_PBMC_human.RData")

DimPlot(data, reduction = "umap", group.by = "batch",cols = col.batch)
DimPlot(data, reduction = "umap", group.by = "nnet2",cols = col) 
DimPlot(data, reduction = "umap") 

save(data,file="data_seu.obj.RData")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_nnet2.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data.integrated, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()


