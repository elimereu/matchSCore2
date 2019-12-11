
load(file="Human/DS_10K/sce.all_classified.technologies.RData")
load(file="col.batch.RData")
load(file="col_PBMC_human.RData")

pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
counts <- sce@assays$data$counts

library(Seurat)
library(ggplot2)

data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data.list <- SplitObject(data, split.by = "batch")

summary(colSums(counts[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 82     759    1117    1382    1740    4397
for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


reference.list <- data.list
data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:16)
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:16)


library(cowplot)

DefaultAssay(data.integrated) <- "integrated"

data.integrated <- ScaleData(data.integrated, verbose = T)
data.integrated <- RunPCA(data.integrated, npcs = 16, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:16)
data <- FindNeighbors(data.integrated, dims = 1:16)
data <- FindClusters(data,resolution = 0.2)
table(data@active.ident)
0    1    2    3    4    5    6    7    8    9   10   11 
4670 4255 2882 2081 1874 1262 1169  528  511  232  214  127 


DimPlot(data, reduction = "umap", group.by = "batch",cols = col.batch)
DimPlot(data, reduction = "umap", group.by = "nnet2",cols = col) 
DimPlot(data, reduction = "umap") 

save(data,file="data_seu.RData")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_nnet2.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()

png(file="UMAP_merged_by_nnet2_leg.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data.integrated, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()

png(file="UMAP_merged_by_protocols_leg.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()


