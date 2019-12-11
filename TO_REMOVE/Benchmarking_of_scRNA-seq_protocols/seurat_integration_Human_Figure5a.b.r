load(file="sce.all_classified.technologies.RData")
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
# 117.0   979.2  1441.0  1854.0  2267.8  6432.0 
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
data <- FindClusters(data,resolution = 0.1)
table(data@active.ident)
0    1    2    3    4    5    6    7 
5213 4201 3490 3255 1145  416  191  123 

load(file="~/Dropbox/HCA_benchmarking_sample/col.batch.RData")
load(file="~/Dropbox/HCA_benchmarking_sample/col_PBMC_human.RData")


DimPlot(data, reduction = "umap", group.by = "batch",cols = col.batch)
DimPlot(data, reduction = "umap", group.by = "nnet2",cols = col) 
DimPlot(data, reduction = "umap") 

save(data,file="data_seu.obj_dim16.RData")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col2,pt.size = 1)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap



png(file="UMAP_merged_by_nnet2.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()


umap1<-DimPlot(object = data.integrated, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()


clus <- data@meta.data$integrated_snn_res.0.2
table(clus,data@meta.data$nnet2)

t <- table(clus,data@meta.data$nnet2)
lev <- colnames(t)[apply(t,1,function(x) which(x==max(x)))]
lev
[1] "HEK cells"         "CD4 T cells"       "CD14+ Monocytes"   "HEK cells"         "CD8 T cells"      
[6] "NK cells"          "B cells"           "CD4 T cells"       "FCGR3A+ Monocytes" "B cells"          
[11] "B cells"           

levels(data@meta.data$nnet2)

col2 <- c("red","green4","aquamarine","brown1","maroon","deepskyblue3","blueviolet","darkolivegreen1","black","violet","pink")

umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "integrated_snn_res.0.2",no.axes = T,cols = col2)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap

png(file="UMAP_merged_by_clusters.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()

levels(data@meta.data$integrated_snn_res.0.2) <- c("0:HEK293T cell","1:CD4+ T-cell","2:CD14+ Monocyte","3:HEK293T cell",
                                                   "4:CD8+ T-cell","5:NK cell","6:B cell","7:CD4+ T cell","8:FCGR3A+ Monocyte","9:B cell","10:B cell")                                                   "9:B cell","10:B cell")


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "integrated_snn_res.0.2",no.axes = T,cols = col2)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap

png(file="UMAP_merged_by_clusters_leg.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()
