
load(file="Mouse/DS_20/sce.all_classified.technologies.RData")
load(file="col.batch.RData")
load(file="Mouse/col_figures_with_Epr.RData")

pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
table(pd$batch)

counts <- sce@assays$data$counts

library(Seurat)
library(ggplot2)

data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data.list <- SplitObject(data, split.by = "batch")


for (i in c(1:length(data.list))) {
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


reference.list <- data.list
data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:6,k.filter = 10,k.score = 10)
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:6)


library(cowplot)

DefaultAssay(data.integrated) <- "integrated"

data.integrated <- ScaleData(data.integrated, verbose = T)
data.integrated <- RunPCA(data.integrated, npcs = 10, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:6)
data <- FindNeighbors(data.integrated, dims = 1:6)
data <- FindClusters(data,resolution = 0.2)


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



umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap

clus <- data@meta.data$integrated_snn_res.0.2
table(clus,data@meta.data$nnet2)

t <- table(clus,data@meta.data$nnet2)
lev <- colnames(t)[apply(t,1,function(x) which(x==max(x)))]
lev
[1] "Transit Amplifying" "Transit Amplifying" "Enterocytes 2"      "Enterocytes 1"      "Fibroblast"        
[6] "Secretory cells"    "Fibroblast"         "Transit Amplifying" "Immuno cells"       "Immuno cells"     


library(ggplot2)
umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "integrated_snn_res.0.2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap


png(file="UMAP_merged_by_clusters.png",width = 3, height = 3,units = 'in', res = 600)
umap+ theme(legend.position = "none")
dev.off()

col2 <- unname(col2) 
levels(data@meta.data$integrated_snn_res.0.2) <- c("0:Transit Amplifying","1:Transit Amplifying","2:Enterocytes 2",
                                                   "3:Enterocytes 1",
                                                   "4:Fibroblast","5:Secretory cell","6:Fibroblast","7:Transit Amplifying",
                                                   "8:Immuno cells",
                                                   "9:Immuno cells")                                                   


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "integrated_snn_res.0.2",no.axes = T,cols = col2)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
umap
png(file="UMAP_merged_by_clusters_leg.png",width = 5, height = 3,units = 'in', res = 600)
umap
dev.off()


