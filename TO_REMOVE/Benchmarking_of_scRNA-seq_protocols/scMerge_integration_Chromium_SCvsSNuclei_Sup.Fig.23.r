library(scater)
require(plyr)

files <- list.files(path="Human", pattern="*.Robj", full.names=T, recursive=FALSE)
files ##### seurat object per protocol (Downsampled to 20K)


p <- grep("Chromium|Nuclei10X",files)

tech <-sapply(files[p],function(x) unlist(strsplit(x,split = "/",fixed = T))[[8]])
tech
tech <-sapply(tech,function(x) unlist(strsplit(x,split = ".",fixed = T))[[1]])

files <- files[p]

############# Merging #######################
tech <- c("Chromium","Chromium (sn)")
sce_list <- list(chromium,chromium_sn)
names(sce_list) <- tech
metadata <- c("nnet2","ident") #### id clustering (ident) and matchSCore2 classification (nnet2)

library(scMerge)

sce<-sce_cbind(sce_list, method = "union", cut_off_batch = 0.001,cut_off_overall = 0.001, exprs = c("counts","logcounts"),
               colData_names = metadata, batch_names = names(sce_list))


sce <- subset(sce,i=1:dim(sce)[1],j=which(colData(sce)$nnet2!="unclassified"))


load(file="SEG_human.RData") ### stable set of genes provided by scMerge
seg_human
batch <- as.factor(sce@colData$batch)
table(batch)
id <- factor(sce@colData$nnet2)


table(batch)
# batch
# Chromium Chromium (sn) 
# 1599           856 


out <- scReplicate(sce,batch = batch, kmeansK = c(rep(6,2)),hvg_exprs="logcounts",exprs = "logcounts",return_all = TRUE,fast_svd = T) 
hvg <- out$HVG
save(hvg,file="HVG_from_merged")


sce <- scMerge(sce, ctl = seg_human, kmeansK = rep(6,length(tech)),hvg_exprs="logcounts",
            exprs = "logcounts",fast_svd = TRUE,return_all_RUV = FALSE,assay_name = "merged")

save(sce,file="sce_merged.RData")


batch <- as.factor(sce@colData$batch)

#### Seurat run after scMerge integration ######


pd <- data.frame(colData(sce))
counts <- sce@assays$data$counts

library(Seurat)


data <- CreateSeuratObject(counts = counts,min.cells = 5,project = "scRNA_vs_snRNA_integration")
dim(data@assays$RNA@data)
#24723  2455


data <- AddMetaData(data,metadata = pd,col.name = names(pd))


data@assays$RNA@data <- sce@assays$data$merged


data <- ScaleData(object = data,features = rownames(data))
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
table(data@active.ident)
0   1   2   3   4   5   6   7   8 
505 329 321 306 298 289 207 159  41 

md <- data@meta.data
table(md$seurat_clusters,md$nnet2)

data <- RunUMAP(data, dims = 1:10)
data <- RunTSNE(data, dims = 1:10)


load(file="col.batch.RData")
load(file="col_PBMC_human.RData")


library(ggplot2)



umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2",no.axes = T,cols = col)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_nnet2.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()


umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "batch",no.axes = T,cols = col.batch)
umap<-umap1+ theme_void()+theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_merged_by_protocols.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()

library(cowplot)
p1 <-DimPlot(data, reduction = "pca", group.by = "batch",cols = col.batch)
p2 <- DimPlot(data, reduction = "pca", group.by = "nnet2",cols = col) 

cowplot::plot_grid(p1,p2)




