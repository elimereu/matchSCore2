setwd("Chromium_V3/")
dir.create("out_Seurat/")


load(file="sce_hsap2_without.viability.RData")
genes=read.table("~/Dropbox/annotation.old/human/gencode_28.tsv",header = T)

pd <- colData(sce_hsap2)
pd$Library <- factor(pd$Library)
counts <- counts(sce_hsap2)
dim(counts)

cd1=data.frame(id=rownames(counts),counts)
genes_cd=merge(genes,cd1,by="id",all.x=F)

rownames(genes_cd)=make.unique(as.character(genes_cd$name))
counts=genes_cd[,6:dim(genes_cd)[2]]

sce_hsap2@assays$data$counts <- as.matrix(counts)

library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)



y <- Matrix::colSums(cd.sparse[,]>0)

q <- quantile(y,probs=0.05) 

data <- CreateSeuratObject(counts = cd.sparse, min.cells = 5, min.features = q, project = "ChromiumV3.Without.Viability")
dim(data@assays$RNA@data)
# 25161  4094

mito.genes <- grep(pattern = "^MT", x = rownames(x = data@assays$RNA@data), value = TRUE)
ribo.genes <- grep(pattern = "^RB", x = rownames(x = data@assays$RNA@data), value = TRUE)
data <- AddMetaData(data,metadata = data.frame(pd),col.name = colnames(pd))
data[["percent.mt"]] <- PercentageFeatureSet(data, features = mito.genes)
data[["percent.rb"]] <- PercentageFeatureSet(data, features = ribo.genes)




VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.rb"),ncol = 3)
VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),ncol = 3)


data <- NormalizeData(object = data,scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)


data <- ScaleData(object = data,vars.to.regress = c("nCount_RNA","percent.mito"))
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:14)
data <- FindClusters(data, resolution = 0.2)
table(data@active.ident)
0    1    2    3    4    5    6 
1487  809  795  394  355  219   35 


data <- RunUMAP(data, dims = 1:14)
data <- RunTSNE(data, dims = 1:14)

DimPlot(data, reduction = "umap")
DimPlot(data, reduction = "tsne")


library(ggplot2)


col=c("blueviolet","aquamarine","green4","orange","black","maroon","coral2","deepskyblue3","red","gold","gray")



tsne1<-DimPlot(data, reduction = "tsne")
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()


markers=FindAllMarkers(data,test.use = "wilcox",only.pos = T)

head(markers)
library(dplyr)
library(matchSCore2)
top10 <-markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10

gene_cl <- cut_markers(levels(markers$cluster),markers,ntop = 100)
gene_cl

png(file="out_Seurat/plots/TSNE_Known_markers.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features = c("CD4","MS4A1","CD79A","CD79B","GNLY", "CD3E","CD3D" , 
                                        "FCGR3A","CD14","LYZ", "NKG7","CD8A"), cols = c("lightgray", "blue"), 
            reduction = "tsne")
dev.off()


VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.rb"),ncol = 3)
VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),ncol = 3)


levels(data@active.ident) <- c("CD4 T cells","CD8 T cells and NK","CD14+ Monocytes","HEK cells","B cells","low-count cells","unclear")


levels(markers$cluster) <- levels(data@active.ident)
names(gene_cl) <- levels(data@active.ident)

save(markers,file="out_Seurat/markers.RData")
save(gene_cl,file="out_Seurat/gene_cl.RData")
save(data,file="out_Seurat/data_seu.obj_res0.2_dim14.RData")

tsne1<-DimPlot(data, reduction = "tsne")
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

png(file="out_Seurat/plots/TSNE_clusters_annotated.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()




