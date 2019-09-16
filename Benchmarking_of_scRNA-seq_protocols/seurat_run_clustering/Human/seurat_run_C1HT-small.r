dir.create("out_Seurat/")

load(file="hsap_filt.RData")

pd <- colData(hsap_filt)
counts <- counts(hsap_filt)

### Filtering out cells with %MT > 25%
mito.genes <- grep(pattern = "^MT-", x = rownames(x = counts), value = TRUE)
percent.mito <- colSums(counts[mito.genes, ])/colSums(counts)
high.mito <- colnames(counts)[which(percent.mito>0.25)]
counts <- counts[,!(colnames(counts) %in% high.mito)]


library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)

summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2407    4672    5910    7511   10786   18444
plot(density(Matrix::colSums(cd.sparse[,]>0)))
 

dim(cd.sparse)
#41790  2380


y <- Matrix::colSums(cd.sparse[,]>0)

q <- quantile(y,probs=0.05) 

data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 5, min.genes = q, project = "C1HT-small")
dim(data@data)
#34285  2261


# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
data <- AddMetaData(object = data, metadata = pd, col.name = colnames(pd))
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")


VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,x.lab.rot = T)

data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.08, 
                          x.high.cutoff = 3, y.cutoff = 0.5)

data <- ScaleData(object = data,vars.to.regress = c("nUMI","percent.mito"))
data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:10, genes.print = 5)
data <- ProjectPCA(object = data, do.print = F)
PCElbowPlot(object = data)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:6, resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
table(data@ident)
0   1   2   3   4   5   6 
622 302 250 144 143  86  85 

data <- RunTSNE(object = data, dims.use = 1:6, do.fast = TRUE,force.recalc=T)

library(ggplot2)

col=c("blueviolet","aquamarine","green4","orange","black","maroon","coral2","deepskyblue3","red","gold","gray")
save(data,file="out_Seurat/data_seu.obj_res0.3_dim8.RData")

tsne1<-TSNEPlot(object = data)
tsne<-tsne1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()


markers=FindAllMarkers(data,test.use = "wilcox",only.pos = T)

head(markers)
library(dplyr)
top10 <-markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10


png(file="out_Seurat/plots/Heatmap_cluster_Top10markers.png",width = 12, height = 12, units = 'in', res = 600)
DoHeatmap(object =data, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "blue",col.mid = "white", col.high = "red",cex.row = 8)
dev.off()

library(matchSCore2)
gene_cl <- cut_markers(levels(data@ident),markers,ntop = 100)


png(file="out_Seurat/plots/TSNE_Known_markers.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("CD4","MS4A1","CD79A","CD79B","GNLY", "CD3E","CD3D" , 
                                              "FCGR3A","CD14","LYZ", "NKG7","CD8A"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()


levels(data@ident) <- c("CD8+ T cells and NK","unclear","CD4+ cells","CD14+ and FCGR3A+ monocytes","HEK cells 1","B cells","HEK cells 2")
