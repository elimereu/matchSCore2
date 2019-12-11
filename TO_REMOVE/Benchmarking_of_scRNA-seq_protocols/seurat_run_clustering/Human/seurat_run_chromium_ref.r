dir.create("out_Seurat/")

load(file="hsap_filt.RData")

pd <- colData(hsap_filt)
counts <- counts(hsap_filt)

### Filtering out cells with %MT > 25%
mito.genes <- grep(pattern = "^MT-", x = rownames(x = counts), value = TRUE)
percent.mito <- colSums(counts[mito.genes, ])/colSums(counts)
high.mito <- colnames(counts)[which(percent.mito>0.25)]
length(high.mito)
#2691
counts <- counts[,!(colnames(counts) %in% high.mito)]


library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)
rm(counts)
summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8     645     758    1503    1020   12192 
plot(density(Matrix::colSums(cd.sparse[,]>0)))
 

dim(cd.sparse)
#39689 40234


q <- quantile(y,probs=0.05) 
q
508
data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 5, min.genes = q, project = "Chromium-ref")
rm(cd.sparse)
dim(data@data)
#32060 38195


# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
data <- AddMetaData(object = data, metadata = pd, col.name = colnames(pd))


data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2,x.high.cutoff = 3, y.cutoff = 0.5)

data <- ScaleData(object = data,vars.to.regress = c("nUMI","percent.mito"))

data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:10, genes.print = 5)
data <- ProjectPCA(object = data, do.print = F)

PCElbowPlot(object = data)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:9, resolution = 0.1, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
table(data@ident)
0     1     2     3     4     5     6     7     8 
13654  6285  5376  5375  3449  1858  1083   670   445 



data <- RunTSNE(object = data, dims.use = 1:9, do.fast = TRUE,force.recalc=T)
data <- RunUMAP(object = data,dims.use = 1:9)
library(ggplot2)


tsne1<-TSNEPlot(object = data)
tsne<-tsne1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne
dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters_annotated.png",width = 10, height = 8,units = 'in', res = 300)
tsne
dev.off()

umap1<-DimPlot(object = data, reduction.use = 'umap')

umap<-umap1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/UMAP_clusters_anno.png",width = 10, height = 8,units = 'in', res = 300)
umap
dev.off()


markers=FindAllMarkers(data,test.use = "wilcox",only.pos = T)



gene_cl <- cut_markers(seq(1:length(levels(markers$cluster)))-1,markers,ntop = 100)



gene_cl <- lapply(gene_cl,function(x) x[!is.na(x)])
seq10x<- gene_cl
names(seq10x) <- levels(data@ident)
save(seq10x,file="out_Seurat/gene_cl_seq10x.RData")


dir.create("out_Seurat/plots")

png(file="out_Seurat/plots/TSNE_Known_markers_umap.png",width = 10, height = 8, units = 'in', res = 600)
FeaturePlot(object = data,features.plot = c("CD3E","CD3D","IL7R","CD4","CD8A","CD8B","MS4A1","CD79A","CD79B","GNLY","NKG7",  
                                              "CD14","LYZ","FCGR3A"), cols.use = c("lightgray", "blue"), 
            reduction.use = "umap")
dev.off()

levels(data@ident) <- c("CD4+ T cells1","CD8+ T cells and NK","HEK cells1","CD14+ and FCGR3A+ Monocytes","B cells","CD4+ T cells2","HEK cells2","HEK cells3","unclear")

