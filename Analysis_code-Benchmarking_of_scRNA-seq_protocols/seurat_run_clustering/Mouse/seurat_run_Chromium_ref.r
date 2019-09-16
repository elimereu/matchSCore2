dir.create("out_Seurat/")


load(file="mmus_filt.RData")

pd <- colData(mmus_filt)
counts <- counts(mmus_filt)


library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)
rm(counts)
summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8     556    1113    1980    3094    9020 
plot(density(Matrix::colSums(cd.sparse[,]>0)))
 

dim(cd.sparse)
#36580 32237


y <- Matrix::colSums(cd.sparse[,]>0)
summary(y)

q <- quantile(y,probs=0.05) 
q
135
data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 5, min.genes = q, project = "Chromium-ref")
rm(cd.sparse)
dim(data@data)
#28925 30603

mito.genes <- grep(pattern = "^Mt", x = rownames(x = data@data), value = TRUE)
percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)

ribo.genes <- grep(pattern = "^RB|RP", x = rownames(x = data@data), value = TRUE)
percent.ribo <- Matrix::colSums(data@raw.data[ribo.genes, ])/Matrix::colSums(data@raw.data)


# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
data <- AddMetaData(object = data, metadata = pd, col.name = colnames(pd))
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
data <- AddMetaData(object = data, metadata = percent.ribo, col.name = "percent.ribo")

VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"),group.by="Library", nCol = 3,x.lab.rot = T)
VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,x.lab.rot = T)



data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.07, 
                          x.high.cutoff = 3, y.cutoff = 0.5)

data <- ScaleData(object = data,vars.to.regress = c("nUMI"))

data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:16, genes.print = 5)
data <- ProjectPCA(object = data, do.print = F)

PCElbowPlot(object = data)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:15, resolution = 0.1, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
table(data@ident)
0     1     2     3     4     5     6     7     8     9    10 
11266  4379  3473  2978  1799  1706  1407  1013   885   849   848 


data <- RunTSNE(object = data, dims.use = 1:15, do.fast = TRUE,force.recalc=T)
data <- RunUMAP(object = data,dims.use = 1:15)
library(ggplot2)


col=c("aquamarine1","coral1","darkolivegreen1","darkmagenta","gray","magenta","rosybrown1","sienna4","tan1","steelblue1","seagreen4")


levels(data@ident) <- c("Enterocytes 1","NIH1","Enteroendocrine cells","Enterocytes 2","Amplifying cells","NIH2","NIH3","unclear","NIH5","NIH6","Immuno cells")

save(data,file="out_Seurat/data_dim15_res0.1.RData")

tsne1<-TSNEPlot(object = data)
tsne<-tsne1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters_annotated_onlycolon.png",width = 10, height = 8,units = 'in', res = 300)
tsne
dev.off()

o <- sapply(levels(data@ident),function(x) mean(data@meta.data$nGenes[which(data@ident==x)]))
# 0         1         2         3         4         5         6         7         8         9        10 
# 1122.5668 4795.3580 2088.7961 1047.6504  219.2617 4774.4183 1142.8127 3861.8479  747.4995 4775.5671 1079.2636 
sort(o)


png(file="out_Seurat/plots/TSNE_clusters_with_manual_annotations.png",width = 8, height = 6,units = 'in', res = 600)
tsne
dev.off()


umap1<-DimPlot(object = data, reduction.use = 'umap')
umap<-umap1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

png(file="out_Seurat/plots/UMAP_clusters_anno.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()


markers=FindAllMarkers(data,test.use = "wilcox",only.pos = T)


png(file="UMAP_Known_markers_onlycolon.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data,features.plot = c("Slc26a3", "Lgr5","Smoc2","Reg4","Atoh1","Muc2","Tff3","Mki67","Pcna","Chga","Chgb"), cols.use = c("lightgray", "blue"), 
            reduction.use = "umap")
dev.off()


png(file="out_Seurat/plots/TSNE_Known_markers2_onlycolon.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data,cells.use = names(colon.clus),features.plot = c("Polr1a","Polr1b","Slc26a3", "Ceacam5","Lgr5","Smoc2","Reg4","Atoh1","Muc2","Tff3","Mki67","Foxm1","Pcna","Chga","Chgb","Lamc2","Fabp1"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()

png(file="out_Seurat/plots/TSNE_Known_ISC.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Sparc","Col1a2","Col1a1","Col3a1","Lgr5","Smoc2","Mki67","Dcn","Epcam"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()


png(file="out_Seurat/plots/TSNE_Known_markers3.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Lyz1","Chgb","Chga","Alpi","Pax4","Neurog2","Krt7","Ang4","Reg4","Lars2","Sox4","Pcna","Ccnb1"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()

png(file="out_Seurat/plots/TSNE_Known_markers.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Tff3","Krt7","Muc2","Agr2","Pax4","Chga","Chgb","Tac1","Tph1","Afp","Ucn3",
                                             "Vgf","Ang4","Lyz1","Cd24a","Krt8","Krt18","Dclk1"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()


png(file="out_Seurat/plots/TSNE_Known_markers5.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Muc3","Atoh1","Ceacam1","Ceacam2","Fabp2"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()

