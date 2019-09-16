dir.create("out_Seurat/")


load(file="mmus_filt.RData")

pd <- colData(mmus_filt)
counts <- counts(mmus_filt)

### Filtering out cells with %MT > 25%
mito.genes <- grep(pattern = "^Mt", x = rownames(x = counts), value = TRUE)
percent.mito <- colSums(counts[mito.genes, ])/colSums(counts)
high.mito <- colnames(counts)[which(percent.mito>0.25)]
length(high.mito)
#0
counts <- counts[,!(colnames(counts) %in% high.mito)]


library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)
rm(counts)
summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#13.0   343.5   686.0  1087.5  1375.0  5792.0 
plot(density(Matrix::colSums(cd.sparse[,]>0)))
 

dim(cd.sparse)
#27347   391


y <- Matrix::colSums(cd.sparse[,]>0)

q <- quantile(y,probs=0.05) 

data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 5, min.genes = q, project = "MARS-Seq")
rm(cd.sparse)
dim(data@data)
#13768   371

mito.genes <- grep(pattern = "^Mt", x = rownames(x = data@data), value = TRUE)
percent.mito <- Matrix::colSums(data@raw.data[mito.genes, ])/Matrix::colSums(data@raw.data)
ribo.genes <- grep(pattern = "^RP|RB", x = rownames(x = data@data), value = TRUE)
percent.ribo <- Matrix::colSums(data@raw.data[ribo.genes, ])/Matrix::colSums(data@raw.data)
# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
data <- AddMetaData(object = data, metadata = pd, col.name = colnames(pd))
data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
data <- AddMetaData(object = data, metadata = percent.ribo, col.name = "percent.ribo")

VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"),group.by="Library", nCol = 3,x.lab.rot = T)
VlnPlot(object = data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,x.lab.rot = T)



data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.45, 
                          x.high.cutoff = 3, y.cutoff = 0.5)

data <- ScaleData(object = data,vars.to.regress = c("nUMI"))

data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:16, genes.print = 5)
data <- ProjectPCA(object = data, do.print = F)

PCElbowPlot(object = data)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
table(data@ident)
0   1   2   3   4   5 
107 106  55  47  38  18  

data <- RunTSNE(object = data, dims.use = 1:10, do.fast = TRUE,force.recalc=T)
data <- RunUMAP(object = data,dims.use = 1:10)
library(ggplot2)


col=c("aquamarine1","coral1","darkolivegreen1","darkmagenta","gray","magenta","rosybrown1","sienna4","tan1","steelblue1","seagreen4")

tsne1<-TSNEPlot(object = data)
tsne<-tsne1+ ggtitle("Clusters") +scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne
dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()


png(file="out_Seurat/plots/TSNE_Known_markers2.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Polr1a","Polr1b","Slc26a3", "Ceacam1","Lgr5","Smoc2","Reg4","Muc2","Agr2","Mki67","Foxm1","Pcna","Lamc2","Fabp2"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()


markers=FindAllMarkers(data,test.use = "wilcox",only.pos = T)



png(file="out_Seurat/plots/TSNE_Known_markers3.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Lyz1","Alpi","Krt7","Ang4","Reg4","Lars2","Sox4","Pcna","Ccnb1"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()

png(file="out_Seurat/plots/TSNE_Known_markers.png",width = 12, height = 12, units = 'in', res = 600)
FeaturePlot(object = data, features.plot = c("Tff3","Krt7","Muc2","Agr2",
                                             "Ang4","Lyz1","Cd24a","Krt8","Krt18","Dclk1"), cols.use = c("lightgray", "blue"), 
            reduction.use = "tsne")
dev.off()



