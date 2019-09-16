setwd("Chromium_V3/")
dir.create("out_Seurat/")


load(file="sce_mmus2_without.viability.RData")
genes=read.table("~/Dropbox/annotation.old/mouse/gencode_M17.tsv.gz",header = T)

pd <- colData(sce_mmus2)
pd$Library <- factor(pd$Library)
counts <- counts(sce_mmus2)
dim(counts)

cd1=data.frame(id=rownames(counts),counts)
genes_cd=merge(genes,cd1,by="id",all.x=F)

rownames(genes_cd)=make.unique(as.character(genes_cd$name))
counts=genes_cd[,6:dim(genes_cd)[2]]

sce_mmus2@assays$data$counts <- as.matrix(counts)

library(Seurat)
library(Matrix)
cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)

summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#128     590    1762    2698    4444   10752
plot(density(Matrix::colSums(cd.sparse[,]>0)))
 

dim(cd.sparse)
#33605  4469


y <- Matrix::colSums(cd.sparse[,]>0)

q <- quantile(y,probs=0.05) 

data <- CreateSeuratObject(counts = cd.sparse, min.cells = 5, min.features = q, project = "ChromiumV3.Without.Viability")
dim(data@assays$RNA@data)
# 23713  4245

mito.genes <- grep(pattern = "^Mt", x = rownames(x = data@assays$RNA@data), value = TRUE)
ribo.genes <- grep(pattern = "^Rb|Rp", x = rownames(x = data@assays$RNA@data), value = TRUE)
data <- AddMetaData(data,metadata = data.frame(pd),col.name = colnames(pd))
data[["percent.mt"]] <- PercentageFeatureSet(data, features = mito.genes)
data[["percent.rb"]] <- PercentageFeatureSet(data, features = ribo.genes)




VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.rb"),ncol = 3)
VlnPlot(object = data, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),ncol = 3)


data <- NormalizeData(object = data,scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)


data <- ScaleData(object = data,vars.to.regress = c("nCount_RNA"))
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:14)
data <- FindClusters(data, resolution = 0.1)
table(data@active.ident)
0    1    2    3    4    5    6    7    8    9   10   11 
1219  777  722  425  418  336   99   91   75   36   24   23 

data <- RunUMAP(data, dims = 1:14)
data <- RunTSNE(data, dims = 1:14)

DimPlot(data, reduction = "umap")
DimPlot(data, reduction = "tsne")


library(ggplot2)


col=c("aquamarine1","coral1","darkolivegreen1","darkmagenta","gray","magenta","rosybrown1","sienna4","tan1","steelblue1","seagreen4","black","red")

tsne1<-DimPlot(data, reduction = "tsne",pt.size = 2.5)
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters.png",width = 8, height = 8,units = 'in', res = 600)
tsne
dev.off()


load(file="~/Dropbox/HCA_benchmarking_sample/Mouse/col_figures_with_Epr.RData")
load(file="~/Dropbox/HCA_benchmarking_sample/Mouse/10X/reference_mod/gene_cl_seq10x_AllStates.RData")
load(file="~/Dropbox/HCA_benchmarking_sample/Mouse/10X/reference_mod/mod_referenceHCA_10X_acc0.87.RData")

library(matchSCore2)
library(nnet)

data <- ScaleData(object = data,features = unlist(seq10x),vars.to.regress = c("nCount_RNA"))

scaled <- data@assays$RNA@scale.data

ids <- identity_map(scale.data = scaled,model = mod,gene_cl.ref = seq10x)
nnet2 <- ids$ids

data@meta.data$nnet2 <- nnet2
save(data,file="out_Seurat/data_seu.obj.RData")

tsne1<-DimPlot(data, reduction = "tsne",group.by="nnet2")
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters_annotated.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()


tsne1<-DimPlot(data, reduction = "tsne",group.by="nnet2",pt.size = 2.5)
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

dir.create("out_Seurat/plots")
png(file="out_Seurat/plots/TSNE_clusters_annotated.png",width = 8, height = 8,units = 'in', res = 600)
tsne
dev.off()

tsne1<-DimPlot(data, reduction = "umap",group.by="nnet2",pt.size = 2.5)
tsne<-tsne1+ scale_color_manual(values=col)+theme_void()+theme(legend.text = element_text(size = 14)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))
tsne

png(file="out_Seurat/plots/UMAP_clusters_annotated.png",width = 8, height = 8,units = 'in', res = 300)
tsne
dev.off()

pd <- data.frame(colData(sce_mmus2))
names(nnet2) <- colnames(data)
head(nnet2)
head(pd)
names <- sapply(rownames(pd),function(x) paste("X",x,sep=""))
annotation <- nnet2[which(names(nnet2) %in% names)]
head(annotation)
rownames(pd) <- names
pd2 <- pd[names,]
head(pd2)

text <- element_text(size = 24)
mapped <- rowSums(pd2[,c(2:4)])
Nreads <- pd2[,1]
pd2$PCTmapped <- mapped/Nreads
pd2$Mapped <- mapped
df <- data.frame(Mapped=mapped[names(annotation)],nGenes=pd2[names(annotation),"nGenes"],annotation)

png(file="Mapped_vs_Ngenes_10K_shape_by_highMito_Mouse_by_annotations.png",width = 8, height = 6, units = 'in', res = 600)
ggplot(df,aes(x=Mapped,y=nGenes))+geom_point(aes(colour=annotation),size=1.5,alpha=0.8)+theme_bw()+scale_color_manual(values = col)+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+scale_y_log10()+labs(y = "nGenes (log-scale)",x="Mapped (log-scale)")+scale_x_log10()+
  theme(panel.grid = element_blank(),axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()

