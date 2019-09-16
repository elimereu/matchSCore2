tech <- c("Chromium","inDrop","C1HT-medium","C1HT-small","CEL-Seq2","ddSEQ","Drop-Seq","ICELL8","MARS-Seq","Chromium(sn)","Quartz-Seq2","SCRB-Seq",
          "Smart-Seq2") ## protocols
sce_list <- list(chromium,indrop,c1ht.m,c1ht.s,celseq,ddseq,dropseq,icell8,marsseq,chromium_sn,quartzseq,scrbseq,smartseq)
names(sce_list) <- tech
metadata <- c("nnet2","ident") ## metadata includes matchSCore2 labels (nnet2) and cluster identities

library(scMerge)

sce<-sce_cbind(sce_list, method = "union", cut_off_batch = 0.001,cut_off_overall = 0.001, exprs = c("counts","logcounts"),
               colData_names = metadata, batch_names = names(sce_list))

sce <- subset(sce,i=1:dim(sce)[1],j=which(colData(sce)$nnet2!="unclassified"))

load(file="SEG_human.RData") ### set of stable genes provided by scMerge tool
seg_human
batch <- as.factor(sce@colData$batch)
id <- factor(sce@colData$nnet2)

sce <- scMerge(sce, ctl = seg_human, kmeansK = rep(6,length(tech)),hvg_exprs="logcounts",
               exprs = "logcounts",fast_svd = TRUE,return_all_RUV = FALSE,assay_name = "merged")

out <- scReplicate(sce,batch = batch, kmeansK = c(rep(6,length(tech))),hvg_exprs="logcounts",
                   exprs = "logcounts",fast_svd = TRUE,
                   return_all = TRUE)
hvg <- out$HVG  ### set of highly variable genes

#### Seurat run #####
library(Seurat)
library(Matrix)

dir.create("out_Seurat/")

pd <- colData(sce)
counts <- sce@assays$data$counts


cd.sparse=Matrix(data=round(as.matrix(counts),digits = 0),sparse=T)

summary(Matrix::colSums(cd.sparse[,]>0))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#117.0   979.2  1441.0  1854.0  2267.8  6432.0
plot(density(Matrix::colSums(cd.sparse[,]>0)))


data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 10, min.genes = 0, project = "HCA Human merged")
dim(data@data)
# 27565 18034
data@data <- sce@assays$data$merged #### The log normalized matrix comes from scMerge normalization
data@var.genes <-hvg ### Highly Variable genes are from scMerge
length(data@var.genes)
5431

data <- ScaleData(object = data)
data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:16, genes.print = 5)
data <- ProjectPCA(object = data, do.print = F)
PCElbowPlot(object = data)
data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:14, resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
table(data@ident)
0    1    2    3    4    5 
5701 3988 3211 2349 1666 1119 

data <- RunTSNE(object = data, dims.use = 1:14, do.fast = TRUE,force.recalc=T)

library(ggplot2)


col1 <- c("blue","red","gray","chartreuse","violet","gold","black","pink")
levels(data@ident) <- c("0:HEK293T cell","1:CD4+ T-cell","2:CD4+ T-cell","3:CD8+ T-cell and NK cell","4:CD14+ Monocyte","5:CD14+ and FCGR3A+ Monocyte","6:B cell","7:CD4+ T-cell")
data <- RunUMAP(object = data,dims.use = 1:14)
umap1<-DimPlot(object = data, reduction.use = 'umap')
umap<-umap1+scale_color_manual(values = col1) +ggtitle("Clusters") + theme(legend.text = element_text(size = 16),axis.text =element_text(size = 24),text = element_text(size = 18)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

png(file="UMAP_cluster.png",width = 9, height = 5,units = 'in', res = 600)
umap
dev.off()

load(file="col.batch.RData")

umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "batch")
umap<-umap1+ ggtitle("Clusters") + scale_color_manual(values=col.batch)+theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

png(file="UMAP_by_technology.png",width = 8, height = 6,units = 'in', res = 300)
umap
dev.off()

load(file="out_Seurat/col.clusters.RData")

umap1<-DimPlot(object = data, reduction.use = 'umap',group.by = "nnet2")
umap<-umap1+ ggtitle("Clusters") + scale_color_manual(values=col)+theme(legend.text = element_text(size = 14),axis.text =element_text(size = 16),text = element_text(size = 20)) + theme(plot.margin = unit(c(0.3,1,1,0), "lines"))

png(file="UMAP_by_cell_type.png",width = 8, height = 6,units = 'in', res = 300)
umap
dev.off()

