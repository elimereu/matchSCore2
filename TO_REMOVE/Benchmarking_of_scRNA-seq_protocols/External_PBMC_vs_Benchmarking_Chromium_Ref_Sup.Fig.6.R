load(file = "data_seu.obj_res0.6_annotated.RData") ### external PBMC dataset
load(file="pbmc_gene_cl.RData") ### external PBMC gene signatures

library(ggplot2)
library(Seurat)

pbmc <- data
rm(data)
pbmc <- UpdateSeuratObject(pbmc)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(data, reduction = "umap")


load(file="col_PBMC_human.RData")
col2 <- unname(col)
umap1<-DimPlot(object = pbmc, reduction.use = 'umap')
umap<-umap1+ theme_void()+scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_clusters_PBMC-external_dataset.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()

save(pbmc,file="pbmc_external_ref.RData")

library(matchSCore2)
clusters <- pbmc@active.ident
scaled <- pbmc@assays$RNA@scale.data
mod <- train_model(scale.data = scaled,clus = clusters,gene_cl.ref = pbmc.gene_cl)

save(mod, file="mod.Rdata")


load(file="10X/Human/out_Seurat/seu/data_dim9_res0.1.RData") #### HCA reference dataset (Chromium)

scale.data <- data@scale.data
ids <- identity_map(scale.data = data@scale.data,model = mod,gene_cl.ref = pbmc.gene_cl)


annotation <- factor(ids$ids)


load(file="../reference_mod/reference3_last_version.RData") ### reference cluster annotations from the benchmarking PBMC+ HEK human data (Chromium)

not_hek <- names(reference3)[grep("^HEK",reference3,invert = T)] ## rm hek from the comparison

DimPlot(data, reduction = "umap",cells = not_hek)

data@meta.data$PBMC_ext_anno <- annotation

umap1<-DimPlot(object = data, reduction.use = 'umap',cells = not_hek,group.by = "PBMC_ext_anno")
umap<-umap1+ theme_void()+scale_color_manual(values=col)+ theme(legend.text = element_text(size = 14))+ theme(plot.margin = unit(c(0.3,1,1,0), "lines"))


png(file="UMAP_clusters_PBMC-HCAref_dataset.png",width = 8, height = 6,units = 'in', res = 600)
umap
dev.off()

load(file="gene_cl.ref3.RData") ### gene signatures from the benchmarking PBMC+ HEK human data (Chromium)
p <- which(names(gene_cl.ref)=="HEK cells") ## remove HEK gene signature

gene_cl.ref <- gene_cl.ref[-p]

png(file="MatchSCore_plot.png",width = 8, height = 6,units = 'in', res = 600)
matchSCore2(gene_cl.ref = pbmc.gene_cl,gene_cl.obs = gene_cl.ref,ylab = "External PBMC data",xlab = "PBMC Benchmarking Sample")
dev.off()

names(annotation) <- colnames(scale.data)
my_ids <- factor(annotation[not_hek])
my_clus <- factor(reference3[not_hek])

gg <- summary_barplot(class.fac = my_ids,obs.fac = my_clus)

gg

load(file="col_PBMC_human.RData")
col2 <- unname(col)
col2 <- c(col,"unclassified"="lightgray")
gg+scale_color_manual(values = col2)

library(matchSCore2)

png(file="Barplots.png",width = 8, height = 6,units = 'in', res = 600)
summary_barplot(class.fac = my_ids,obs.fac = my_clus)+scale_fill_manual(values = col2)
dev.off()
