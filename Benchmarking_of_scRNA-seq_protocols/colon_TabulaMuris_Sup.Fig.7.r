setwd("TM") ## Tabula Muris (TM) colon data (Smart-seq2)
load(file="tiss.Robj") ### colon data from TM


t <- table(ident=tiss@ident,anno=tiss@meta.data$free_annotation)
lev <- apply(t,1,function(x) colnames(t)[which(x==max(x))]) 
levels(tiss@ident) <- lev
table(tiss@ident)

# lev
# c("Amplifying undifferentiated cell","Lgr5+ undifferentiated cell","Enterocyte (Proximal)","Lgr5+ undifferentiated cell",
# "Lgr5- undifferentiated cell","Goblet cell (Distal)","Enterocyte (Distal)","Lgr5+ undifferentiated cell","Goblet cell (Proximal)",
# "Lgr5+ undifferentiated cell","Goblet cell, top of crypt (Distal)","Tuft cell","Chromaffin Cell")


tiss <- UpdateSeuratObject(tiss)

tiss <- RunUMAP(tiss, dims = 1:20)


table(tiss@active.ident,tiss@meta.data$cell_ontology_class)

markers=FindAllMarkers(tiss,test.use = "wilcox",only.pos = T)
# 
head(markers)
library(dplyr)
top10 <-markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10
# 
save(markers,file="markers_TM.RData")



library(matchSCore2)

gene_cl <- cut_markers(levels(tiss@active.ident),markers,ntop = 100)
save(gene_cl,file="gene_cl_TM2.RData")

load(file="gene_cl_reference.mouse.RData") ### gene signatures from the Mouse reference HCA benchmarking data 
seq10x <- seq10x[which(!names(seq10x) %in% c("Immuno cells","Fibroblast"))]

png(file="MatchSCore_overlap_signatures2.png",width = 10, height = 6, units = 'in', res = 600)
matchSCore2(gene_cl.ref = gene_cl,gene_cl.obs = seq10x,ylab = "TM",xlab = "Colon benchmarking sample")
dev.off()

##### TM model ######

mod <- train_model(scale.data = tiss@assays$RNA@scale.data,clus = tiss@active.ident,gene_cl.ref = gene_cl,prop = 0.75,p.threshold = 0.5)
0.92



load(file="../../Mouse/10X/out_Seurat/seu/data_dim15_res0.1.RData") ### Mouse reference HCA benchmarking data 

ids <- identity_map(scale.data = data@scale.data,model = mod,gene_cl.ref = gene_cl)
table(ids$ids)

data@meta.data$TM_projections <- ids$ids
data <- UpdateSeuratObject(data)
load(file="../../Mouse/10X/reference_mod/reference_clus.RData")
load(file="../../Mouse/10X/col_figures_with_Epr.RData")

data@active.ident <- ref_anno
DimPlot(data, reduction = "umap",cells = names(ref_anno),cols = col)
names(ref_anno)

cells <- names(ref_anno)[grep("Immuno cells|Fibroblast",ref_anno,invert = T)]
png(file="UMAP_Colon_benchmark_TM_projections.png",width = 7, height = 6,units = 'in', res = 600)
DimPlot(data, reduction = "umap",cols = col_pools,group.by = "TM_projections",cells = cells)+theme_void()+theme(legend.position = "none")
dev.off()

png(file="UMAP_Colon_benchmark_TM_projections2.png",width = 8, height = 6,units = 'in', res = 600)
DimPlot(data, reduction = "umap",cols = col_pools,group.by = "TM_projections",cells = cells)+theme_void()+theme(legend.text = element_text(size = 14))
dev.off()




png(file="UMAP_anno_clusters_TM3.png",width = 7, height = 6,units = 'in', res = 600)
DimPlot(tiss, reduction = "umap",cols = col_pools)+theme_void()+theme(legend.position = "none")
dev.off()


markers=FindAllMarkers(tiss,test.use = "wilcox",only.pos = T)
# 
head(markers)
library(dplyr)
top10 <-markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top10
# 

load(file="metadata_seu.obj_colon_mouse.bench.ref.RData")
classification <- ids$ids
save(classification,file="classification_TM_matchSCore2_model.RData")
head(ref_anno)

head(cells)
t <- table(tm=classification[cells],ref=factor(ref_anno[cells]))
# ref
# tm                                   Transit Amplifying Enterocytes 1 Enteroendocrine Enterocytes progenitor Immuno cells Stem cells Enterocytes 2 Fibroblast
# Amplifying undifferentiated cell                  634             0               4                      0          432        471            45        548
# Chromaffin Cell                                   715             6             208                      3           17         26            96        104
# Enterocyte (Distal)                               515          1827               2                    271            4          0             1          6
# Enterocyte (Proximal)                             561           329              10                     43           10          3          2478         17
# Goblet cell (Distal)                              134             1              45                      4            3         14             1         28
# Goblet cell (Proximal)                            273             0              15                      1           28         26             3        187
# Goblet cell, top of crypt (Distal)                 56             9               2                      7            2          1             0          2
# Lgr5- undifferentiated cell                      1244             0              16                      0          120         17             6        292
# Lgr5+ undifferentiated cell                      1843           678              22                    572           78        998           176        118
# Tuft cell                                          78             0               2                      0          113         65             9         54
# unclassified                                      189            41               5                      6           41         10            15         51
# ref
# tm                                   Secretory cells
# Amplifying undifferentiated cell                68
# Chromaffin Cell                                 29
# Enterocyte (Distal)                              1
# Enterocyte (Proximal)                            2
# Goblet cell (Distal)                          1268
# Goblet cell (Proximal)                         630
# Goblet cell, top of crypt (Distal)              13
# Lgr5- undifferentiated cell                     10
# Lgr5+ undifferentiated cell                    573
# Tuft cell                                        6
# unclassified                                    62

ref=factor(ref_anno[cells])
table(ref)

n <- colSums(t)
n

percentage.cell_type <- sapply(names(n),function(x) round(t[,x]/n[x]*100,digits = 2))
colSums(percentage.cell_type)


df2 <- data.frame(percentage.cell_type)
head(df2)
library(reshape2)
y <- melt(df2)
y$tm <- df$tm
y$ref <- df$ref
y$Freq <- df$Freq
names(y)[2] <- "Cell.type.percent"

df <- data.frame(t)
head(df)
save(col_pools,file="col_pools_named.RData")
load(file="col_pools_named.RData")
png(file="Cumulative_bar_plots_cell.type_classifications.png",width = 7, height = 6,units = 'in', res = 600)
ggplot(y,aes(x=ref,y=Cell.type.percent)) + geom_bar(aes(y = Cell.type.percent, x = ref, fill = tm),stat="identity")+theme_bw()+
scale_fill_manual(values = col_pools)+
theme(axis.title = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
      axis.text.x = element_text(colour = "black",size = 14,angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(colour = "black",size = 16),
      axis.line = element_line(colour = "black",size = 1),panel.grid.major = element_blank(),panel.border = element_blank())
 # theme(plot.margin = unit(c(0.3,1,-1,1), "lines"))
dev.off()


t <- table(tm=classification[cells],ref=factor(ref_anno[cells]))
t
library(matchSCore2)

png(file="Cumulative_bar_plots_cell.type_classifications.png",width = 7, height = 6,units = 'in', res = 600)
summary_barplot(class.fac = classification[cells],obs.fac = factor(ref_anno[cells]))+scale_fill_manual(values = col_pools)+theme(axis.text.x = element_text(angle = 70))
dev.off()

