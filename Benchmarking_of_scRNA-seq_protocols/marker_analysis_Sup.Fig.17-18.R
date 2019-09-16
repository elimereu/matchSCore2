load(file="Downsampled_analysis/Human/DS_20K/merged/sce.all_classified.technologies.RData")
load(file="col.batch.RData")
load(file="col_PBMC_human.RData")

pd <- data.frame(colData(sce))
pd$batch <-as.factor(pd$batch)
counts <- sce@assays$data$counts

library(Seurat)
library(ggplot2)

data <- CreateSeuratObject(counts,min.cells = 5,meta.data = pd)
data.list <- SplitObject(data, split.by = "batch")


for (i in 1:length(data.list)) {
  data.list[[i]]@active.ident <- factor(data.list[[i]]@meta.data$nnet2)
  names(data.list[[i]]@active.ident) <- colnames(data.list[[i]])
  data.list[[i]]@active.ident <- data.list[[i]]@active.ident[grep("Dendritic|Megakaryocytes",data.list[[i]]@active.ident,invert = T)]
  markers <- FindAllMarkers(data.list[[i]],test.use = "wilcox",only.pos = T)
  file_name <- paste("markers_",names(data.list)[i],".RData",sep="")
  save(markers,file=file_name)
}

files <- list.files(path=".", pattern="*.RData", full.names=F, recursive=FALSE) #### files containing the data.frame with cluster-specific markers computed by FindAllMarkers

load(file=files[1])
df <- data.frame(markers,protocol=rep("C1HT-medium",dim(markers)[1])) ## file 1
for(f in seq(2:14)){
  rm(markers)
  name_data <- unlist(strsplit(files[f],split = ".",fixed = T))[1]
  name_data <- unlist(strsplit(name_data,split = "_",fixed = T))[2]
  load(file=files[f])
  print(dim(markers))
  print(name_data)
  print(f)
  df <- rbind(df,data.frame(markers,protocol=rep(name_data,dim(markers)[1])))
}

source(file="ji_markers.r")
n <- length(files)


J <- lapply(p,function(y) sapply(p, function(x) jaccard_index(df[grep(y,df$protocol),],df[grep(x,df$protocol),],ntop = 100))) ##Jaccard Index
J
dim(J)
91 13
names(J) <- levels(df$protocol)
save(J,file="Jaccard_Index_markers.RData")

y <- sapply(J,function(x) x["B cells",])
df <- data.frame(y)
colnames(df) <- rownames(df)

library(matchSCore2)


png(file="Overlap_markers_Bcells.png",width = 7.5, height = 7, units = 'in', res = 600)
summary_ggplot(y,ylab = "",xlab = "")+scale_fill_gradient(low = "white", high = "blueviolet",name="Jaccard Index \n")+
  theme(legend.title = element_text(size=14),legend.text = element_text(size = 14,colour = "black"))+theme(legend.position = "none")+
  theme(axis.text.x=element_text(colour = "black",angle=80,hjust = 1,size=14),axis.text.y = element_text(colour = "black",size = 16))
dev.off()

y <- sapply(J,function(x) x["HEK cells",])
df <- data.frame(y)
colnames(df) <- rownames(df)



png(file="Overlap_markers_HEKcells.png",width = 8, height = 7, units = 'in', res = 600)
summary_ggplot(y,ylab = "",xlab = "")+scale_fill_gradient(low = "white", high = "blueviolet",name="Jaccard Index \n")+
  theme(legend.title = element_text(size=14),legend.text = element_text(size = 14,colour = "black"))+theme(legend.position = "none")+
  theme(axis.text.x=element_text(colour = "black",angle=80,hjust = 1,size=14),axis.text.y = element_text(colour = "black",size = 16))
dev.off()

png(file="Overlap_markers_HEKcells_legend.png",width = 8, height = 7, units = 'in', res = 600)
summary_ggplot(y,ylab = "",xlab = "")+scale_fill_gradient(low = "white", high = "blueviolet",name="Jaccard Index \n")+
  theme(legend.key.size = unit(1,"cm"),legend.position = "top",legend.title = element_text(size=14),legend.text = element_text(size = 14,colour = "black"))+
  theme(axis.text.x=element_text(colour = "black",angle=80,hjust = 1,size=14),axis.text.y = element_text(colour = "black",size = 16))
dev.off()


y <- sapply(J,function(x) x["CD14+ Monocytes",])
df <- data.frame(y)
colnames(df) <- rownames(df)


png(file="Overlap_markers_CD14_Monocytes.png",width = 8, height = 7, units = 'in', res = 600)
summary_ggplot(y,ylab = "",xlab = "")+scale_fill_gradient(low = "white", high = "blueviolet",name="Jaccard Index \n")+
  theme(legend.title = element_text(size=14),legend.text = element_text(size = 14,colour = "black"))+theme(legend.position = "none")+
  theme(axis.text.x=element_text(colour = "black",angle=80,hjust = 1,size=14),axis.text.y = element_text(colour = "black",size = 16))
dev.off()

#### Marker accuracy - Reference vs Protocol ########

p <- levels(df$protocol)


load(file="10X/reference_mod/markers.ref3.RData") #### Markers from the reference HCA benchmarking dataset (Chromium ref)
head(markers)
markers <- markers[grep("Ambiguous",markers$cluster,invert = T),]
markers <- rbind(df,data.frame(markers,protocol=rep("Reference",dim(markers)[1])))
df.all <- rbind(df,markers)
p <- levels(df.all$protocol)

p <- p[-14]
df <- df.all
J <- sapply(p,function(x) jaccard_index(df[grep("Reference",df$protocol),],df[grep(x,df$protocol),],ntop = 200))
J

J <- J[c("B cells","CD14+ Monocytes","HEK cells"),]
j <- apply(J,2,function(x) mean(x))
marker_accuracy <- j

save(marker_accuracy,file="marker_accuracy.RData")

load(file="Downsampled_analysis/Human/summary/clus.accuracy.RData")
clus_accuracy <- clus_accuracy[names(marker_accuracy)]

df <- data.frame(clus_accuracy,marker_accuracy,protocol=names(marker_accuracy))



load(file="col.batch.RData")

png(file="Clustering_accuracy_vs_marker_sensitivity_by_protocols.png",width = 8, height = 6,units = 'in', res = 600)
ggplot(df,aes(x=clus_accuracy,y=marker_accuracy))+geom_point(aes(colour=protocol),size=5)+ 
  scale_color_manual(values = col.batch)+
  xlab("Clustering Accuracy")+ylab("Marker Sensitivity")+ theme_bw()+
  theme(panel.grid = element_blank(),axis.line = element_line(size=1),legend.text = element_text(size = 14),axis.text =element_text(size = 24),text = element_text(size = 20)) 
dev.off()

