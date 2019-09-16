setwd("~/Dropbox/HCA_benchmarking_sample/10X/")
load(file="10X8x10K.hsap.full.SCE.Robj")
sce <-full.SCE.hsap 
rm(full.SCE.hsap)
pd=as.data.frame(colData(sce))
names(pd)
[1] "nTReads"                        "nExonReads"                     "nIntronReads"                   "nIntergenicReads"              
[5] "nUnmappedReads"                 "nAmbiguityReads"                "nUMIs"                          "nGenes"                        
[9] "Species"                        "Library"                        "is_cell_control"                "total_features_by_counts"      
[13] "log10_total_features_by_counts" "total_counts"                   "log10_total_counts"             "pct_counts_in_top_50_features" 
[17] "pct_counts_in_top_100_features" "pct_counts_in_top_200_features" "pct_counts_in_top_500_features"
table(pd$Species)
# Doublet   Human   Mouse     Dog 
# 1067    3520   43148   32265

mapped <- rowSums(pd[,c(2:4,6)],na.rm = T)
Nreads <- rowSums(pd[,c(2:6)],na.rm = T)
PCTmapped <- mapped/Nreads*100              

pd<-data.frame(Mapped=mapped,Nreads,PCTmapped,pd)
head(pd)

summary(pd$Mapped[which(pd$Species=="Human")])

col=c("cornflowerblue","aquamarine3","chocolate1","lightslategray","darkorchid","firebrick1","magenta","brown1","turquoise4","hotpink","royalblue","green4")
col1=c("firebrick1","blue","turquoise4","green4","royalblue","violet")

library(ggplot2)

text <- element_text(size = 20)
species <- pd$Species

png(file="Mapped_vs_pMapped.png",width = 6, height = 6, units = 'in', res = 300)
ggplot(pd,aes(x=mapped,y=PCTmapped))+geom_point(aes(colour = species),size=2,alpha=0.8)+scale_color_manual(values = col)+theme_bw()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_log10()+
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()


intronic <- pd$nIntronReads/mapped*100
exonic <- pd$nExonReads/mapped*100
pd<-data.frame(exonic=exonic,intronic=intronic,pd)

png(file="Mapped_vs_Intronic.pct.png",width = 6, height = 6, units = 'in', res = 300)
ggplot(pd,aes(x=mapped,y=intronic))+geom_point(aes(colour=species),size=2,alpha=0.8)+ theme_bw()+
  scale_x_log10()+scale_color_manual(values = col)+
  #scale_color_gradientn(colours = rainbow(6))  +
  labs(x = "Mapped",y="Intronic content (%)")+ 
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()

png(file="Mapped_vs_Exonic.pct.png",width = 6, height = 6, units = 'in', res = 300)
ggplot(pd,aes(x=mapped,y=exonic))+geom_point(aes(colour=species),size=2,alpha=0.8)+ theme_bw()+
  scale_x_log10()+scale_color_manual(values = col)+
  #scale_color_gradientn(colours = rainbow(6))  +
  labs(x = "Mapped",y="Exonic content (%)")+ 
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()

sce@colData<-DataFrame(pd)

sce <- subset(sce,i=1:dim(sce)[1],j=which(sce@colData$Species=="Human"))

species <- pd$Species


pd=as.data.frame(sce@colData)
x=log10(pd$Mapped)
summary(x)
plot(density(x))


plot(density(x),main="Density Min library size")

summary(x)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.301   3.007   3.262   3.450   3.903   5.403 
med=median(x)
mad=mad(x,center = med, constant=3)
cut.off1=med-mad
cut.off1
# 3.24

by_mapping=colnames(sce)[which((x)< cut.off1 | pd$PCTmapped<65)] ### low mapping cells


pd=colData(sce)

l=length(colData(sce))
l
low_map_cell=ifelse(colnames(sce) %in% by_mapping,TRUE,FALSE)
colData(sce)[l+1]=low_map_cell 
names(colData(sce))[l+1]="is_low_quality"

df=as.data.frame(colData(sce))
is_low_quality=df$is_low_quality

table(is_low_quality)
is_low_quality
FALSE  TRUE 
42925   223 
n=length(is_low_quality)


sce
hsap_filt <- subset(sce,i=1:dim(sce)[1],j=which(is_low_quality==FALSE))


counts <- counts(hsap_filt)
cells_names=colnames(counts)

genes=read.table("~/Dropbox/annotation/human/gencode_28.tsv",header = T)


cd1=data.frame(id=rownames(counts),counts)
genes_cd=merge(genes,cd1,by="id",all.x=F)
rownames(genes_cd)=make.unique(as.character(genes_cd$name))
counts=genes_cd[,6:dim(genes_cd)[2]]

hsap_filt@assays$data$counts <- counts
rownames(hsap_filt) <- rownames(counts)
save(hsap_filt,file="hsap_filt.RData")


