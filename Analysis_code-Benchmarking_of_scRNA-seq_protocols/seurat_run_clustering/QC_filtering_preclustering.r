load(file="ProtocolName.hsap/mmus.full.SCE.Robj")
sce <- full.SCE.hsap 
pd=as.data.frame(colData(sce))
names(pd)

mapped <- rowSums(pd[,c("nExonReads","nIntronReads","nIntergenicReads","nAmbiguityReads")])
Nreads <- pd$nTReads
PCTmapped <- mapped/Nreads*100              

pd<-data.frame(Mapped=mapped,PCTmapped,pd)


sce@colData<-DataFrame(pd)
sce
sce <- subset(sce,i=1:dim(sce)[1],j=which(sce@colData$Species=="Human"))


x=log10(pd$Mapped)
summary(x)
plot(density(x))


by_mapping=colnames(sce)[which((x)< 4 | pd$PCTmapped<65)] ### low-mapping cells



l=length(colData(sce))
l
low_map_cell <- ifelse(colnames(sce) %in% by_mapping,TRUE,FALSE)
colData(sce)[l+1] <- low_map_cell 
names(colData(sce))[l+1]<- "is_low_quality"

df=as.data.frame(colData(sce))
is_low_quality <- df$is_low_quality


sce
hsap_filt <- subset(sce,i=1:dim(sce)[1],j=which(is_low_quality==FALSE))


counts <- counts(hsap_filt)
head(counts[,1:4])

cells_names=colnames(counts)

genes=read.table("~/Dropbox/annotation/human/gencode_28.tsv",header = T)

cd1=data.frame(id=rownames(counts),counts)
genes_cd=merge(genes,cd1,by="id",all.x=F)


rownames(genes_cd)=make.unique(as.character(genes_cd$name))
counts=genes_cd[,6:dim(genes_cd)[2]]
save(counts,file="counts_hsap_filtered.RData")

hsap_filt@assays$data$counts <- counts
rownames(hsap_filt) <- rownames(counts)
save(hsap_filt,file="hsap_filt.RData")


col=c("cornflowerblue","aquamarine3","chocolate1","lightslategray","darkorchid","firebrick1","magenta","brown1","turquoise4","hotpink","royalblue","green4")
col1=c("firebrick1","blue","turquoise4","green4","royalblue","violet")

library(ggplot2)

text <- element_text(size = 20)

png(file="Mapped_vs_pMapped.png",width = 6, height = 6, units = 'in', res = 600)
ggplot(pd,aes(x=mapped,y=PCTmapped))+geom_point(aes(colour = Species),size=2,alpha=0.8)+scale_color_manual(values = col)+theme_bw()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_log10()+
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()


intronic <- pd$nIntronReads/mapped*100
exonic <- pd$nExonReads/mapped*100

pd<-data.frame(exonic=exonic,intronic=intronic,pd)

png(file="Mapped_vs_Intronic.pct.png",width = 6, height = 6, units = 'in', res = 600)
ggplot(pd,aes(x=mapped,y=intronic))+geom_point(aes(colour=species),size=2,alpha=0.8)+ theme_bw()+
  scale_x_log10()+scale_color_manual(values = col)+
  labs(x = "Mapped",y="Intronic content (%)")+ 
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()

png(file="Mapped_vs_Exonic.pct.png",width = 6, height = 6, units = 'in', res = 300)
ggplot(pd,aes(x=mapped,y=exonic))+geom_point(aes(colour=species),size=2,alpha=0.8)+ theme_bw()+
  scale_x_log10()+scale_color_manual(values = col)+
  labs(x = "Mapped",y="Exonic content (%)")+ 
  theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()

png(file="Mapped_vs_pMapped_by_HighQualitycells.png",width = 6, height = 6, units = 'in', res = 300)
ggplot(df,aes(x=df$Mapped,y=df$PCTmapped))+geom_point(aes(colour = is_low_quality),size=2,alpha=0.8)+ theme_bw()+
  scale_x_log10()+
  labs(x = "Mapped",y="Mapped (%)")+scale_color_manual(values=c("gray","blue"),labels=c( "FALSE","TRUE"))+
  theme(axis.line.x = element_line(color="black", size = 1,l),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()