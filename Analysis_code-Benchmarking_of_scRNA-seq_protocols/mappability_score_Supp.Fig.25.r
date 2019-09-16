load(file="Merged/sce.all_classified.technologies.RData")

cells <- colnames(sce)
save(cells,file="cells.RData")
tech <- factor(sce@colData$batch)
table(tech)
# tech
# C1HT-medium    C1HT-small      CEL-Seq2      Chromium Chromium (sn)         ddSEQ      Drop-Seq        ICELL8 
# 2216          1606          1074          1604          1399          2109          2261          1926 
# inDrop      MARS-Seq   Quartz-Seq2      SCRB-Seq    Smart-Seq2 
# 686          1480          1326          1397           721 

metadata <- data.frame(sce@colData)

library(Seurat)

x <- sapply(levels(tech),function(x) which(tech==x))

counts <- counts(sce)

raw <- sapply(x,function(y) counts[,y])
meta <- lapply(x,function(y) metadata[y,])

data <- sapply(raw,function(y) CreateSeuratObject(raw.data = y, min.cells = 10, min.genes = 0, project = "DS20K"))
data <- sapply(data,function(y) NormalizeData(object = y, normalization.method = "LogNormalize",scale.factor = 10000))
data <- sapply(data,function(y) ScaleData(object = y,vars.to.regress = "nUMI")) 

load(file = "10X/reference_mod/mod_referenceHCA_10X_40K_acc0.96_reference3.RData") #### matchSCore2 model trained by using the HCA benchmarking reference (Chromium) 
load(file = "10X/reference_mod/gene_cl.ref.RData") #### gene signature from the HCA reference

library(Matrix)
library(nnet)

var.test2 <- lapply(data, function(y) sapply(gene_cl.ref, function(z) colSums(y@scale.data[which(rownames(y@scale.data) %in% z),])))
var.test2 <- lapply(var.test2, function(y) t(apply(y,1, function(z) (z-min(z))/(max(z)-min(z)))))
fitted.results <- lapply(var.test2, function(y) predict(mod, newdata = data.frame(y), "probs"))
fit <- sapply(fitted.results,function(y) apply(y,1,function(z) colnames(y)[which(z==max(z)[1])]))

nnet <- sce@colData$nnet2
res <- lapply(levels(tech), function(x) data.frame(nnet=nnet[which(tech==x)],pred=fit[[x]],prob=fitted.results[[x]]))
names(res) <- levels(tech)
head(res)

bcells <- lapply(res, function(x) data.frame(row.names = rownames(x)[which(x$nnet=="B cells")],prob=x$prob.B.cells[which(x$nnet=="B cells")]))
hek <- lapply(res, function(x) data.frame(row.names = rownames(x)[which(x$nnet=="HEK cells")],prob=x$prob.HEK.cells[which(x$nnet=="HEK cells")]))
monocytes <- lapply(res, function(x) data.frame(row.names = rownames(x)[which(x$nnet=="CD14+ Monocytes")],prob=x$prob.CD14..Monocytes[which(x$nnet=="CD14+ Monocytes")]))

b.prob <- unlist(bcells)
tech.b <- sapply(names(b.prob),function(y) unlist(strsplit(y,split = ".",fixed = T))[1]) 
b.df <- data.frame(prob=b.prob,cell_type=rep("B",length(b.prob)),method=tech.b)
ave.b <- sapply(levels(tech),function(x) mean(b.df$prob[which(b.df$method==x)]))

hek.prob <- unlist(hek)
tech.hek <- sapply(names(hek.prob),function(y) unlist(strsplit(y,split = ".",fixed = T))[1]) 
hek.df <- data.frame(prob=hek.prob,cell_type=rep("HEK cells",length(hek.prob)),method=tech.hek)
ave.hek <- sapply(levels(tech),function(x) mean(hek.df$prob[which(hek.df$method==x)]))

mono.prob <- unlist(monocytes)
tech.mono <- sapply(names(mono.prob),function(y) unlist(strsplit(y,split = ".",fixed = T))[1]) 
mono.df <- data.frame(prob=mono.prob,cell_type=rep("CD14+ Monocytes",length(mono.prob)),method=tech.mono)
ave.mono <- sapply(levels(tech),function(x) mean(mono.df$prob[which(mono.df$method==x)]))

df20 <- list(b=b.df,hek=hek.df,mono=mono.df)
save(df20,file="df20_allsamples.RData")




ave.mapping <- data.frame(cbind(ave.b,ave.hek,ave.mono))
names(ave.mapping) <- c("B cells","HEK cells","Monocytes")
ave.mapping <- data.frame(mapping=ave.mapping,method=rownames(ave.mapping))
names(ave.mapping)[1:3] <- c("B cells","HEK cells","Monocytes")


library(reshape2)
load(file="col.batch.RData")
mappability<- melt(ave.mapping)
mappability <- data.frame(mappability,depth=rep("20K",dim(mappability)[1]))

names(mappability)[2:3] <- c("cell_type","prob")
m <- factor(mappability$method) 
mappability$method <- m

dir.create("mappability")
png(file="mappability/mappability.png",width = 8, height = 6,units = 'in', res = 600)
ggplot(mappability,aes(x=cell_type,y=prob,fill=method))+geom_point(aes(colour=method),size=4)+ 
  geom_jitter(aes(colour=method))+scale_colour_manual(values = col.batch)+theme(legend.text = element_text(size = 14),
  axis.text =element_text(size = 16),text = element_text(size = 20)) + 
  xlab("")+ylab("Classification probability")
dev.off()  

map20k <- mappability
save(map20k,file="mappability/ave_map20k.RData") 

