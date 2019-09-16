load(file="data_seu.obj.RData") ### data: seurat object after integration across protocols. Downsampled data to 20K reads. 

load(file="col.batch.RData")
load(file="col_PBMC_human.RData")


tech <- levels(factor(data@meta.data$batch)) ## tech: vector of protocols
nnet <- factor(data@meta.data$nnet2) ## nnet: vector of cell identities assigned by matchSCore2
cell_types <- levels(nnet)[c(1:4,7,9)] ## cell_type: vector of cell types to be considered (B cell, CD4+ T-cells, CD8+ T-cells, NK cells, HEK293T cells)
cell_types



load(file="acc_by_ct.RData") ### clustering accuracy
load(file="mixing_scores.RData") ### mixability scores



scores <- scores[names(acc_by_ct)]
df <- data.frame(acc=acc_by_ct,mix=scores)
tech <- rownames(df) ## protocols


library(ggplot2)
png(file="cell_type_sep_vs_mixability.png",width = 8, height = 5,units = 'in', res = 600)
ggplot(df,aes(x=acc,y=mix))+geom_point(aes(colour=tech),size=5)+labs(colour="Method")+theme_bw()+
  xlim(c(0.5,0.8))+ylim(c(0.6,0.9))+
  xlab("Cell Type Separation")+ylab("Mixability")+ scale_color_manual(values = col.batch)+
  theme(panel.grid = element_blank(),axis.line = element_line(size=1),legend.text = element_text(size = 14),axis.text =element_text(size = 24),text = element_text(size = 20)) 
dev.off()
