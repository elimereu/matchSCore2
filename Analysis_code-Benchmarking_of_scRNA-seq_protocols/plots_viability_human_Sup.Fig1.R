load(file="10XV3Rev.hsap.full.SCE.20Kcells.Robj")

library(scater)

sce <-filter(full.SCE.hsap,Species=="Human")
load(file="mito.genes_ENS.ids.RData")
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = which(rownames(sce) %in% mito.genes)))


pd <- data.frame(colData(sce))

mapped <- rowSums(pd[,c(2:4)]) ### exons,introns,intergenic reads counts
Nreads <- pd[,1] ### number of total reads
pd$PCTmapped <- mapped/Nreads
pd$Mapped <- mapped

Library <- pd$Library
table(Library)
Library
# AN4492 AN4491 AN4493 
# 6061   7072   9754 


load(file="ChromV3_EmptyDroplets_selected_cells_human_and_mouse.RData")


pd1 <- pd[which(rownames(pd) %in% selected.cells$Lib91),] 
dim(pd1)
# 7031   48
pd2 <- pd[which(rownames(pd) %in% selected.cells$Lib92),] 
dim(pd2)
#6047   48


library(ggplot2)

text <- element_text(size = 24)


pd_new <- rbind(pd1,pd2)
pd_new$Library <- factor(pd_new$Library)
levels(pd_new$Library) <- rev(c("Viability","Without.Viability"))
table(pd_new$Library)
# Without.Viability         Viability 
# 6047              7031 


high_mito <- rownames(pd_new)[which(pd_new$pct_counts_Mt>25)]
is_high_mito <- ifelse(rownames(pd_new) %in% high_mito,TRUE,FALSE)
lib <- pd_new$Library

x <- log10(pd_new$Mapped)
plot(density(x))
abline(v=3.25,col="red",lty=2)
by_mapping=rownames(pd_new)[which((x)< 3.25 & pd_new$PCTmapped< 65)] ### low mapping cells
low_map_cell=ifelse(rownames(pd_new) %in% by_mapping,TRUE,FALSE)
low_quality <- ifelse(rownames(pd_new) %in% union(by_mapping,high_mito),TRUE,FALSE)


png(file="nGenesVS_MTcontent_shape_by_lowquality.png",width = 8, height = 6, units = 'in', res = 600)
ggplot(pd_new,aes(x=nGenes,y=pct_counts_Mt))+geom_jitter(aes(colour=Library,shape=low_quality),alpha=0.6,size=1.5)+theme_bw()+scale_color_manual(values = c("red","blue","green4"))+
  labs(x ="nGenes (log-scale)",y="% MT content")+scale_x_log10()+
  theme(panel.grid = element_blank(),axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 1), 
        axis.text= text,axis.title = text,legend.text = text,legend.title = text)
dev.off()


png(file="Mapped_vs_Ngenes_20K_shape_by_highMito.png",width = 8, height = 6, units = 'in', res = 600)
ggplot(pd_new,aes(x=Mapped,y=nGenes))+geom_point(aes(colour = Library,shape=is_high_mito),size=1.5,alpha=0.5)+geom_vline(xintercept = c(1778),linetype="dotted")+theme_bw()+scale_color_manual(values = rev(c("blue","red")))+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+scale_y_log10()+labs(y = "nGenes (log-scale)",x="Mapped (log-scale)")+scale_x_log10()+
  theme(panel.grid = element_blank(),axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.6), axis.text= text,axis.title = text)
dev.off()


