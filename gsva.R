colnames(Epithelial@meta.data)

#载入R包
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(SingleR)
library(GSVA)
library(GSEA)
library(GSEABase)
library(clusterProfiler)
library(limma)
library(pheatmap)
library(xtable)
library(reshape2)
library(ComplexHeatmap)

expr <- as.data.frame(Epithelial@assays$RNA@data)
kegggmt <- read.gmt("/home/data/chendemeng/gmtfiles/mm.c2.cp.kegg.v7.2.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
expr=as.matrix(expr)
kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)

write.table(kegg2, '/home/data/maomao/bca_anticd276_matrix/harmony/Epithelial cells/gsva figures.xls', row.names=T, col.names=NA, sep='\t')

colnames(basal.sub@meta.data)
meta <- basal.sub@meta.data[,c("RNA_snn_res.0.2","orig.ident")]

##每个细胞类别与功能相关热图
meta <- meta %>%arrange(meta$RNA_snn_res.0.2)
data <- kegg2[,rownames(meta)]
group <- factor(meta[,"RNA_snn_res.0.2"],ordered = F)
data1 <-NULL
for(i in 0:(length(unique(group))-1)){
  ind <-which(group==i)
  dat <- apply(data[,ind], 1, mean)
  data1 <-cbind(data1,dat)
}
colnames(data1) <-c("C1","C2","C3","C4","C5","C6","C7")
result<- t(scale(t(data1)))
result[result>2]=2
result[result<-2]=-2
library(pheatmap)
p <- pheatmap(result[1:30,],
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("blue", "white","red"))(100),
              cellwidth = 10, cellheight = 15,
              fontsize = 10)
ggsave("Basalgsva_celltype.png",plot=p, width = 7,height = 7)


