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
#创建文件夹
dir.create("/home/data/maomao/bca_anticd276_matrix/harmony")
setwd("/home/data/maomao/bca_anticd276_matrix/harmony")
getwd()

DMSO1 <- read.table("Bca-WT1_matrix.tsv.gz")
DMSO2 <- read.table("Bca-WT2_matrix.tsv.gz")
antiCD276_1 <- read.table("Bca-CD276-1_matrix.tsv.gz")
antiCD276_2 <- read.table("Bca-cd276I_matrix.tsv.gz")

DMSO1 <- CreateSeuratObject(counts = DMSO1, project = "DMSO1", min.cells = 3, min.features = 200)
DMSO2 <- CreateSeuratObject(counts = DMSO2, project = "DMSO2", min.cells = 3, min.features = 200)
antiCD276_1 <- CreateSeuratObject(counts = antiCD276_1, project = "antiCD276_1", min.cells = 3, min.features = 200)
antiCD276_2 <- CreateSeuratObject(counts = antiCD276_2, project = "antiCD276_2", min.cells = 3, min.features = 200)

##########计算质控指标#计算细胞中线粒体基因比例
DMSO1[["percent.mt"]] <- PercentageFeatureSet(DMSO1, pattern = "^mt-")
DMSO2[["percent.mt"]] <- PercentageFeatureSet(DMSO2, pattern = "^mt-")
antiCD276_1[["percent.mt"]] <- PercentageFeatureSet(antiCD276_1, pattern = "^mt-")
antiCD276_2[["percent.mt"]] <- PercentageFeatureSet(antiCD276_2, pattern = "^mt-")

#计算红细胞比例
HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(DMSO1@assays$RNA)) 
HB.genes <- rownames(DMSO1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
DMSO1[["percent.HB"]]<-PercentageFeatureSet(DMSO1, features=HB.genes) 

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(DMSO2@assays$RNA)) 
HB.genes <- rownames(DMSO2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
DMSO2[["percent.HB"]]<-PercentageFeatureSet(DMSO2, features=HB.genes) 

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(antiCD276_1@assays$RNA)) 
HB.genes <- rownames(antiCD276_1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
antiCD276_1[["percent.HB"]]<-PercentageFeatureSet(antiCD276_1, features=HB.genes) 

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(antiCD276_2@assays$RNA)) 
HB.genes <- rownames(antiCD276_2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
antiCD276_2[["percent.HB"]]<-PercentageFeatureSet(antiCD276_2, features=HB.genes) 

##数据质控###每个sample单独run一次
DMSO1 <- subset(DMSO1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
DMSO2 <- subset(DMSO2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
antiCD276_1 <- subset(antiCD276_1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
antiCD276_2 <- subset(antiCD276_2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)


###########手动去doubets
DMSO1 <- NormalizeData(DMSO1)
DMSO1 <- FindVariableFeatures(DMSO1, selection.method = "vst", nfeatures = 2000)
DMSO1 <- ScaleData(DMSO1, features = VariableFeatures(object = DMSO1))
DMSO1 <- RunPCA(DMSO1, verbose = FALSE)
DMSO1 <- RunUMAP(DMSO1, dims = 1:30)
sweep.res.list <- paramSweep_v3(DMSO1, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- DMSO1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(DMSO1@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
DMSO1 <- doubletFinder_v3(DMSO1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DMSO1 <- doubletFinder_v3(DMSO1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(DMSO1@meta.data)[6], sct = FALSE)
colnames(DMSO1@meta.data)
DMSO1 <- subset(DMSO1, subset = DF.classifications_0.25_0.04_983 == "Singlet")


DMSO2 <- NormalizeData(DMSO2)
DMSO2 <- FindVariableFeatures(DMSO2, selection.method = "vst", nfeatures = 2000)
DMSO2 <- ScaleData(DMSO2, features = VariableFeatures(object = DMSO2))
DMSO2 <- RunPCA(DMSO2, verbose = FALSE)
DMSO2 <- RunUMAP(DMSO2, dims = 1:30)
sweep.res.list <- paramSweep_v3(DMSO2, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- DMSO2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(DMSO2@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
DMSO2 <- doubletFinder_v3(DMSO2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DMSO2 <- doubletFinder_v3(DMSO2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(DMSO2@meta.data)[6], sct = FALSE)
colnames(DMSO2@meta.data)
DMSO2 <- subset(DMSO2, subset = DF.classifications_0.25_0.22_525 == "Singlet")


antiCD276_1 <- NormalizeData(antiCD276_1)
antiCD276_1 <- FindVariableFeatures(antiCD276_1, selection.method = "vst", nfeatures = 2000)
antiCD276_1 <- ScaleData(antiCD276_1, features = VariableFeatures(object = antiCD276_1))
antiCD276_1 <- RunPCA(antiCD276_1, verbose = FALSE)
antiCD276_1 <- RunUMAP(antiCD276_1, dims = 1:30)
sweep.res.list <- paramSweep_v3(antiCD276_1, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- antiCD276_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(antiCD276_1@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
antiCD276_1 <- doubletFinder_v3(antiCD276_1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
antiCD276_1 <- doubletFinder_v3(antiCD276_1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(antiCD276_1@meta.data)[6], sct = FALSE)
colnames(antiCD276_1@meta.data)
antiCD276_1 <- subset(antiCD276_1, subset = DF.classifications_0.25_0.28_461 == "Singlet")


antiCD276_2 <- NormalizeData(antiCD276_2)
antiCD276_2 <- FindVariableFeatures(antiCD276_2, selection.method = "vst", nfeatures = 2000)
antiCD276_2 <- ScaleData(antiCD276_2, features = VariableFeatures(object = antiCD276_2))
antiCD276_2 <- RunPCA(antiCD276_2, verbose = FALSE)
antiCD276_2 <- RunUMAP(antiCD276_2, dims = 1:30)
sweep.res.list <- paramSweep_v3(antiCD276_2, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- antiCD276_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(antiCD276_2@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
antiCD276_2 <- doubletFinder_v3(antiCD276_2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
antiCD276_2 <- doubletFinder_v3(antiCD276_2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(antiCD276_2@meta.data)[6], sct = FALSE)
colnames(antiCD276_2@meta.data)
antiCD276_2 <- subset(antiCD276_2, subset = DF.classifications_0.25_0.02_407 == "Singlet")


##merge+harmony去批次
bca.integrated <- merge(x = DMSO1, y = list (DMSO2,antiCD276_1,antiCD276_2))
bca.integrated
head(bca.integrated@meta.data)
library("harmony")
bca.integrated <- NormalizeData(bca.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
bca.integrated <- FindVariableFeatures(bca.integrated, selection.method = "vst", nfeatures = 2000)
bca.integrated<-ScaleData(bca.integrated)%>%RunPCA(verbose=FALSE)
system.time({bca.integrated <- RunHarmony(bca.integrated, group.by.vars = "orig.ident")})

bca.integrated<-FindNeighbors(bca.integrated, reduction = "harmony",dims = 1:30)
bca.integrated<-FindClusters(bca.integrated, reduction = "harmony",resolution = 0.1)
bca.integrated<-RunUMAP(bca.integrated,dims = 1:30,reduction = "harmony")
bca.integrated<-RunTSNE(bca.integrated,dims = 1:30,reduction = "harmony")

DimPlot(bca.integrated, reduction = "umap",label = T, pt.size = 0.5)
bca.integrated$new_group <- ifelse(bca.integrated$orig.ident %in% c("DMSO1","DMSO2"),"DMSO","antiCD276") 


############去线粒体核糖体基因
mt_rp0<-read.table("/home/data/wangxiaochen/mm_excludeGene.txt",header = T,sep = "\t",stringsAsFactors = F)
length(unique(mt_rp0[,1]))==nrow(mt_rp0)
length(rownames(bca.integrated))
grep("^mt-",rownames(bca.integrated))
length(grep("^mt-",rownames(bca.integrated)))
rownames(bca.integrated)[grep("^mt-",rownames(bca.integrated))]
grep("^Rp[ls]",rownames(bca.integrated))
length(grep("^Rp[ls]",rownames(bca.integrated)))
rownames(bca.integrated)[grep("^Rp[ls]",rownames(bca.integrated))]
mt_rp<-unique(c(mt_rp0[,1],(rownames(bca.integrated)[grep("^mt-",rownames(bca.integrated))]),(rownames(bca.integrated)[grep("^Rp[ls]",rownames(bca.integrated))])))
mt_rp
keepgene<-setdiff(rownames(bca.integrated),mt_rp)
bca.integrated<-subset(bca.integrated,features=keepgene)
bca.integrated


######Normalization
bca.integrated <- NormalizeData(bca.integrated)
bca.integrated <- FindVariableFeatures(bca.integrated, selection.method = "vst", nfeatures = 2000)
bca.integrated <- ScaleData(bca.integrated, features = VariableFeatures(object = bca.integrated))
bca.integrated <- RunPCA(bca.integrated, features = VariableFeatures(object = bca.integrated))
ElbowPlot(bca.integrated,ndims = 50)
bca.integrated <- FindNeighbors(bca.integrated, dims = 1:20,k.param = 20)
bca.integrated <- FindClusters(bca.integrated, resolution = 0.2)
bca.integrated <- RunUMAP(bca.integrated, dims = 1:20)



#####图像和另外一种去doublet的方法，每去一次细胞都要重新normalized
DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)

DimPlot(seurat0,reduction = "umap",label = TRUE,pt.size = 0.5)
bca.integrated
seurat0
seurat0<-bca.integrated
co_exp1<-subset(seurat0,subset=Epcam > 0 & Ptprc > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp1)]
co_exp2<-subset(seurat0,subset=Epcam > 0 & Cd3d > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp2)]
co_exp3<-subset(seurat0,subset=Cd79a > 0 & Cd3d > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp3)]
co_exp4<-subset(seurat0,subset=Col1a1 > 0 & Ptprc > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp4)]
co_exp5<-subset(seurat0,subset=Cdh5 > 0 & Ptprc > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp5)]
co_exp6<-subset(seurat0,subset=Cdh5 > 0 & Epcam > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp6)]
co_exp7<-subset(seurat0,subset=Cdh5 > 0 & Col1a1 > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp7)]
co_exp8<-subset(seurat0,subset=Flt1 > 0 & Myl9 > 0)
seurat0<-seurat0[,!colnames(seurat0) %in% colnames(co_exp8)]
bca.integrated <- seurat0



####singleR####IGD for immune cells#########
load("/home/data/SingleR_ref/ref_Mouse_imm.RData")
singler_ref <- ref_IGD
for (i in 0:13) {
  subset_data<-subset(bca.integrated,idents = i)
  normlized_data <- as.matrix(subset_data@assays$RNA@data)
  pred.cluster <- SingleR(test = normlized_data, ref = singler_ref, labels=singler_ref$label.main)
  assign(paste("cluster_",i,sep = ""), table(pred.cluster$first.labels))
}

for (i in 0:13) {
  cluster_table <- as.data.frame(get(paste("cluster_",i,sep = "")))
  names(cluster_table) <- c("Cell_Type",paste0("cluster",i))
  if (i == 0) final_table <- cluster_table else final_table <- merge(final_table, cluster_table, by = "Cell_Type", all.x =T, all.y = T)
}

final_table[is.na(final_table)] <- 0

write.csv(final_table, file = "identify_clusters_immune.csv", row.names = F, quote = F)


############根据umap的分布来去双细胞，每去一次细胞都要重新normalized
table(Idents(bca.integrated), bca.integrated$seurat_clusters)
bca.integrated$keep <- ifelse(bca.integrated$seurat_clusters %in% c("13","14"),"not","keep")
bca.integrated <- subset(bca.integrated, subset = keep == "keep")
dim(bca.integrated)

################分群#########
bca.integrated@meta.data$new.lables<-"Epithelial cells"
#bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(14))]<-"Stromal cells"#
bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(5,8,10,11))]<-"Fibroblast"#
bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(9))]<-"Endothelial cells"#
#bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(5,8))]<-"EMT cells"#
bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(13))]<-"Lymphoid cells"#
bca.integrated@meta.data$new.lables[which(bca.integrated@meta.data$seurat_clusters %in% c(2,12))]<-"Myeloid cells"
Idents(bca.integrated) <- "new.lables"
colnames(bca.integrated@meta.data)
#########作图
p1 <- DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5)
getwd()
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_umap.pdf", plot = p1, width = 8, height = 6, dpi = 1000) 

DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5, group.by = "orig.ident")

DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5, group.by = "new_group")

DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5,split.by = "seurat_clusters")

DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)

colnames(bca.integrated@meta.data)
# Specify genes  气泡图
genes_to_check = c("Krt19","Krt5","Krt14","Krt8","Dcn","Col1a2","Col3a1","Col1a1","Cdh5","Cd34","Aqp1","Lrg1",
                   "C1qc","Cd68","Cd74","Cd3d","Rac2","Trac","Cd3g")
# All on Dotplot 
p <- DotPlot(bca.integrated, features = genes_to_check,  cols = c("lightgrey", "red")) + coord_flip()
p
ggsave2("气泡图.tiff", plot = p, width = 12, height = 8, dpi = 300)
ggsave2("气泡图.pdf", plot = p, width = 12, height = 8, dpi = 300)

getwd()
DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)
p1<-DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5)
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_umap.pdf", plot = p1, width = 8, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_umap.png", plot = p1, width = 8, height = 6, dpi = 1000) 

p2<-DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5, group.by = "orig.ident")
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_orig.ident.pdf", plot = p2, width = 8, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_orig.ident.png", plot = p2, width = 8, height = 6, dpi = 1000) 

p3<-DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5, group.by = "new_group")
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_new_group.pdf", plot = p3, width = 8, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_new_group.png", plot = p3, width = 8, height = 6, dpi = 1000) 

p4<-DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5, group.by = "RNA_snn_res.0.2")
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_clusters.pdf", plot = p4, width = 8, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_clusters.png", plot = p4, width = 8, height = 6, dpi = 1000) 

p5<-DimPlot(bca.integrated,reduction = "umap",label = F,pt.size = 0.5, split.by = "new_group")
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_split_new_group.pdf", plot = p5, width = 12, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_split_new_group.png", plot = p5, width = 12, height = 6, dpi = 1000) 


DimPlot(bca.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)

############################cell cycle analysis#############
colnames(bca.integrated@meta.data)
DefaultAssay(bca.integrated) <- "RNA"
bca.integrated<- CellCycleScoring(object = bca.integrated, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = bca.integrated@meta.data)
p6<-DimPlot(bca.integrated,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_phase.pdf", plot = p6, width = 8, height = 6, dpi = 1000) 
ggsave("/home/data/maomao/bca_anticd276_matrix/harmony/all_phase.png", plot = p6, width = 8, height = 6, dpi = 1000) 



#########find markers for each clusters算每个群的特异表达基因##########
bca.integrated.markers <- FindAllMarkers(bca.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
bca.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(bca.integrated.markers,"bca.integrated.markers.csv",col.names=T,row.names=T,sep=",")

top5 <- bca.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DefaultAssay(bca.integrated) <- "RNA"
bca.integrated@assays$RNA@scale.data <- scale(bca.integrated@assays$RNA@data, scale = TRUE)
DoHeatmap(bca.integrated, features = top5$gene) + NoLegend()



top2 <- bca.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DoHeatmap(bca.integrated, features = top10$gene)+NoLegend()
FeaturePlot(bca.integrated, features = c("Nat10","Pcif1","Sox2","Bmi1","Axin2","Lamc2"),reduction = "umap")
FeaturePlot(bca.integrated, features = c("Pcif1"),reduction = "umap")
FeaturePlot(bca.integrated, features = c("Mettl3"),reduction = "umap")
FeaturePlot(bca.integrated, features = c("Myl9"),reduction = "umap")
FeaturePlot(bca.integrated, features = c("Flt1"),reduction = "umap")

# Specify genes  气泡图
genes_to_check = c("Krt14","Krt17","Krt6a","Il17a","Rgs1","Cd3g","Dcn","Gsn","Col1a1","Il1b","Cxcl2","Cd74","Mylk","Myl9","Myh11","Flt1","Pecam1","Eng")
p7<-DotPlot(bca.integrated, features = genes_to_check) 
ggsave("/home/data/wangxiaochen/scRNAseq/06072021/all_markers_dotblot.pdf", plot = p7, width = 12, height = 6, dpi = 1000) 
ggsave("/home/data/wangxiaochen/scRNAseq/06072021/all_markers_dotblot.png", plot = p7, width = 12, height = 6, dpi = 1000) 

##保存数据
saveRDS(bca.integrated,'bca.integrated.rds')
bca.integrated <- readRDS('bca.integrated.rds')
