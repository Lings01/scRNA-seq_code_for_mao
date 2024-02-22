devtools::install_github("sqjin/CellChat")
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
# ⚠️要根据自己python3的路径来设置，可以在终端使用which python3来查看 
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)


cellchat <- createCellChat(object=bca.integrated_1,group.by = "new.lables")
#cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "cell_type")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.mouse
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object


## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.mouse) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")



#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
# save as TIL/net_number_strength.pdf
dev.off()
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_number_individual.pdf

## 运行上述代码出现报错 Error in plot.new() : figure margins too large
# par("mar")
## [1] 5.1 4.1 4.1 2.1
# par(mar=c(1,1,1,1))
# 重新运行上面的代码

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_strength_individual.pdf

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_strength_individual.pdf
table(bca.integrated_1$new.lables)

cellchat@netP$pathways  #查看都有哪些信号通路
pathways.show <- c("TGFb") 

levels(cellchat@idents)    # show all celltype
# [1] "Naive CD4 T"  "Memory CD4 T" "CD14+ Mono"   "B"            "CD8 T"       
# [6] "FCGR3A+ Mono" "NK"           "DC"           "Platelet"    
vertex.receiver = c(1,2,4,6) # define a numeric vector （淋系细胞）giving the index of the celltype as targets
#par(mar=c(5.1,4.1,4.1,2.1))
dev.off()
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy.pdf

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as TIL/CXCL_circle.pdf

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# save as TIL/CXCL_chord.pdf

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# save as TIL/CXCL_heatmap.pdf

#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对TGFb贡献最大）
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #提取对TGFb有贡献的所有配体受体 
# save as TIL/CXCL_LR_contribution.pdf

#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.TGFb[1,] 
vertex.receiver = c(1,2,4,6) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy2.pdf
dev.off()
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# save as TIL/CXCL_circle2.pdf

netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# save as TIL/CXCL_chord2.pdf

# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步
dir.create("all_pathways_com_circle") #创建文件夹保存批量画图结果
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../")

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), 
                     targets.use = c(1,2,4,6), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.png", p, width = 8, height = 12,dpi = 1000) #髓系对淋巴的调节
# save as TIL/Mye_Lymph_bubble.pdf

#比如制定CCL和CXCL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), targets.use = c(1,2,4,6), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(3,6,8), targets.use = c(1,4,5), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)
## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "TGFb")
ggsave("TGFb_GeneExpression_vln.png", p, width = 8, height = 8,dpi = 1000)
p = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot")
ggsave("TGFb_GeneExpression_dot.png", p, width = 8, height = 6,dpi = 1000)

#计算网络中心性权重
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


#通过计算每个细胞群的网络中心性指标，识别每类细胞在信号通路中的角色/作用C（发送者、接收者、调解者和影响者）

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 15, height = 6, font.size = 10)
# # save as TIL/SNA_CXCL_signalingRole.pdf

#识别细胞的信号流模式

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
ht1 + ht2
# save as TIL/SNA_SignalingPattern.pdf

#计算分解成几个因子(pattern)比较合适（这一步运行比较慢 。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
install.packages('path.package')

selectK(cellchat, pattern = "outgoing")
# save as TIL/NMF_outgoing_selectK.pdf


#热图

nPatterns = 2 # 挑选曲线中第一个出现下降的点（从2就开始下降了）
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
# save as TIL/NMF_outgoing_comPattern_heatmap.pdf


###########################################################cellchat不同组别之间的分析
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
options(stringsAsFactors = FALSE)
setwd('/home/data/caihua2/ITGB6/cellchat/wt_vs_ko')
setwd('/home/data/caihua2/ITGB6/cellchat/wt_pd1_vs_ko')
setwd('/home/data/caihua2/ITGB6/cellchat/hnsc_1')

cellchat_1_3<-cellchat
cco.list_1_3<-cco.list


bca.integrated_3<-readRDS('/home/data/caihua/ITGB6/hnsc_1.rds')
bca.integrated_3<-readRDS('/home/data/caihua/ITGB6/bca.integrated_3.rds')
bca.integrated_2<-readRDS('/home/data/caihua/ITGB6/bca.integrated_2.rds')
bca.integrated_1<-readRDS('/home/data/caihua2/ITGB6/bca.integrated_1.rds')

table(bca.integrated$new.lables)
bca.integrated<-subset(bca.integrated,new.lables%in%c('Epithelial cells','Endothelial cells','Fibroblasts',
                                                      'Macrophages','Neutrophils',"T cells","DCs"))




Idents(bca.integrated) <- 'new_group'
table(bca.integrated$new_group)
sco.wt_3 <- subset(bca.integrated, new_group%in%c('cKO'))
sco.ko_3 <- subset(bca.integrated, new_group%in%c('WT'))
sco.ko_3 <- createCellChat(sco.ko_3@assays$RNA@data, meta = sco.ko_3@meta.data, group.by = "new.lables")
sco.wt_3 <- createCellChat(sco.wt_3@assays$RNA@data, meta = sco.wt_3@meta.data, group.by = "new.lables")
#save(cco.til, cco.pbmc, file = "cco.rda")

#2.2 分析样本cco.pbmc的细胞通讯网络
#⚠️：cellchat不太稳定，identifyOverExpressedGenes常出错（不出现进度条就是出错了），重启Rstudio后再运算
cellchat <- sco.wt_3
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
sco.wt_3 <- cellchat
#saveRDS(sco.wt_2, "sco.wt_2.rds")

#2.3 分析样本cco.til的细胞通讯网络
cellchat <- sco.ko_3
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
sco.ko_3 <- cellchat
#saveRDS(sco.ko_2, "sco.ko_2.rds")


cco.list <- list(WT=sco.wt_3, KO=sco.ko_3)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p
ggsave("互作数量和强度.png", p, width = 6, height = 4,dpi = 1000)
#左图展示通讯数量之间的差异，右图展示通讯强度之间的差异。本例中信号通路强度weight值过低，导致显示时均为0（实际上有数值的，只是过小，显示为0）

#数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save as Diff_number_strength_net.pdf
#红色是case相对于control上调的，蓝色是下调的。

#数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
p<-h1+h2
p
# save as Diff_number_strength_heatmap.pdf
#case和control对比，红色是上调，蓝色是下调。


#细胞互作数量对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf
#左图是control，右图是case，可以直接对比数量变化。

#3.2 指定细胞互作数量对比网络图
table(bca.integrated$new.lables)
par(mfrow = c(1,1))
s.cell <- c("Macrophages","Epithelial cells")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

#3.3 保守和特异性信号通路的识别与可视化
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
p
ggsave("保守和特异性信号通路.png", p, width = 20, height = 12,dpi = 1000)
#左图最下面5个信号通路是case组独有的
data.raw<-as.data.frame(cellchat@var.features$WT$features.info)




#3.5 细胞信号模式对比
library(ComplexHeatmap)
#总体信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
p<-ht1 + ht2
p
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6
ggsave("Compare.png", p, width = 20, height = 12,dpi = 1000)


#输出信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6


#输入信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6


#3.6 特定信号通路的对比
#网络图
table(cco.list$KO@netP$pathways)
pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

#热图
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

#3.7 配体-受体对比分析
#气泡图展示所有配体受体对的差异
levels(cellchat@idents$joint)
library(CellChat)
p <- netVisual_bubble(cellchat, sources.use = c(5),
                      targets.use = c(1:7),  comparison = c(1, 2), angle.x = 45)#signaling = c('TGFb')
p
p <- netVisual_bubble(cellchat_1_3, sources.use = c(2), targets.use = c(1,7,8,10),  comparison = c(1, 2), angle.x = 45)
p

ggsave("信号通路.png", p, width = 20, height = 18)
#case和control配对的图，文章中常见。cellphoneDB找到的结果，也可以用这种方式呈现。

p1 <- netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(1,7,8,10), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(1,7,8,10), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
pc
ggsave("信号通路上下调分开.png", pc, width = 20, height = 18)

saveRDS(cco.list,'cco.list.rds')
saveRDS(cellchat,'cellchat.rds')