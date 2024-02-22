
##建立工作文件夹##
dir.create("/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB")
setwd("/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB")
getwd()

##加载R包##
library(biomaRt)
library(tidyverse)
library(DT)

##读取数据##
bca.integrated <- readRDS("/home/data/maomao/bca_anticd276_matrix/20220418/bca.integrated.rds")
table(bca.integrated$new.lables)
DefaultAssay(bca.integrated) <- "RNA"
Epithelial <- subset(bca.integrated, new.lables=="Epithelial cells")
Macrophages <- subset(bca.integrated, new.lables=="Macrophages")
dim(Macrophages)
dim(Epithelial)
Epithelial$new.lables <- ifelse(Epithelial$new_group %in% c("cKO"),"Epithelial_WT","Epithelial_cKO") 
Macrophages$new.lables <- ifelse(Macrophages$new_group %in% c("cKO"),"Macrophages_WT","Macrophages_cKO")
merge_object<-merge(Epithelial,Macrophages)
dim(merge_object)
table(merge_object$new.lables)

#从源头转小鼠基因到人类基因##转大表格##
human = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl(biomart="ensembl", dataset = "mmusculus_gene_ensembl")
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('csv', 'excel', 'pdf'),
                               lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
}

##读取features.tsv文件(这里的文件指的单细胞测序回来的数据) #同一批样本中的不同样本的这个基因对照是一样的
dat<-merge_object@assays$RNA@meta.features
gene_list <- as.data.frame(x = rownames(dat) , y= 1)
write.csv(gene_list,"/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDBgene_list.csv")
gene_list<-read.csv("/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDBgene_list.csv")
create_dt(gene_list)
dim(gene_list)

##转小鼠基因名为人类基因名
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene_list$`rownames(dat)`, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
genesV2 <- genesV2[!duplicated(genesV2[,2]),] 
create_dt(genesV2)
dim(genesV2)

##保存文件
getwd()
write.table(genesV2, "/home/data/maomao/bca_anticd276_matrix/harmony/cellphoneDB/allcrosstalk/mouse_to_human_genes.txt", sep="\t", row.names=F, quote=F)

install.packages('homologene')
#加载homologene这个R包
library(homologene)
#这里以小鼠的三个基因为例
#更多基因方法是一样的
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
genesV2<-homologene(gene_list$rownames.dat., inTax = 10090, outTax = 9606)

genesV2<-genesV2[,c(1,2)]


#应用转换后的小鼠到人类基因表
##读取文件
#genesV2 <- read.table("/home/data/maomao/bca_anticd276_matrix/harmony/cellphoneDB/allcrosstalk/mouse_to_human_genes.txt", sep="\t", header=T)
dim(genesV2)

### 准备小鼠数据集
library(Seurat)
installed.packages("SeuratData")
install.packages("SeuratData")
library(SeuratData)


##读取自己的小鼠数据##转完大表格后从这一步开始##
#seurat_object <- merge_object
#sub_cells <- subset(seurat_object, cells = sample(Cells(seurat_object), 100))
#sub_cells@assays$RNA[1:4,1:4]

## Extract Expression Data

sp_all_counts <- as.matrix(merge_object@assays$RNA@data) 
dim(merge_object)
dim(sp_all_counts)
sp_all_counts <- data.frame(gene=rownames(sp_all_counts), sp_all_counts, check.names = F)
dim(merge_object)
dim(sp_all_counts)


sp_all_counts[1:4, 1:4]

# 转小鼠基因名为人类基因名
sp_all_counts$Gene <- genesV2[match(sp_all_counts$gene, genesV2[,1]),2]
dim(merge_object)
dim(sp_all_counts)
#sp_all_counts <- subset(sp_all_counts, Gene!='NA')
sp_all_counts <- dplyr::select(sp_all_counts, Gene, everything())
dim(merge_object)
dim(sp_all_counts)
sp_all_counts <- sp_all_counts[, !(colnames(sp_all_counts) %in% 'gene')]
dim(merge_object)
dim(sp_all_counts)
sp_all_counts[1:4,1:4]



## 保存文件
#write.table(sp_all_counts, "/home/data/maomao/maomao/cellphoneDB/merge_object_counts_human.txt", row.names=F, sep='\t', quote=F)


##输出cellphonedb需要的文件##这一步不用修改名称
dim(merge_object)
MC_Epi_counts <- as.matrix(merge_object@assays$RNA@data)
MC_Epi_counts <- data.frame(Gene=sp_all_counts$Gene, MC_Epi_counts)
dim(MC_Epi_counts)

MC_Epi_meta <- data.frame(Cell=rownames(merge_object@meta.data), cell_type=merge_object@meta.data$new.lables)
dim(MC_Epi_meta)
#head(MC_Epi_counts)
write.table(MC_Epi_counts, "/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB/merge_object_counts.txt", row.names=F, quote = FALSE,sep='\t')
write.table(MC_Epi_meta, "/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB/merge_object_meta.txt", row.names=F, quote = FALSE,sep='\t')
getwd()
#MC_Epi_meta<-read.table("/home/data/maomao/maomao/cellphoneDB/merge_object_meta.txt")
#cellphonedb method statistical-analysis merge_object_meta.txt merge_object_counts.txt --threads=50 --counts-data=gene_name
#cellphonedb plot heatmap-plot merge_object_meta.txt
#cellphonedb可视化
##读取cellphonedb结果文件
setwd("/home/data/maomao/bca_anticd276_matrix/harmony/cellphoneDB/allcrosstalk")
new_path <- 'out'
mypvals <- read.delim(file.path(new_path,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(file.path(new_path, "means.txt"), check.names = FALSE)
mysigmeans <- read.delim(file.path(new_path, "significant_means.txt"), check.names = FALSE)
dim(mymeans)
mymeans[1:4, 1:4]

#调整受体配体对顺序，要求配体在前，受体在后
##Rearrange data column sequence
library(dplyr)
order_sequence <- function(df){
  da <- data.frame()
  for(i in 1:length(df$gene_a)){
    sub_data <- df[i, ]
    if(sub_data$receptor_b=='False'){
      if(sub_data$receptor_a=='True'){
        old_names <- colnames(sub_data)
        my_list <- strsplit(old_names[-c(1:11)], split="\\|")
        my_character <- paste(sapply(my_list, '[[', 2L), sapply(my_list, '[[', 1L), sep='|')
        new_names <- c(names(sub_data)[1:4], 'gene_b', 'gene_a', 'secreted', 'receptor_b', 'receptor_a', "annotation_strategy", "is_integrin", my_character)
        sub_data = dplyr::select(sub_data, new_names)
        # print('Change sequence!!!')
        names(sub_data) <- old_names
        da = rbind(da, sub_data) 
      }
    }else{
      da = rbind(da, sub_data)
    }
  }
  return(da)
}
# 受体配体对，要求一个为受体，一个为配体
df <- subset(mymeans, receptor_a == 'True' & receptor_b == 'False' | receptor_a == 'False' & receptor_b == 'True')

# 筛选单基因的受体配体对，待商榷！
df <- df %>% dplyr::mutate(na_count = rowSums(is.na(df) | df == "")) %>% subset(na_count == 0) %>% dplyr::select(-na_count)
dim(df)


# 运行函数，重排受体配体顺序，并生成新列
means_order <- order_sequence(df) %>% tidyr::unite(Pairs, gene_a, gene_b)
pvals_order <- order_sequence(mypvals) %>% tidyr::unite(Pairs, gene_a, gene_b)

# 保存文件
# write.table(means_order %>% dplyr::select(-interacting_pair), file.path(new_path,"means_order.txt"), sep = '\t', quote=F, row.names=F)
# write.table(pvals_order %>% dplyr::select(-interacting_pair), file.path(new_path,"pvalues_order.txt"), sep = '\t', quote=F, row.names=F)

##合并表达量文件和pvalue文件
means_sub <- means_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
pvals_sub <- pvals_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
means_gather <- tidyr::gather(means_sub, celltype, mean_expression, names(means_sub)[-1])
pvals_gather <- tidyr::gather(pvals_sub, celltype, pval, names(pvals_sub)[-1])
mean_pval <- dplyr::left_join(means_gather, pvals_gather, by = c('Pairs', 'celltype'))
create_dt(mean_pval)

#保存文件
# write.table(mean_pval, file.path(new_path,"mean_pval_order.txt"), sep = '\t', quote=F, row.names=F)

##提取显著性表达的受体配体对
# 至少在一组细胞类型两两组合中，pvalue显著的受体配体对，认为是显著性表达的受体配体对
a <- mean_pval %>% dplyr::select(Pairs, celltype, pval) %>% tidyr::spread(key=celltype, value=pval)
sig_pairs <- a[which(rowSums(a<=0.05)!=0), ]
dim(sig_pairs)
## [1] 72 17

# 保存显著性表达的受体配体对
mean_pval_sub <- subset(mean_pval, Pairs %in% sig_pairs$Pairs)
# mean_pval_sub_pic <- mean_pval_sub
# mean_pval_sub_pic <- filter(mean_pval_sub_pic,mean_pval_sub_pic$celltype==
#                               c("CD4+ effector cells|Epithelial cells","CD8+ effector cells|Epithelial cells","NKT|Epithelial cells",
#                                 "Treg|Epithelial cells","Ki67+ T cells|Epithelial cells","Th17|Epithelial cells",
#                                 "rδT|Epithelial cells","Navie T cells|Epithelial cells","ILC2|Epithelial cells","ILC1|Epithelial cells",
#                                 "Epithelial cells|CD4+ effector cells","Epithelial cells|CD8+ effector cells","Epithelial cells|NKT",
#                                 "Epithelial cells|Treg","Epithelial cells|Ki67+ T cells","Epithelial cells|Th17",
#                                 "Epithelial cells|rδT","Epithelial cells|Navie T cells","Epithelial cells|ILC2","Epithelial cells|ILC1"))



#保存文件
write.table(mean_pval_sub, file.path("/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB/","mean_pval_sig_subtumor.txt"), sep = '\t', quote=F, row.names=F)
mean_pval_sub <- read.table("/home/data/maomao/bca_anticd276_matrix/20220418/cellphoneDB/cellphoneDB1.txt", sep="\t", header=T)


###可视化显著性表达的受体配体对
#library(RColorBrewer)
#library(scales)
#library(ggplot2)
#library(cowplot)

# 比较均值和P值数据变换前后的分布
p1 <- mean_pval_sub %>% ggplot(aes(x=mean_expression)) + geom_density(alpha=0.2, color='red')
p2 <- mean_pval_sub %>% ggplot(aes(x=log2(mean_expression))) + geom_density(alpha=0.2, color='red')
mean_distribution <- plot_grid(p1, p2, labels = "AUTO", nrow=1)

p1 <- mean_pval_sub %>% ggplot(aes(x=pval)) + geom_density(alpha=0.2, color='red')
p2 <- mean_pval_sub %>% ggplot(aes(x=-log10(pval+1*10^-3))) + geom_density(alpha=0.2, color='red')
pval_distribution <- plot_grid(p1, p2, labels = "AUTO", nrow=1)
plot_grid(mean_distribution, pval_distribution, labels = "AUTO", ncol=1)

# 展示部分（10对）显著性表达的受体配体对
p <- mean_pval_sub %>% dplyr::arrange(Pairs) %>% head(100*100) %>% 
  ggplot(aes(y=Pairs, x=celltype )) +
  #geom_point(aes(color=log2(mean_expression), size=pval)) +
  # scale_size(trans = 'reverse') +
  geom_point(aes(color=log2(mean_expression), size=-log10(pval+1*10^-3)) ) +
  guides(colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
  labs(x='', y='') +
  scale_color_gradientn(name='Expression level \n(log2 mean expression \nmolecule1, molecule2)', colours = rainbow(100)) +
  #scale_color_gradient2('Expression level \n(log2 mean expression \nmolecule1, molecule2)', low = 'blue' , mid = 'yellow', high = "red") +
  theme(axis.text.x= element_text(angle=45, hjust=1)) +
  # coord_flip() +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background   = element_blank(),
    axis.title.x       = element_blank(),
    axis.title.y       = element_blank(),
    axis.ticks         = element_blank()
    # plot.title         = element_text(hjust = 0.5),
    # legend.position = 'bottom' # guides(fill = guide_legend(label.position = "bottom"))
    # legend.position    = "bottom"
    # axis.text.y.right  = element_text(angle=270, hjust=0.5)
  ) +
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        # legend.key.height = unit(1, 'cm'), #change legend key height
        # legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size

p
DimPlot(Immune,reduction = "umap",label = F,pt.size = 0.5,split.by = "new_group")




DefaultAssay(Epithelial) <- "RNA"
DefaultAssay(Immune) <- "RNA"
Idents(Epithelial) <- "new_group"
Idents(Immune) <- "new_group"
VlnPlot(Immune, features = c("C1qb","Arg1","Top2a","Mki67","Trac","Nkg7","Lgkc","Ccr7","Ccl17"),pt.size = 0.0)
VlnPlot(Epithelial, features = c("Anxa1"),pt.size = 0.0,ncol = 1,y.max = 10.0 )+stat_compare_means(comparisons = list( c("WT", "KO")),label = "p.signif",method = "t.test")
?VlnPlot
VlnPlot(Immune, features = c("Cd74"),pt.size = 0.0,ncol = 1,y.max =12.0 )+stat_compare_means(comparisons = list( c("WT", "KO")),label = "p.signif",method = "t.test")



table(Immune$new.lables)
Macrophage1 <- subset(Immune, new.lables=="Arg1+ Macrophage")
Macrophage2 <- subset(Immune, new.lables=="Cycle Macrophage")
Macrophage <- merge(Macrophage1,Macrophage2)
VlnPlot( Macrophage2, features = c("Cd74"),pt.size = 0.0,ncol = 1,y.max = 10.0 )+stat_compare_means(comparisons = list( c("DMSO", "antiCD276")),label = "p.signif",method = "t.test")
