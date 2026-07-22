rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)
#revise
plan(multisession, workers = 10)  
options(future.globals.maxSize = 100000 * 1024^4)
setwd("Epithelial/") 
dir.create("CNV")
setwd("CNV")
set.seed(12345)
sce=readRDS( "../sce.all_int.rds")


resolution1 <- "SCT_snn_res.0.1"
DimPlot(sce, reduction = "umap",
        group.by = resolution1,label = T) #+
sce <- sce[,!(sce$SCT_snn_res.0.1 %in% c(6))]

# ####### 提取名为 'umap' 的降维结果的坐标
# umap_coords <- Embeddings(sce, reduction = "umap")
# #####删除“异常”#########
# id1 <- rownames(umap_coords)[umap_coords[,1] > -2.5]
# #id2 <- rownames(umap_coords)[umap_coords[,2] < -1 & umap_coords[,1] < -5]
# sce2 <- sce[,colnames(sce) %in% c(id1)]

p1 <- DimPlot(sce, group.by = "SCT_snn_res.0.1",label = T)
p1
table(sce$SCT_snn_res.0.1)
#sce <- sce2[,!(sce2$SCT_snn_res.0.1 %in% c(7,8))]
#p1 <- DimPlot(sce, group.by = "SCT_snn_res.0.1",label = T)
#p1
#table(sce$SCT_snn_res.0.1)
#unique(sce$SCT_snn_res.0.1)
sce$SCT_snn_res.0.1 <- as.numeric(as.character(sce$SCT_snn_res.0.1))
sce$SCT_snn_res.0.1[sce$SCT_snn_res.0.1 == 7] <- 6
unique(sce$SCT_snn_res.0.1)

p1 <- DimPlot(sce, group.by = "SCT_snn_res.0.1",label = T)
p1
sce$SCT_snn_res.0.1 <- factor(sce$SCT_snn_res.0.1, levels = c(0,1,2,3,4,5,6))
unique(sce$SCT_snn_res.0.1)
resolution1 <- "SCT_snn_res.0.1"
Idents(sce) <- resolution1

epi <- sce

umap =DimPlot(epi, reduction = "umap",
              group.by = resolution1,label = T) #+
umap
ggsave(paste0(resolution1,'_umap.pdf'),width = 9,height = 7)

cancer1 <- c("ITGA6", "KRT1", "UCHL1", "SERPINE1", "CCND1")
#Single-cell profiling of response to neoadjuvant chemo-immunotherapy in surgically resectable esophageal squamous cell carcinoma
DotPlot(epi, features = cancer1,group.by = resolution1) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,       # 调整字体大小
      vjust = 0.5     # 垂直对齐微调
    ),
    plot.margin = margin(10, 10, 10, 30)  # 增加底部边距（防止标签被裁剪）
  )
# | 状态         | 标志组合                                            |
#   | ---------- | ----------------------------------------------- |
#   | **正常上皮**   | EpCAM+, Lgr5+/Alpi+/Muc2+, Mki67-               |
#   | **癌前/可塑性** | EpCAM+, Krt8+/Krt5+, Gprc5a-/-背景                |
#   | **癌上皮**    | EpCAM+, Krt8+/Krt18+, Mki67+, Ccnd1+, Serpine1+ |
#   
# cancer1 <- c("EPCAM", "TSTA3", "NECTIN4", "SFN", 
#               "KRT5", "KRT16", "KRT18", "KRT17", "KRT19", "KRT6A", "KRT6B")
# DotPlot(epi, features = cancer1,group.by = resolution1) +
#   theme(
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1,
#       size = 8,       # 调整字体大小
#       vjust = 0.5     # 垂直对齐微调
#     ),
#     plot.margin = margin(10, 10, 10, 30)  # 增加底部边距（防止标签被裁剪）
#   )
###CNV分析确证####
#1.准备输入文件开始分析之前，需要准备三个文件，https://github.com/broadinstitute/inferCNV/wiki/File-Definitions
#1.1表达矩阵文件。counts矩阵，每列为一个细胞，每行为一个基因。
library(infercnv)
dir.create("CNV")
setwd("CNV") 


###吃亏上当一百年，如果不加前缀，后面0-20直接变成了1-21
load("../../../T_NK//finaltype/t_meta.Rdata")
tcell <- rownames(t_meta)
tcell <- sample(tcell, size = 2000, replace = FALSE)

load("../../../Bcell/finaltype/t_meta.Rdata")
bcell <- rownames(t_meta)
bcell <- sample(bcell, size = 2000, replace = FALSE)

sce <- readRDS("../../../../3-auto_annotate/sce.clean.rds")
epi_meta <- data.frame(epi@meta.data)
epi_meta$rowname <- rownames(epi_meta)
# epi_meta2 <-  epi_meta %>%
#   group_by(SCT_snn_res.0.1) %>%
#   sample_n(size = min(500, n()), replace = FALSE) %>%
#   ungroup()
epi_meta2 <- epi_meta
epi_meta2$SCT_snn_res.0.1 <- paste0("C",epi_meta2$SCT_snn_res.0.1)
epi_n <- epi_meta2$rowname
sce <- sce[,colnames(sce) %in% c(tcell,bcell,epi_n)]
sce <- sce[rowSums(sce) > 0,]
dim(sce)

groupinfo=data.frame(v1= c(tcell,bcell,epi_n),
                     v2=c(rep('tcell',length(tcell)),
                           rep('bcell',length(bcell)),
                     epi_meta2$SCT_snn_res.0.1))
rownames(groupinfo) <- groupinfo$v1
groupinfo <- groupinfo[colnames(sce),]
sum(colnames(sce) != groupinfo$v1)
counts <- GetAssayData(sce, layer = "counts")



library(AnnoProbe)
geneInfor=annoGene(rownames(counts),"SYMBOL",'human')
geneInfor <- geneInfor %>% select(SYMBOL,chr,start,end)
geneInfor <- geneInfor[geneInfor$SYMBOL %in% rownames(counts),]

z <- as.data.frame(table(geneInfor$SYMBOL))
z <- z$Var1[z$Freq > 1]
geneInfor <- geneInfor[!duplicated(geneInfor$SYMBOL), ]
z %in% geneInfor$SYMBOL
unique(geneInfor$chr)
## 1) remove mitochondrial genes (names starting with "MT-")
keep <- !grepl("^MT-|^mt-", rownames(geneInfor))
# 2) remove ribosomal genes (RPS/RPL) and hemoglobin genes
keep <- keep & !grepl("^RPS|^RPL|^HBA|^HBB", rownames(geneInfor), ignore.case = TRUE)
keep <- keep & !(geneInfor$chr %in% c("X","Y","chrM"))
geneInfor <- geneInfor[keep,]
unique(geneInfor$chr)
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
counts <- counts[rownames(counts) %in% geneInfor$SYMBOL,]
counts <- counts[match(geneInfor$SYMBOL,rownames(counts)),]

library(infercnv)
rownames(groupinfo) <- groupinfo$v1
groupinfo <- groupinfo %>% dplyr::select(-v1)
rownames(geneInfor) <- geneInfor$SYMBOL
geneInfor <- geneInfor %>% dplyr::select(-SYMBOL)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts),
                                    annotations_file=groupinfo,
                                    delim="\t",
                                    gene_order_file= geneInfor,
                                    ref_group_names=c('tcell',
                                                      'bcell'))  ## 这个取决于自己的分组信息里面的
unlink("infercnv_output/",recursive = T)
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=TRUE ,   # cluster
                              #write_expr_matrix = T, #评分文件
                              output_format = "pdf",
                              num_threads = 5,
                              hclust_method="ward.D2", 
                              plot_steps=TRUE)
# infercnv::plot_cnv(
#   infercnv_obj2,
#   output_filename = "infercnv_recolor",
#   output_format   = "pdf",
#   plot_chr_scale  = TRUE,
#   # 关键：包成函数，而不是直接给向量
#   custom_color_pal = colorRampPalette(c("#053061", "white", "#67001F"))
# )
#save(groupinfo,geneInfor,dat,epi_sce,file = 'infercnv.Rdata')

#应用infercnv的结果
infer_CNV_obj<-readRDS('infercnv_output/run.final.infercnv_obj')
#library(RColorBrewer)
# infercnv::plot_cnv(infer_CNV_obj, #上两步得到的infercnv对象
#                    plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
#                    output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
#                    custom_color_pal = color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色

expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)

meta = sce@meta.data
###如果从CNV图实在分不清，也可以看看CNVscore
#run.final.infercnv_obj对象文件

tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$tcell]
tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$bcell]
tmp= cbind(tmp1,tmp2)
down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd)) ## 定义拷贝数减少的阈值
up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd)) ### 定义拷贝数增加的阈值
oneCopy=up-down
#正态分布：假设基因表达数据在正常细胞中服从正态分布，那么大部分（约95%）的正常细胞的表达值会落在平均值加减两倍标准差的范围内。
#阈值设定：通过设定 down 和 up 阈值，可以区分出显著低于或高于正常范围的表达值，从而推断拷贝数的减少或增加。
#单位拷贝变化：oneCopy 的计算是为了量化拷贝数变化的程度，便于后续对拷贝数变异进行分类和评分
oneCopy
a1= down- 2*oneCopy
a2= down- 1*oneCopy
down;up
a3= up +  1*oneCopy
a4= up + 2*oneCopy 
  
cnv_score_table<-infer_CNV_obj@expr.data
cnv_score_table[1:4,1:4]
cnv_score_mat <- as.matrix(cnv_score_table)
  
# Scoring
cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
# Check
table(cnv_score_table[,1])

# Replace with score 
cnv_score_table_pts <- cnv_score_mat
rm(cnv_score_mat)
# 
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
cnv_score_table_pts[1:4,1:4]
str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
colnames(cell_scores_CNV) <- "cnv_score" 
head(cell_scores_CNV) 
score=cell_scores_CNV
sum(rownames(score) != rownames(groupinfo))

score$group = groupinfo$v2
colnames(score)
library(ggthemes)
library(ggplot2)
score$group <- gsub("tcell","Tcell",score$group)
score$group <- gsub("bcell","Bcell",score$group)
p <- ggplot(score, aes(x = group, y = cnv_score, fill = group)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  labs(x = "", y = "CNV Score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1))
p

ggsave(p, file = "CNV_score.pdf", width = length(unique(score$group)) / 1.5, height = 5)
###推测有可能是正常上皮
data1 <- score %>%
  group_by(group) %>%                # 按分组列分组
  summarise(median_value = median(cnv_score)) %>%  # 计算每组中位数
  mutate(is_gt_1500 = median_value > 1000)    # 判断中位数是否大于1500
idt <- data1$group[data1$is_gt_1500]
setwd('../')
#####细胞生物学命名
length1 <- length(unique(epi$SCT_snn_res.0.1))
celltype=data.frame(ClusterID= 0:c(length1 - 1),
                    celltype=  0:c(length1 - 1)) 
# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
#idt <- gsub("C","",idt)

#推测3，5，6，7，8为正常细胞
epi$subtype <-  as.character(epi$SCT_snn_res.0.1)
epi@meta.data <- epi@meta.data %>%
  mutate(subtype = ifelse(subtype %in% c("0","1","2"), "Normal cell", "Malignant cell"))
umap =DimPlot(epi, reduction = "umap",
              group.by = "subtype",label = T) +
  labs(title = "")
umap
ggsave(umap,file = paste0('subtype.pdf'),width = 9,height = 7)
##############样品比较#############################
data1 <- data.frame(table(epi$sample,epi$subtype))
data1 <- data1 %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = 100 * Freq / `sum(Freq)`)

meta <- read.delim("../../../meta.txt")

data2 <- data1 %>% left_join(meta,by = c("Var1" = "id"))

dir.create("compare2")
setwd("compare2/")
order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
data2$group <- factor(data2$group, levels = order1)
color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")

data3 <- data2
data3 <- data3 %>%
  arrange(number, group)
#data3是用来配对的
comp1 <- list(c("Post_R","Pre_R"),c("Post_NR","Pre_NR"))
comp2 <- list(c("Post_R","Post_NR"),c("Pre_R","Pre_NR"))
library(ggpubr)
p1 <- ggplot(data3, aes(x = group, y = percent)) +
  geom_boxplot(aes(color = group), show.legend = F, 
               width = 0.6,outlier.colour = "black",
               outlier.size = 0.5) +  #箱线图
  scale_fill_manual(values = color1) +  #设置颜色
  #geom_point(size = 3,color='red') +
  geom_point(size = 0.5) +  #绘制散点
  geom_line(aes(group = number), color = 'gray', lwd = 0.1) +
  facet_wrap( ~ Var2,scales = "free_y",ncol = 3) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        panel.background = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'Relative abundance (%)')+
  stat_compare_means(paired = T,
                     method = "wilcox.test",comparisons = comp1,vjust = 1.5)
p1
n1 <- length(unique(data3$Var2))
n1 <- sqrt(n1)
ggsave(p1, filename = paste0("four_paired.png"),width = n1 * 2 + 4,height = n1 * 3)
ggsave(p1, filename = paste0("four_paired.pdf"),width = n1 * 2 + 4,height = n1 * 3)

comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
p1 <- ggplot(data2, aes(x = group, y = percent)) +
  geom_boxplot(aes(color = group), show.legend = F, 
               width = 0.6,outlier.colour = "black",
               outlier.size = 0.5) +  #箱线图
  scale_fill_manual(values = color1) +  #设置颜色
  #geom_point(size = 3,color='red') + 
  geom_point(size = 0.5) +  #绘制散点
  geom_line(aes(group = number), color = 'gray', lwd = 0.1) +
  facet_wrap( ~ Var2,scales = "free_y",ncol = 3) +
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'Relative abundance (%)')+
  stat_compare_means(paired = F,
                     method = "wilcox.test",comparisons = comp1,vjust = 1.5)
p1
n1 <- length(unique(data2$Var2))
n1 <- sqrt(n1)
ggsave(p1, filename = paste0("four_unpaired.png"),width = n1 * 2 + 4,height = n1 * 3)
ggsave(p1, filename = paste0("four_unpaired.pdf"),width = n1 * 2 + 4,height = n1 * 3)

# 假设你的数据框名为 df
data4 <- data2 %>%
  group_by(group, Var2) %>%
  summarise(percent = mean(percent, na.rm = TRUE))
# 2. 绘制堆积柱状图
p1 <- ggplot(data4, aes(x = group, y = percent, fill = Var2)) +
  geom_col(position = "stack") +  # position = "stack" 是堆积柱状图的关键
  labs(
    #title = "各实验组细胞类型比例分布",
    #x = "实验组",
    #y = "细胞比例 (%)",
    fill = "Celltype"  # 修改图例标题
  ) +
  #scale_fill_brewer(palette = "Set2") +  # 使用更美观的颜色方案
  theme(
    panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 0.5)
  )
ggsave(p1, filename = paste0("four_group.png"),width = 6,height = 4.5)

#####pdcd1 picture########################
pd1 <- c("CD274","PDCD1LG2")
sce.all2 <- epi
sce.all2$group <- meta$group[match(sce.all2$orig.ident,meta$id)]
DotPlot(sce.all2,features = "PDCD1",group.by = "group")
group1 <- unique(meta$group)
for (k in pd1) {
  result1  <- list()
  result2  <- list()
  for (i in group1) {
    print(i)
    sce1 <- sce.all2[,sce.all2$group == i]
    v1 <-  AverageExpression(sce1, 
                             group.by = "subtype",
                             assays = "SCT")  #
    v1 <- v1$SCT
    v1 <- as.matrix(v1[c(k,"CD8A"),])
    v1 <- as.data.frame(v1[1,,drop = F])
    rownames(v1)[1] <- i
    
    colnames(v1)[colnames(v1) == "T-NK"] <- "T_NK"
    result1[[i]] <- v1
    
    v2 <- v1
    colnames(v1)
    for (j in colnames(v1)) {
      print(j)
      if (j == "T-NK"){
        j <- "T_NK"
      }
      if (v1[,j] == 0){
        b <- 0
      }else{
        sce2 <- sce1[,sce1$subtype == j]
        b <- sum(sce2@assays$SCT$data[k,] > 0) / dim(sce2)[2]
        b <- b * 100
      }
      v2[1,j] <- b
    }
    result2[[i]]  <- v2
  }
  abundance <- bind_rows(result1)
  percent <- bind_rows(result2)
  
  library(reshape2)
  long_df1 <- melt(as.matrix(abundance), 
                   varnames = c("group", "celltype"),
                   value.name = "expression")
  #long_df1$celltype <- recode(long_df1$celltype, "T-NK" = "T_NK")
  long_df2 <- melt(as.matrix(percent), 
                   varnames = c("group", "celltype"),
                   value.name = "percent")
  data1 <- merge(long_df1, long_df2, 
                 by = c("group", "celltype"),
                 all = TRUE)  # all=TRUE表示全连接，保留所有行
  # 假设您的数据框叫df，包含group, celltype, expression, percent四列
  
  #data1$celltype[data1$celltype == "Tcell"] <- "T/NK"
  #data1$celltype <- recode(data1$celltype, "Tcell" = "T/NK")
  data1$group <- factor(data1$group, levels = rev(c("Pre_NR","Post_NR","Pre_R","Post_R")))
  p1 <- ggplot(data1, aes(x = celltype, y = group)) +
    geom_point(aes(size = percent, color = expression)) +
    scale_size_continuous(range = c(1, 10), name = "Expression\nPercentage") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Expression\nLevel") +
    labs(x = "", y = "", 
         title = paste0(k," expression")) +
    #theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 横坐标字体调大
          axis.text.y = element_text(size = 12),  # 纵坐标字体调大
          panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
          legend.key = element_rect(fill = "transparent", color = NA)  # 图例键背景透明，无边框
    )
  ggsave(p1,file = paste0(k," pdcd1.png"),width = 8, height = 5)
}
#############################


###################marker################################################
setwd("../")
t_meta <- epi@meta.data
save(t_meta,file = "t_meta.Rdata")

#   | **正常上皮**   | EpCAM+, Lgr5+/Alpi+/Muc2+, Mki67-               |
#   | **癌前/可塑性** | EpCAM+, Krt8+/Krt5+, Gprc5a-/-背景                |
#   | **癌上皮**    | EpCAM+, Krt8+/Krt18+, Mki67+, Ccnd1+, Serpine1+ |

Cancer <- c(  "MKI67",     # Ki-67，增殖标志物
              "PCNA",      # 增殖细胞核抗原
              "CCND1",     # Cyclin D1，细胞周期失调
              "CDKN2A",    # p16INK4a，细胞周期失调
              "TP53",      # p53（突变型），抑癌基因突变
              "CTNNB1",    # β-catenin，Wnt通路激活
              "BCL2",      # 抗凋亡蛋白
              "BIRC5",     # Survivin，抗凋亡蛋白
              "SERPINE1",  # PAI-1，侵袭和转移
              "MMP2",      # 基质金属蛋白酶2，侵袭
              "MMP9",      # 基质金属蛋白酶9，侵袭
              "VEGFA",     # 血管内皮生长因子A，血管生成
              "SOX2",      # 干细胞因子，去分化
              "POU5F1",    # OCT4，干细胞因子
)
Normal <- c(  "CDH1",      # E-cadherin，细胞粘附
              "DSP",       # Desmoplakin，细胞连接
              "TJP1",      # ZO-1，紧密连接
              "ITGA6",     # Integrin α6，基底膜连接
              "ITGB4",     # Integrin β4，基底膜连接
              "LGR5",      # 肠道干细胞标志物
              "MUC2",      # 杯状细胞分泌蛋白
              "ALPI",      # 肠道碱性磷酸酶
              "VILL",      # Villin，刷状缘蛋白
              "KRT5",      # 角蛋白5，复层上皮基底细胞
              "KRT14"      # 角蛋白14，复层上皮基底细胞
              )
marker1 <- list(Cancer = Cancer,
                Normal = Normal)
p1 <- DotPlot(epi, features = marker1,group.by = "subtype") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,       # 调整字体大小
      vjust = 0.5     # 垂直对齐微调
    ),
    plot.margin = margin(10, 10, 10, 30)  # 增加底部边距（防止标签被裁剪）
  )
p1
p1 <- DotPlot(epi, features = marker1,group.by = resolution1) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,       # 调整字体大小
      vjust = 0.5     # 垂直对齐微调
    ),
    plot.margin = margin(10, 10, 10, 30)  # 增加底部边距（防止标签被裁剪）
  )
ggsave(p1, file = "final_marker.pdf",height = 4, width = 10) 

############final infercnv#########################################
# groupinfo$v2[groupinfo$v2 %in% paste0("C",c(0:2,4:7))] <- "Malignant cell"
# groupinfo$v2[groupinfo$v2 %in% paste0("C",c(3))] <- "Normal cell"
# groupinfo$v2[groupinfo$v2 %in% ("bcell")] <- "Bcell"
# groupinfo$v2[groupinfo$v2 %in% ("tcell")] <- "Tcell"
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
#                                     annotations_file=groupinfo,
#                                     delim="\t",
#                                     gene_order_file= geneInfor,
#                                     ref_group_names=c('Tcell',
#                                                       'Bcell'))  ## 这个取决于自己的分组信息里面的
# unlink("infercnv_output2/",recursive = T)
# # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                               out_dir= "infercnv_output2",  # dir is auto-created for storing outputs
#                               cluster_by_groups=TRUE ,   # cluster
#                               #write_expr_matrix = T, #评分文件
#                               output_format = "pdf",
#                               num_threads = 5,
#                               hclust_method="ward.D2", 
#                               plot_steps=TRUE)