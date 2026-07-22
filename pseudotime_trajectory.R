library(Seurat)
library(ggplot2)
library(harmony)
library(tidyverse)
library(clustree)
library(cowplot)
library(data.table)
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(tidyverse)
setwd("T_NK/finaltype/tex/")
dir.create("pseudotime")
setwd("pseudotime/")
library(monocle3)
#######mast##
sce <-  readRDS("../sce.all_int.rds")
load("../t_meta.Rdata")
unique(t_meta$subtype)

sce <- sce[,rownames(t_meta)]
a <- sum(rownames(t_meta) != colnames(sce))
if (a > 0){
  print("something was wrong")
}
sce@meta.data <- t_meta
p1 <- DimPlot(sce,group.by = "subtype",label = T)
p1

meta <- read.delim("../../../../../meta.txt")
sce$group <- meta$effect[match(sce$orig.ident,meta$id)]
p1 <- DimPlot(sce,group.by = "group",label = F)
p1


sce1 <- JoinLayers(sce,assay = "RNA")

expression_matrix <- sce1[["RNA"]]$counts
cell_metadata <- sce1@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
cds <- preprocess_cds(cds, num_dim = 20)

plot_pc_variance_explained(cds)  #observe whether num_dim is enough

# 或者手动过滤低表达基因
# expressed_genes <- rowSums(counts(cds) > 0) > 10  # 至少在10个细胞中表达
# cds <- cds[expressed_genes, ]
#unique(cds@colData$sample)
cds <- align_cds(cds, alignment_group = "sample") #adjust batch 

cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
# 获取Seurat的UMAP坐标
reducedDims(cds)$UMAP <- Embeddings(sce1, "umap")
# 使用Seurat的聚类结果
cds@clusters$UMAP$clusters <- sce1$subtype
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "subtype")

#2.1 Monocyte-like (early / progenitor-like) macrophages
#These are the most common roots
#Human
#c("S100A8", "S100A9", "FCN1", "VCAN", "LYZ", "CTSS", "LGALS3")
#mouse
#c("S100a8", "S100a9", "Ly6c2", "Lyz2", "Ctss", "Fcer1g")

# ciliated_genes <- c("S100A8", "S100A9", "FCN1", "VCAN", "LYZ", "CTSS", "LGALS3")
# 
# 
# plot_cells(cds,
#            genes=ciliated_genes,
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE) + 
#   scale_color_gradientn(
#     colours = c("blue", "yellow", "red"),
#     values = scales::rescale(c(0, 0.1, 0.6, 1)),  # 调整颜色分布位置
#     name = "Expression"
#   )
#从UMAP图中看不出来
# 1. 提取表达数据
# gene_exp_data <- FetchData(sce, vars = c("subtype","S100A8", "S100A9", "FCN1", "VCAN", "LYZ", "CTSS", "LGALS3"))
# # 2. 转换数据为长格式，便于ggplot2绘图
# exp_long <- pivot_longer(gene_exp_data,
#                          cols = c("S100A8", "S100A9", "FCN1", "VCAN", "LYZ", "CTSS", "LGALS3"),
#                          names_to = "Gene",
#                          values_to = "Expression")
# 
# # 3. 绘制箱线图 (基础版本)
# p <- ggplot(exp_long, aes(x = subtype, y = Expression, fill = subtype)) +
#   geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +  # 按亚型填充颜色
#   facet_wrap(~ Gene, ncol = 2, scales = "free_y") +  # 按基因分面
#   labs(x = "Cell Subtype", 
#        y = "Expression Level",
#        title = "Expression of KIT and TPSB2 across Cell Subtypes") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 旋转X轴标签
# 
# print(p)
# ggsave(p, file = "start.png", width = 6, height = 4)

# # 计算M1评分
# sce1 <- AddModuleScore(
#   object = sce1,
#   features = list(ciliated_genes),
#   name = "start_score"
# )
# # 1. 提取表达数据
# gene_exp_data <- FetchData(sce1, vars = c("subtype","start_score1"))
# # 3. 绘制箱线图 (基础版本)
# p <- ggplot(gene_exp_data, aes(x = subtype, y = start_score1, fill = subtype)) +
#   geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +  # 按亚型填充颜色
#   #facet_wrap(~ Gene, ncol = 2, scales = "free_y") +  # 按基因分面
#   # labs(x = "Cell Subtype", 
#   #      y = "Expression Level",
#   #      title = "Expression of KIT and TPSB2 across Cell Subtypes") +
#   theme_classic() +
#   scale_fill_manual(values = colors_set3_named) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 旋转X轴标签
# 
# print(p)
# ggsave(p, file = "start2.pdf", width = 6, height = 4)


plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds, 
            use_partition = FALSE,  #使用 use_partition = TRUE（默认）的情况：
            #当你的细胞群体包含多个完全独立的谱系或细胞类型，它们之间没有直接的发育转换关系时。
            #例如：你的数据同时含有肥大细胞、B细胞、T细胞。它们属于不同的造血分支，不应被连成一条轨迹。此时，让每个partition独立学习轨迹是正确的。
            learn_graph_control = list(
              minimal_branch_len = 5,  # 设置很大的值来抑制分支
              #ncenter = 100,            # 增加中心点使轨迹更平滑
              prune_graph = TRUE
              )
            )
plot_cells(cds,
           color_cells_by = "subtype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)  #Manual node selection

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 trajectory_graph_segment_size = 0.5,
                 # 控制轨迹线的颜色
                 trajectory_graph_color = "black",      # 改变线的颜色
                 alpha = 1)
p1
ggsave(p1, file = "time1.png",width = 6,height = 4)


p2 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 trajectory_graph_segment_size = 0,
                 rasterize  = T,
                 graph_label_size = 0,  #隐藏所有图标签
                 label_branch_points=FALSE) +
  scale_color_gradientn(colors = c("#4361EE", "#9D4EDD", "#FF70A6","red")) +
  labs(color = "Pseudotime")
ggsave(p2, file = "time2.png",width = 6,height = 4)
ggsave(p2, file = "time2.pdf",width = 6,height = 4)

# 寻找随拟时序变化的基因
deg_cds <- graph_test(cds, 
                      neighbor_graph = "principal_graph",
                      cores = 4)
# 查看最显著的基因
deg_cds <- deg_cds[order(deg_cds$morans_I, decreasing = TRUE), ]
head(deg_cds)
save(deg_cds,file = "deg.Rdata")
# 提取前6个基因进行可视化
top_genes <- row.names(deg_cds)[1:30]
p1 <- plot_genes_in_pseudotime(cds[top_genes,], 
                               color_cells_by = "pseudotime",
                               min_expr = 0.5,
                               ncol = 2)
ggsave(p1, file = "top6.png",width = 7, height = 6)

p1 <- plot_genes_in_pseudotime(cds[c(top_genes),], 
                               color_cells_by = "subtype",
                               min_expr = 0.5,
                               ncol = 2)
p1
ggsave(p1, file = "top6_2.png",width = 7, height = 20)
save(cds,file = "time.Rdata")

#load("time.Rdata")
pseudotime_df <- data.frame(
  cell_id = colnames(cds),
  pseudotime = pseudotime(cds),
  sample_id = colData(cds)$sample,  # 根据实际的列名调整
  row.names = colnames(cds)
)
pseudotime_df$group <- meta$group[match(pseudotime_df$sample_id,meta$id)]

library(ggplot2)

# 假设数据框结构：pseudotime_df 包含 sample_id, pseudotime, group 等列
ggplot(pseudotime_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Group",
       x = "Group",
       y = "Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

load("../tex_gsva_data.Rdata")
gsva <- gsva_data$gsva
identical(colnames(gsva),pseudotime_df$cell_id)

pseudotime_df$Progenitor_exhaustion <- gsva["Progenitor_exhaustion",]

ggplot(pseudotime_df, aes(x = pseudotime, y = Progenitor_exhaustion, color = group)) +
  geom_point() +
  theme_classic() +
  labs(
       y = "Progenitor_exhaustion",
       x = "Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
