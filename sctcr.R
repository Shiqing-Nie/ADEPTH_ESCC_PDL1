library(Seurat)
library(tidyverse)
library(scRepertoire)
dir1 <- basename(list.dirs("../TCR_Ranger/", recursive = FALSE))
rm1 <- c("XCH_pre","LGY_pre","WDF_pre")
dir1 <- dir1[!(dir1 %in% rm1)]
sce_t <- readRDS("../4-subtype/T_NK/sce.annote.rds")
load("../4-subtype/T_NK/finaltype/t_meta2.Rdata")
t_meta <- t_meta[t_meta$subtype != "NK",]
sce_t <- sce_t[,rownames(t_meta)]
a <- sum(rownames(sce_t@meta.data) != rownames(t_meta))
if (a > 0){
  print("something was wrong")
}
sce_t@meta.data <- t_meta
DimPlot(sce_t,group.by = "subtype",label = T)

z <- as.data.frame(table(sce_t$orig.ident))


meta <- read.delim("../meta.txt")

z$id <- meta$id1[match(z$Var1,meta$id)]

#dir1[!(dir1 %in% z$id)]
#z$id[!(z$id %in% dir1)]
#union1 <- intersect(dir1,z$id)
#dir1 <- dir1[dir1 %in% union1]

contig.list <- list()
for (i in dir1) {
  contig.list[[i]] <- read.csv(paste0("../TCR_Ranger/",i,"/outs/filtered_contig_annotations.csv"))
}

c <- c()
for (i in 1:length(contig.list)) {
  z <- contig.list[[i]]
  z$barcode <- paste0(names(contig.list)[i],"_",z$barcode)
  z1 <- unique(z$barcode)
  a <- sum((z1 %in% colnames(sce_t)))
  c <- c(c,a)
  print(a)
  z2 <- unique(z$high_confidence)
  if (length(z2) > 1){
    print("something was wrong")
  }
}
z <- data.frame(name = names(contig.list), n = c)
rm2 <- z$name[z$n < 50]
contig.list[rm2] <- NULL

dir1 <- names(contig.list)
id1 <- meta$id1[match(dir1,meta$id)]
sce_t <- sce_t[,sce_t$orig.ident %in% dir1]


meta <- meta[meta$id %in% unique(sce_t$orig.ident),]
meta <- unique(meta)
rownames(meta) <- meta$id
z <- names(table(meta$number))[table(meta$number) > 1]

meta2 <- meta[meta$number %in% z,]

table(meta2$effect)

a <- sum(dir1 !=  names(contig.list))
if (a > 0){
  print("something was wrong")
}

combined.TCR <- combineTCR(
  contig.list,
  samples = dir1,
  ID = NULL,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE,
  filterNonproductive = TRUE
)


#length(unique(combined.TCR[["LQY_post"]]$barcode))

sugery <- meta$surgery[match(names(combined.TCR),meta$id)]
id <- meta$id1[match(names(combined.TCR),meta$id)]
effect <- meta$effect[match(names(combined.TCR),meta$id)]
group <- paste0(sugery,"_",effect)
library(stringr)
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "surgery", 
                            variables = sugery)
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "id", 
                            variables = id)
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "effect", 
                            variables = effect)
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "group", 
                            variables = group)

for (i in 1:length(combined.TCR)) {
  id1 <- names(combined.TCR)[i]
  id2 <- meta$id1[match(id1,meta$id)]
  z <- combined.TCR[[id1]]
  names(combined.TCR)[i] <- id2
  z$sample <- id2
  combined.TCR[[id2]] <- z
}
combined.TCR <- combined.TCR[order(names(combined.TCR))]
##########伪bulk TCR##################################
dir.create("bulk")
setwd("bulk/")
p1 <- clonalQuant(combined.TCR, 
                  cloneCall="strict", 
                  chain = "both", 
                  scale = TRUE)
width1 <- length(dir1)
ggsave(p1,file = "unique_clone.png",height = 6,width = width1)
p1 <- clonalQuant(combined.TCR, cloneCall = "gene", group.by = "group", scale = TRUE)
ggsave(p1,file = "unique_clone_group.png",height = 6,width = 6)

p1 <- clonalProportion(combined.TCR, 
                       cloneCall = "gene") 
ggsave(p1,file = "clonalProportion.png",height = 6,width = width1)

p1 <- clonalDiversity(combined.TCR, 
                      cloneCall = "gene")
p1
tcr_bulk_diversity <- p1@data
save(tcr_bulk_diversity,file = "tcr_bulk_diversity.Rdata")
ggsave(p1,file = "diversity.png",height = 6,width = 10)

clonalOverlap(combined.TCR, 
              cloneCall = "gene", 
              method = "morisita") +
  theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                   hjust = 1,    # 右对齐
                                   vjust = 1))   # 向下微调
#z <- combined.TCR$P05_T_Pre$CTstrict[(combined.TCR$P05_T_Pre$CTstrict %in% combined.TCR$P05_T_Post$CTstrict)]
#z <- length(unique(z))
p1 <- clonalOverlap(combined.TCR, 
                    cloneCall = "strict", 
                    group.by = NULL,  #可以改成细胞？
                    order.by = NULL,
                    method = "morisita")
ggsave(p1,file = "clonalOverlap.png",height = 6,width = 10)


setwd("../")

###############bulk-choice##################################
for (i in 1:length(combined.TCR)) {
  tcr_data <- combined.TCR[[i]]
  # 找出在 sce 中存在的 barcode
  keep <- tcr_data$barcode %in% colnames(sce_t)
  combined.TCR[[i]] <- tcr_data[keep, ]
}
# 删除行数为 0 的元素
combined.TCR <- Filter(function(x) nrow(x) > 0, combined.TCR)

all_barcodes <- unlist(lapply(combined.TCR, function(x) x$barcode))

sce_t <- sce_t[,all_barcodes]
dir.create("bulk-choice")
setwd("bulk-choice/")

p1 <- clonalQuant(combined.TCR, 
                  cloneCall="strict", 
                  chain = "both", 
                  scale = TRUE)
width1 <- length(dir1)
ggsave(p1,file = "unique_clone.png",height = 6,width = width1)
p1 <- clonalQuant(combined.TCR, cloneCall = "gene", group.by = "group", scale = TRUE)
ggsave(p1,file = "unique_clone_group.png",height = 6,width = 6)

p1 <- clonalProportion(combined.TCR, 
                       cloneCall = "gene") 
ggsave(p1,file = "clonalProportion.png",height = 6,width = width1)

p1 <- clonalDiversity(combined.TCR, 
                      cloneCall = "gene")
ggsave(p1,file = "diversity.png",height = 6,width = 10)

clonalOverlap(combined.TCR, 
              cloneCall = "strict", 
              method = "raw")

p1 <- clonalOverlap(combined.TCR, 
                    cloneCall = "strict", 
                    group.by = NULL,  #可以改成细胞？
                    order.by = NULL,
                    method = "morisita") + 
  theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                   hjust = 1,    # 右对齐
                                   vjust = 1))   # 向下微调
ggsave(p1,file = "clonalOverlap_morisita.png",height = 6,width = 10)

p2 <- clonalOverlap(combined.TCR, 
                    cloneCall = "strict", 
                    group.by = NULL,  #可以改成细胞？
                    order.by = NULL,
                    method = "raw") + 
  theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                   hjust = 1,    # 右对齐
                                   vjust = 1))   # 向下微调
ggsave(p2,file = "clonalOverlap_raw.png",height = 6,width = 10)


# library(stringr)
# data_t <- p2$data %>%
#   #filter(str_detect(Var1, "_Pre$")) %>%
#   #filter(str_detect(Var2, "_Post$")) %>%
#   mutate(Var1 = as.character(Var1)) %>%
#   mutate(Var2 = as.character(Var2)) %>%
#   pivot_wider(names_from = Var2, values_from = value) %>%
#   as.data.frame() %>%
#   column_to_rownames(var = "Var1") %>%
#   select(ends_with("_Pre"))
# data_t <- data_t[str_detect(rownames(data_t), "_Post$"),]
# data_t <- data_t[,order(colnames(data_t))]
# data_t <- data_t[order(rownames(data_t)),]
# library(pheatmap)
# pheatmap(data_t,na_col = "white",
#          display_numbers = TRUE,
#          number_color    = "white",
#          cluster_rows = F,
#          cluster_cols = F)


combined.TCR2 <- combined.TCR
sce_t$seurat_clusters <- sce_t$subtype
for (i in 1:length(combined.TCR2)) {
  temp1 <- combined.TCR2[[i]]
  temp1$cluster <- sce_t@meta.data$seurat_clusters[match(temp1$barcode,rownames(sce_t@meta.data))]
  temp1$cluster <- as.character(temp1$cluster)
  #temp1$cluster <- paste0("C",temp1$cluster)
  combined.TCR2[[i]] <- temp1
}


# temp1$barcode[1]
# match(temp1$barcode[1],rownames(sce@meta.data))
# rownames(sce@meta.data)[9420]
methods_vec <- c("raw", "overlap", "jaccard", "cosine", "morisita")
for (k in methods_vec) {
  p1 <- clonalOverlap(combined.TCR2, 
                      cloneCall = "strict", 
                      group.by = NULL,  #可以改成细胞？
                      order.by = NULL,
                      method = k) + 
    theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                     hjust = 1,    # 右对齐
                                     vjust = 1))   # 向下微调
  ggsave(p1,file = paste0("clonalOverlap_",k,".png"),height = 6,width = 10)
  #“raw”给绝对计数；jaccard/overlap 看“存在”；cosine/morisita 把丰度也拉进来。
}
setwd("../")


###############single##################################
dir.create("single")
setwd("single")


load("../../4-subtype/subid.Rdata")

sce2_s <- combineExpression(combined.TCR2, 
                          sce_t, 
                          cloneCall="strict", #按照strict列计算clonsize,默认为strict
                          #肿瘤或感染研究：
                          #若关注精细克隆亚群，建议 cloneCall="nt" 或 "strict"。
                          group.by = "sample", #group.by = NULL
                          #默认不分组，即全局统计克隆频率（所有细胞一起计算）。此时,clonalFrequency 表示该克隆在 所有样本中的总细胞数。 
                          proportion = F,  #Whether to proportion (TRUE) or total frequency (FALSE) of the clone based on the group.by variable.
                          cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
z <- sce2_s@meta.data[sce2_s$cloneSize == "Hyperexpanded (100 < X <= 500)",]
z1 <- z[z$CTstrict == z$CTstrict[1],]
table(z1$SCT_snn_res.1) #多个亚群共享同一克隆.

z <- sce2_s@meta.data
sum(z$clonalProportion)

table(sce2_s@meta.data$cloneSize)

n1 <- 1:nrow(z)
n1 <- n1[z$clonalFrequency > 100]
a <- 1873
n3 <- z$clonalFrequency[a]
n3
z$orig.ident[a]
n2 <- sum((z$CTstrict == z$CTstrict[a]) & (z$orig.ident ==  z$orig.ident[a]))
n2

sce2_gene <- combineExpression(combined.TCR2, 
                            sce_t, 
                            cloneCall="gene", #按照strict列计算clonsize,默认为strict
                            #肿瘤或感染研究：
                            #若关注精细克隆亚群，建议 cloneCall="nt" 或 "strict"。
                            group.by = "sample", #group.by = NULL
                            #默认不分组，即全局统计克隆频率（所有细胞一起计算）。此时,clonalFrequency 表示该克隆在 所有样本中的总细胞数。 
                            proportion = F,  #Whether to proportion (TRUE) or total frequency (FALSE) of the clone based on the group.by variable.
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
z <- sce2_gene@meta.data
n3 <- z$clonalFrequency[a]
n3
z$orig.ident[a]
n2 <- sum((z$CTgene == z$CTgene[a]) & (z$orig.ident ==  z$orig.ident[a]))
n2
n2 == n3  #手动计算clonsize是否正确

sce3_agene <- combineExpression(combined.TCR2, 
                               sce_t, 
                               cloneCall="gene", #按照strict列计算clonsize,默认为strict
                               #肿瘤或感染研究：
                               #若关注精细克隆亚群，建议 cloneCall="nt" 或 "strict"。
                               #group.by = "sample", 
                               group.by = NULL,
                               #默认不分组，即全局统计克隆频率（所有细胞一起计算）。此时,clonalFrequency 表示该克隆在 所有样本中的总细胞数。 
                               proportion = F,  #Whether to proportion (TRUE) or total frequency (FALSE) of the clone based on the group.by variable.
                               cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
z <- sce3_agene@meta.data
n3 <- z$clonalFrequency[a]
n3
z$orig.ident[a]
n2 <- sum((z$CTgene == z$CTgene[a]) & z$orig.ident == z$orig.ident[a])
n2
n2 == n3  #手动计算clonsize是否正确  #好像还是按样品分组了。

sce3_astrict <- combineExpression(combined.TCR2, 
                                sce_t, 
                                cloneCall="strict", #按照strict列计算clonsize,默认为strict
                                #肿瘤或感染研究：
                                #若关注精细克隆亚群，建议 cloneCall="nt" 或 "strict"。
                                #group.by = "sample", #group.by = NULL
                                #默认不分组，即全局统计克隆频率（所有细胞一起计算）。此时,clonalFrequency 表示该克隆在 所有样本中的总细胞数。 
                                proportion = F,  #Whether to proportion (TRUE) or total frequency (FALSE) of the clone based on the group.by variable.
                                cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
z <- sce3_astrict@meta.data
n3 <- z$clonalFrequency[a]
n3
z$orig.ident[a]
n2 <- sum((z$CTgene == z$CTgene[a]))
n2
n2 == n3  #手动计算clonsize是否正确


z_gene <- as.data.frame(table(sce3_agene@meta.data$SCT_snn_res.1,sce3_agene@meta.data$cloneSize))


data1 <- list(strict = sce3_astrict, gene = sce3_agene)
data1_g <- list(strict = sce2_s, gene = sce2_gene)
type2 <- c("subtype","celltype")

dim(sce_t)
dim(sce3_agene)
sce_tc <- sce3_agene[,sce3_agene$subtype == "CD8 Tc"]
tc_meta <- sce_tc@meta.data
tc_meta <- tc_meta %>% select(1,25,27:33)
sce.all.int <- readRDS("../../4-subtype/T_NK/finaltype/cd8_tc/sce.all_int.rds")
#unique(sce.all.int$SCT_snn_res.0.2)
#DimPlot(sce.all.int, group.by = "SCT_snn_res.0.2")
#tc_meta$subtype2 <- sce.all.int$RNA_snn_res.0.2[match(rownames(tc_meta),rownames(sce.all.int@meta.data))]
#tc_meta$group <- meta$group1[match(tc_meta$orig.ident,meta$id)]
# 统计每个样本的唯一序列数
#result <- tc_meta %>%
#  group_by(orig.ident) %>%
#  summarise(unique_sequences = n_distinct(CTnt))
#result$group <- meta$group1[match(result$orig.ident,meta$id)]
#comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
#comp1 <- list(c("Pre_R","Post_R"),c("Pre_NR","Post_NR"))
#comp2 <- list(c("Post_R","Post_NR"))
#order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
#color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
#result$group <- factor(result$group,levels = order1)
#library(scales)
#library(ggpubr)
#p1 <- ggplot(result, aes(x = group, y = unique_sequences)) +
#  geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
#  scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
#  #geom_point(size = 3,color='red') +
#  #geom_line(aes(group = id), color = 'gray', lwd = 0.5) +  #配对样本间连线
#  geom_point() +  #绘制散点
#  scale_y_continuous(
#    labels = label_number(accuracy = 1)  # accuracy=1 表示整数
#  ) +
#  theme(panel.grid = element_blank(),
#        axis.line = element_line(colour = 'black', linewidth = 1),
#        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
#        panel.background = element_blank(),
#        plot.title = element_text(size = 15, hjust = 0.5),
#        plot.subtitle = element_text(size = 15, hjust = 0.5),
#        axis.text = element_text(size = 15, color = 'black'),
#        axis.title = element_text(size = 15, color = 'black')) +
#  labs(x = '', y = 'Richness')+
#  stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
#p1
#ggsave(p1,file = "Tc_richness.png",width = 6, height = 5)

#result <- tc_meta %>% filter(subtype2 == 1) %>%
#  group_by(orig.ident) %>%
#   summarise(unique_sequences = n_distinct(CTnt))
# result$group <- meta$group1[match(result$orig.ident,meta$id)]
# result$group <- factor(result$group,levels = order1)
# library(scales)
# library(ggpubr)
# p1 <- ggplot(result, aes(x = group, y = unique_sequences)) +
#   geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
#   scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
#   #geom_point(size = 3,color='red') +
#   #geom_line(aes(group = id), color = 'gray', lwd = 0.5) +  #配对样本间连线
#   geom_point() +  #绘制散点
#   scale_y_continuous(
#     labels = label_number(accuracy = 1)  # accuracy=1 表示整数
#   ) +
#   theme(panel.grid = element_blank(),
#         axis.line = element_line(colour = 'black', linewidth = 1),
#         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
#         panel.background = element_blank(),
#         plot.title = element_text(size = 15, hjust = 0.5),
#         plot.subtitle = element_text(size = 15, hjust = 0.5),
#         axis.text = element_text(size = 15, color = 'black'),
#         axis.title = element_text(size = 15, color = 'black')) +
#   labs(x = '', y = 'Richness')+
#   stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
# p1
# ggsave(p1,file = "Tc_c1_richness.png",width = 6, height = 5)

# z <- contig.list
# for (i in names(z)) {
#    z1 <- z[[i]]
#    z1$barcode <- paste0(i,"_",z1$barcode)
#    z[[i]] <- z1
# }
# z <- bind_rows(z)
# z <- z[z$barcode %in% rownames(tc_meta),]
# tc_meta <- tc_meta %>% dplyr::select(orig.ident,subtype,subtype2,group)
# tc_meta$barcode <- rownames(tc_meta)
# tc_meta2 <- full_join(z, tc_meta, by = c("barcode" = "barcode"))
# 
# save(tc_meta2,file = "tc_meta.Rdata")

library(ggpubr)
i <- names(data1)[1]
for (i in names(data1)){
  unlink(i,recursive = T)
  dir.create(i)
  setwd(i)
  data2 <- data1[[i]]
  sce2 <- data2
  j <- type2[1]
  for (j in type2) {
    dir.create(j)
    setwd(j)
#    if (j == "subtype"){
#      sce2$cluster <- sce2@meta.data[,j]
    # }else{
    #   sce2$cluster <- paste0("C",sce2@meta.data[,j])
    # }
    sce2$cluster <- sce2@meta.data[,j]
    color1 <- c(
      "Hyperexpanded (100 < X <= 500)" = "#FF8888",  # 浅红色，更柔和
      "Large (20 < X <= 100)" = "#FFD166",           # 
      "Medium (5 < X <= 20)" = "#A3A500",            # 浅金黄色
      "Small (1 < X <= 5)" = "#66BB99",              # 浅青绿色
      "Single (0 < X <= 1)" = "#6699CC"              # 浅蓝色
    )
    
    p1 <- DimPlot(data2,group.by = "cloneSize",raster = T,pt.size = 2) +
      labs(title = "") + theme_void() +
      scale_color_manual(values = color1)
    p1  
    ggsave(p1, file = "whole_umap.pdf", width = 8, height = 5)
    ggsave(p1, file = "whole_umap.png", width = 8, height = 5)
    
    pa1 <-  as.data.frame(table(sce2$orig.ident,sce2$cloneSize))
    pa2 <- pa1 %>% 
      group_by(Var1) %>%
      mutate(percent = Freq / sum(Freq) * 100) %>%
      ungroup() %>% filter(Var2 != "None ( < X <= 0)")
    ggplot(pa2, aes(x = Var1, y = percent)) +
      geom_bar(aes(fill = Var2) , stat = "identity") +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
            # 旋转x轴标签45度
            text = element_text(color = "black"),          # 所有文本
            title = element_text(color = "black"),         # 所有标题
            legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = " ", fill = NULL)+labs(y = 'Relative proportion (%)') +
      scale_fill_manual(values = color1) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))
    saveRDS(pa2,file = "patient_expand.rds")
    
    data3 <-  as.data.frame(table(data2@meta.data[,j],data2@meta.data$cloneSize,data2@meta.data$orig.ident))
    colnames(data3) <- c("cluster","class","id","freq")
    data3 <- data3[data3$class != "None ( < X <= 0)",]
    # if (j != "subtype"){
    #   data3$cluster <- paste0("C",data3$cluster)
    # }
    data3_2 <- data3
    data3 <- data3 %>%  select(-id) %>%
      group_by(cluster) %>%
      mutate(percent = freq / sum(freq) * 100) %>%
      ungroup()
    p1 <- ggplot(data3, aes(x = cluster, y = percent)) +
      geom_bar(aes(fill = class) , stat = "identity") +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
            # 旋转x轴标签45度
            text = element_text(color = "black"),          # 所有文本
            title = element_text(color = "black"),         # 所有标题
            legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = " ", fill = NULL)+labs(y = 'Relative proportion (%)') +
      scale_fill_manual(values = color1) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
    
    p1

    l1 <- length(unique(data3$cluster))
    ggsave(p1, file = "whole_cluster.pdf",width = l1 / 2 + 3, height = 5)
    ggsave(p1, file = "whole_cluster.png",width = l1 / 2 + 3, height = 5)
    z <- p1@data
    library(writexl)
    write_xlsx(z, "whole_cluster.xlsx")
    
    data3_2$group <- meta$group[match(data3_2$id,meta$id)]
    p_data <- list()
    group3 <- unique(meta$group)
    
    sce2$group <- meta$group[match(sce2$orig.ident,meta$id)]
    for (g1 in group3) {
      #data3_3 <- data3_2[data3_2$group == g1,]
      data3_3 <- data3_2
      #data3_3 <- data3_2[data3_2$cluster %in% c("CD4 Treg","CD8 Tc","CD8 Tpex","CD8 Tex-term"),]
      data3_4 <- data3_3 %>% group_by(cluster,group,class) %>%
        summarise(
          freq = sum(freq, na.rm = TRUE),  # 求和，忽略NA值
          .groups = 'drop'  # 取消分组，防止后续操作出错
        )
      
      data3_4 <- data3_4 %>%
        group_by(cluster,group) %>%
        mutate(percent = freq / sum(freq) * 100) %>%
        ungroup()
      
      data33_4 <- data3_4 #合并治疗前后
      data33_4$group[data33_4$group %in% c("Pre_R","Post_R")] <- "R"
      data33_4$group[data33_4$group %in% c("Pre_NR","Post_NR")] <- "NR"
      data33_4 <- data33_4 %>%
        group_by(cluster,group) %>%
        mutate(percent = freq / sum(freq) * 100) %>%
        ungroup()
      p1 <- ggplot(data33_4, aes(x = cluster, y = percent)) +
        geom_bar(aes(fill = class) , stat = "identity") +
        facet_wrap(~group, scales = "free_x") +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
              # 旋转x轴标签45度
              text = element_text(color = "black"),          # 所有文本
              title = element_text(color = "black"),         # 所有标题
              legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "",y = "Relative abundance (%)", fill = NULL) +
        scale_fill_manual(values = color1) +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间

      p1
      ggsave(p1, file = paste0("whole_cluster_two.pdf"),width = l1 / 1.5 + 3, height = 5)
      ggsave(p1, file = paste0("whole_cluster_two.png"),width = l1 / 1.5 + 3, height = 5)
      
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      data3_4$group <- factor(data3_4$group, levels = order1)
      p1 <- ggplot(data3_4, aes(x = class, y = cluster)) +
        geom_tile(aes(fill = percent), color = "white", linewidth = 0.5) +
        # 添加数值标签
        geom_text(
          aes(label = sprintf("%.1f", percent)),  # 显示1位小数
          color = "white",  # 标签颜色
          size = 3,         # 字体大小
          fontface = "bold"
        ) +
        facet_wrap(~group) +
        scale_fill_viridis_c(
          name = "Relative abundance (%)",
          option = "viridis",
          direction = -1,
          limits = c(0, 100),
          na.value = "grey90"
        ) +
        theme_void() +
        theme(
          panel.grid = element_blank(),
          text = element_text(color = "black"),
          title = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "right",
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(0.5, "cm")
        ) +
        labs(x = "", y = "", fill = "Relative abundance (%)")
      p1
      
      
      l1 <- length(unique(data3_4$cluster))
      ggsave(p1, file = paste0("whole_cluster_four.pdf"),width = l1 / 1.5 + 3, height = 5)
      ggsave(p1, file = paste0("whole_cluster_four.png"),width = l1 / 1.5 + 3, height = 5)
      z <- p1@data
      #write_xlsx(z, "whole_cluster_four.xlsx")
      sce2_t <- sce2[,sce2$group == g1]
      sce2_t <- sce2_t[,sce2_t$cluster %in% c("CD4 Treg","CD8 Tc","CD8 Tpex","CD8 Tex-term")]
      #sce2_t$cloneSize <- factor(sce2_t$cloneSize, levels = level1)
      # color1 <- c(
      #   "Hyperexpanded (100 < X <= 864)" = "#F8766D",
      #   "Large (20 < X <= 100)" = "#00BA38", 
      #   "Medium (5 < X <= 20)" = "#619CFF",
      #   "Small (1 < X <= 5)" = "#F564E3",
      #   "Single (0 < X <= 1)" = "#00BFC4"
      # )
      color1 <- c(
        "Hyperexpanded (100 < X <= 500)" = "#FF8888",  # 浅红色，更柔和
        "Large (20 < X <= 100)" = "#FFD166",           # 
        "Medium (5 < X <= 20)" = "#A3A500",            # 浅金黄色
        "Small (1 < X <= 5)" = "#66BB99",              # 浅青绿色
        "Single (0 < X <= 1)" = "#6699CC"              # 浅蓝色
      )
      
      
      p1 <- DimPlot(sce2_t,group.by = "cloneSize") + labs(title = "") +
        scale_color_manual(
          values = color1,
          drop = FALSE
        )
      p1
      ggsave(p1, file = paste0("whole_umap_",g1,".pdf"), width = 8, height = 5)
      ggsave(p1, file = paste0("whole_umap_",g1,".png"), width = 8, height = 5)
    }
    
    # meta3 <-  sce2@meta.data %>% select(subtype,cluster)
    # meta3 <- unique(meta3)
    # data3$type <- meta3$subtype[match(data3$cluster,meta3$cluster)]
    # data3 <- data3 %>% select(class,freq,type)
    # 
    # data4 <- data3 %>%
    #   group_by(class, type) %>%
    #   summarise(value = sum(freq, na.rm = TRUE),  # 累加value列
    #             .groups = "drop")  # 自动解除分组
    # data4 <- data4 %>%
    #   group_by(type) %>%
    #   mutate(percent = value / sum(value) * 100) %>%
    #   ungroup()
    # p1 <- ggplot(data4, aes(x = type, y = percent)) +
    #   geom_bar(aes(fill = class) , stat = "identity") +
    #   theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
    #         # 旋转x轴标签45度
    #         legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
    #         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    #   labs(y = " ", fill = NULL)+labs(y = 'Relative proportion (%)') +
    #   scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
    # 
    # p1
    # l1 <- length(unique(data4$type))
    # ggsave(p1, file = "whole_type.pdf",width = (l1 / 2) + 4, height = 6)
    
    library(ggraph)
    #No Identity filter
    # p1 <- clonalNetwork(sce2_s,
    #                     reduction = "umap",
    #                     group.by = "subtype",
    #                     filter.clones = NULL,
    #                     filter.identity = NULL,
    #                     cloneCall = "aa")
    # ggsave(p1,file = "cloneshare1.png",height = 6,width = 8)
    #Examining Cluster 4 only
    # clonalNetwork(sce2, 
    #               reduction = "umap", 
    #               group.by = "seurat_clusters",
    #               filter.identity = 0,
    #               cloneCall = "aa")
    # shared.clones <- clonalNetwork(sce2, 
    #                                reduction = "umap", 
    #                                group.by = "cluster",
    #                                cloneCall = "aa", 
    #                                exportClones = TRUE)
    
    Idents(sce2) <- j
    methods_vec <- c("raw", "overlap", "jaccard", "cosine", "morisita")
    dir.create("overlap")
    for (k in methods_vec) {
      p1 <- clonalOverlap(sce2, 
                          cloneCall = "strict", 
                          group.by = NULL,  #可以改成细胞？
                          order.by = NULL,
                          method = k) + 
        theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                         hjust = 1,    # 右对齐
                                         vjust = 1))   # 向下微调
      ggsave(p1,file = paste0("overlap/all_clonalOverlap_",k,".pdf"),height = 6,width = 10)
      #“raw”给绝对计数；jaccard/overlap 看“存在”；cosine/morisita 把丰度也拉进来。
    }
    
    for (k in methods_vec) {
      p1 <- clonalOverlap(sce2[,sce2$group %in% c("Post_R","Pre_R")], 
                          cloneCall = "strict", 
                          group.by = NULL,  #可以改成细胞？
                          order.by = NULL,
                          method = k) + 
        theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                         hjust = 1,    # 右对齐
                                         vjust = 1))   # 向下微调
      ggsave(p1,file = paste0("overlap/R_clonalOverlap_",k,".pdf"),height = 6,width = 10)
      #“raw”给绝对计数；jaccard/overlap 看“存在”；cosine/morisita 把丰度也拉进来。
    }
    
    for (k in methods_vec) {
      p1 <- clonalOverlap(sce2[,sce2$group %in% c("Post_NR","Pre_NR")], 
                          cloneCall = "strict", 
                          group.by = NULL,  #可以改成细胞？
                          order.by = NULL,
                          method = k) + 
        theme(axis.text.x = element_text(angle = 45,    # 倾斜角度
                                         hjust = 1,    # 右对齐
                                         vjust = 1))   # 向下微调
      ggsave(p1,file = paste0("overlap/NR_clonalOverlap_",k,".pdf"),height = 6,width = 10)
      #“raw”给绝对计数；jaccard/overlap 看“存在”；cosine/morisita 把丰度也拉进来。
    }
    
    id2 <- unique(sce2@meta.data$orig.ident)
    library(stringr)
    id2 <- str_split_fixed(id2, "_", 2)[, 1]
    z <- table(id2)
    id2 <- names(z)[z == 2]
    sce2$shared <- "NA"
    data_share <- list()
    for (k in id2) {
      pa1 <- paste0(k,"_","pre")
      pa2 <- paste0(k,"_","post")
      sce2$shared[(sce2$orig.ident == pa1) & (sce2$CTnt %in% sce2$CTnt[sce2$orig.ident == pa2])] <- "Shared"
      z <- as.data.frame(table(sce2$orig.ident,sce2$shared))
      sce2$shared[(sce2$orig.ident == pa2) & (sce2$CTnt %in% sce2$CTnt[sce2$orig.ident == pa1])] <- "Shared"
      z <- as.data.frame(table(sce2$orig.ident,sce2$shared))
      
      ####判断共享####################
      z <- sce2@meta.data
      z1 <- z[z$orig.ident == pa1,]
      z2 <- z[z$orig.ident == pa2,]
      data_share2 <- data.frame(source = "a", target = "b", CTstrict = "c",
                                number = "d",effect = "e")
      for (k1 in 1:nrow(z1)) {
        id3 <- z1$CTstrict[k1]
        z3 <- z2[z2$CTstrict == id3,]
        if (nrow(z3) > 0){
          for (k2 in 1:nrow(z3)) {
            a <- z1$subtype[k1]
            b <-  z3$subtype[k2]
            c <- id3
            d <- k
            e <- meta$effect[match(pa1,meta$id)]
            data_share2[nrow(data_share2) + 1,] <- c(a,b,c,d,e)
          }
        }
      }
      data_share2 <- data_share2[-1,]
      data_share[[k]] <- data_share2
    }
    length(unique(sce2$CTnt[(sce2$orig.ident == pa1) & (sce2$shared == "Shared") ]))
    #17个，和bulk-choice中的raw能完全对应上。
    sce2$shared[sce2$shared == "NA"] <- "Not Shared"
    sce2$shared2 <- sce2$shared
    sce2$surgery <- str_extract(sce2$orig.ident, "[^_]+$")
    #sce2$shared2[(sce2$shared2 == "Shared") & (sce2$surgery == "pre")] <- "Shared_Pre"
    #sce2$shared2[(sce2$shared2 == "Shared") & (sce2$surgery == "post")] <- "Shared_Post"
    regex <- paste0("^(", paste(id2, collapse = "|"), ")")
    sce2_shared <- sce2[,grepl(regex,sce2$orig.ident)]
    
    col_s <- c(Shared = "#FF7043",  # 深橙色（温暖醒目）
               'Not Shared' =   "gray"   # 明亮蓝色
    )
    p1 <- DimPlot(sce2_shared,group.by = "shared2",alpha = 1) +
       scale_color_manual(
        values = col_s
      )  +
      labs(title = "")
    p1
    #"P04" "P05" "P06" "P08" "P18" "P20" "P22"只包括这7对
    ggsave(p1, file = "shared_umap.pdf", width = 6.5, height = 5)
    ggsave(p1, file = "shared_umap.png", width = 6.5, height = 5)
    
    pa1 <-  as.data.frame(table(sce2_shared$orig.ident,sce2_shared$shared2))
    pa2 <- pa1 %>% 
      group_by(Var1) %>%
      mutate(percent = Freq / sum(Freq) * 100)
    ggplot(pa2, aes(x = Var1, y = percent)) +
      geom_bar(aes(fill = Var2) , stat = "identity") +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
            # 旋转x轴标签45度
            text = element_text(color = "black"),          # 所有文本
            title = element_text(color = "black"),         # 所有标题
            legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = " ", fill = NULL)+labs(y = 'Relative proportion (%)') +
      scale_fill_manual(values = col_s) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))
    saveRDS(pa2,file = "patient_shared.rds")
    #############share流动############################
    data_share3 <- bind_rows(data_share)
    # 绘制桑基图
    table(data_share3$number,data_share3$effect)
    library(ggsankey)
    # 删除任何列中包含"CD4"开头的行
    data_share3 <- data_share3 %>%
      filter(!if_any(everything(), ~str_detect(., "^CD4"))) %>%
      dplyr::rename(Pre = source, Post = target)
    plot_data <- data_share3 %>%
      make_long(Pre, Post)
    colors_set3_named <- c(
      `CD4 Tcm-KLF2` = "#8DD3C7",      # 青绿色
      `CD4 Tcm-TALAM1` = "#FFFFB3",    # 浅黄色（新增）
      `CD4 Tem` = "#BEBADA",          # 浅紫色（新增）
      `CD4 Tfh` = "#E58606",          # 橙红色
      `CD4 Treg` = "#FDB462",         # 橙色
      `CD8 Proliferation` = "#80B1D3", # 天蓝色
      `CD8 Tc` = "#ECCBAE",           # 米色
      `CD8 Tcm-ANXA1` = "#B3DE69",    # 黄绿色
      `CD8 Tcm-NFKBIA` = "#FCCDE5",   # 粉色
      `CD8 Tem-HSPA1A` = "#D9D9D9",   # 浅灰色
      `CD8 Tem-IFI6` = "#BC80BD",     # 紫红色（调整位置）
      `CD8 Temra` = "#FFED6F",        # 金黄色（新增）
      `CD8 Tex-term` = "#CCEBC5",     # 淡绿色
      `CD8 Tpex` = "#FB8072",         # 珊瑚红
      `NK` = "#8DA0CB"                # 淡蓝色
    )
    
    p1 <- ggplot(plot_data, aes(x = x, 
                                next_x = next_x, 
                                node = node, 
                                next_node = next_node,
                                fill = node,
                                label = node)) +
      scale_fill_manual(values = colors_set3_named) +
      geom_sankey(flow.alpha = 0.6, node.color = "gray", smooth = 6) +
      geom_sankey_label(size = 3.2, color = "black") +
      #scale_fill_viridis_d(option = "plasma") +
      theme_sankey(base_size = 16) +
      labs(x = NULL, 
           title = "") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(p1, file = paste0("share_flow.pdf"),width = 6,height = 10)
    ggsave(p1, file = paste0("share_flow.png"),width = 6,height = 10) 
    
    plot_data_tpex <- data_share3 %>% filter(Pre == "CD8 Tpex", Post %in% c("CD8 Tpex","CD8 Tex-term")) %>%
      make_long(Pre, Post)
    p1 <- ggplot(plot_data_tpex, aes(x = x, 
                                next_x = next_x, 
                                node = node, 
                                next_node = next_node,
                                fill = node,
                                label = node)) +
      scale_fill_manual(values = colors_set3_named) +
      geom_sankey(flow.alpha = 0.6, node.color = "gray", smooth = 6) +
      geom_sankey_label(size = 3.2, color = "black") +
      #scale_fill_viridis_d(option = "plasma") +
      theme_sankey(base_size = 16) +
      labs(x = NULL, 
           title = "") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(p1, file = paste0("share_flow_tpex.pdf"),width = 6,height = 5)
    ggsave(p1, file = paste0("share_flow_tpex.png"),width = 6,height = 5) 
    
    
    group1  <- c("R","NR")
    for (k3 in group1) {
      data_share4 <- data_share3[data_share3$effect == k3,]
      plot_data <- data_share4 %>%
        make_long(Pre, Post)
      p1 <- ggplot(plot_data, aes(x = x, 
                                  next_x = next_x, 
                                  node = node, 
                                  next_node = next_node,
                                  fill = node,
                                  label = node)) +
        geom_sankey(flow.alpha = 0.6, node.color = "black", smooth = 6) +
        geom_sankey_label(size = 3.2, color = "black") +
        #scale_fill_viridis_d(option = "plasma") +
        scale_fill_manual(values = colors_set3_named) +
        theme_sankey(base_size = 16) +
        labs(x = NULL, 
             title = "") +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      ggsave(p1, file = paste0("share_flow_",k3,".pdf"),width = 6,height = 10)
      ggsave(p1, file = paste0("share_flow_",k3,".png"),width = 6,height = 10) 
    }
    data_share3 <- data_share3 %>%
      filter(!if_any(everything(), ~str_detect(., "CD8 Tem|CD8 Tcm-|CD8 Proliferation")))
    
    class1 <- list(all = c("R","NR"),
                   R = "R",
                   NR = "NR")
    for (g1 in names(class1)) {
      data_share3_t <- data_share3[data_share3$effect %in% class1[[g1]],]
      mat <- data_share3_t %>%
        dplyr::count(Pre, Post) %>%
        pivot_wider(
          names_from = Pre,
          values_from = n,
          values_fill = 0  # 将缺失值填充为0
        ) %>% as.data.frame()
      rownames(mat) <- mat$Post
      mat <- mat %>% select(-Post)
      library(circlize)
      # 绘制和弦图
      pdf(paste("share_flow_circle_",g1,".pdf"), 
          width = 3, 
          height = 5
          #family = "Helvetica",      # 字体
          #pointsize = 12,            # 字号
          #paper = "a4",              # 纸张大小
          #bg = "white"
      )              # 背景色
      chordDiagram(as.matrix(mat), 
                   grid.col = colors_set3_named,
                   transparency = 0.3,
                   directional = 1,
                   #direction.type = "arrows",
                   link.arr.type = "big.arrow"
                   )
      dev.off()
      
      png(paste("share_flow_circle_",g1,".png"), 
          width = 6, 
          height = 10, 
          units = "in", 
          res = 300)  # 设置分辨率
      chordDiagram(as.matrix(mat), 
                   transparency = 0.3,
                   directional = 1,
                   grid.col = colors_set3_named,
                   #direction.type = "arrows",
                   link.arr.type = "big.arrow")
      dev.off()
    }
    
    ########整体shared############################
    sce2_shared$group <- meta$group[match(sce2_shared$orig.ident,meta$id)]
    sce2_shared2 <- sce2_shared[,sce2_shared$group %in% c("Post_R","Post_NR")]
    z <- as.data.frame(table(sce2_shared2$shared,sce2_shared2$cluster))
    z1 <- z %>%
      group_by(Var2) %>%
      mutate(percent = Freq / sum(Freq) * 100) %>%
      ungroup()
    col_s <- c(Shared = "#FF7043",  # 深橙色（温暖醒目）
               'Not Shared' =   "#42A5F5"   # 明亮蓝色
    )
    p1 <- ggplot(z1, aes(x = Var2, y = percent)) +
      geom_bar(aes(fill = Var1) , stat = "identity") +
      scale_fill_manual(values = col_s) +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 0.25),
            # 旋转x轴标签45度
            text = element_text(color = "black"),          # 所有文本
            title = element_text(color = "black"),         # 所有标题
            axis.text = element_text(color = "black"),     # 坐标轴刻度标签
            legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(y = 'Relative proportion (%)',x = "", fill = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
    l1 <- length(unique(z1$Var2))
    p1
    ggsave(p1, file = "whole_shared.pdf",width = l1 / 2 + 3, height = 5)
    ggsave(p1, file = "whole_shared.png",width = l1 / 2 + 3, height = 5)
    z_whole <- p1@data

    z <- as.data.frame(table(sce2_shared$shared,sce2_shared$cluster,sce2_shared$cloneSize))
    z <- z[z$Var2 %in% c("CD8 Tc","CD8 Tpex","CD8 Tex-term","CD8 Proliferation"),]
    z1 <- z %>%
      group_by(Var1,Var2) %>%
      mutate(percent = Freq / sum(Freq) * 100) %>%
      ungroup()
    p1 <- ggplot(z1, aes(x = Var1, y = percent)) +
      geom_bar(aes(fill = Var3) , stat = "identity") +
      scale_fill_manual(values = color1) +
      facet_wrap(~ Var2,nrow = 1,scales = "free") +
      
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 0.25),
            # 旋转x轴标签45度
            text = element_text(color = "black"),          # 所有文本
            title = element_text(color = "black"),         # 所有标题
            axis.text = element_text(color = "black"),     # 坐标轴刻度标签
            axis.title = element_text(color = "black"),     # 坐标轴标题
            legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(x = "", fill = NULL)+labs(y = 'Relative proportion (%)') +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
    p1
    l1 <- length(unique(z1$Var2))
    ggsave(p1, file = paste0("whole_shared2.pdf"),width = l1 / 2 + 6, height = 4.5)
    ggsave(p1, file = paste0("whole_shared2.png"),width = l1 / 2 + 6, height = 4.5)
    
    ###########只看治疗后共享#######################################
    group2 <- unique(meta$group)[c(3,4)]
    data4_list <- list()
    for (g1 in group2) {
      sce3_shared <- sce2_shared[,sce2_shared$group == g1]
      p1 <- DimPlot(sce3_shared,group.by = "shared") + 
        scale_color_manual(values = col_s) +
        labs(title = "")
      p1
      ggsave(p1, file = paste0("shared_umap_",g1,".pdf"), width = 6.5, height = 5)
      ggsave(p1, file = paste0("shared_umap_",g1,".png"), width = 6.5, height = 5)
      
      ########整体shared############################
      z <- as.data.frame(table(sce3_shared$shared,sce3_shared$cluster))
      z1 <- z %>%
        group_by(Var2) %>%
        mutate(percent = Freq / sum(Freq) * 100) %>%
        ungroup()
      p1 <- ggplot(z1, aes(x = Var2, y = percent)) +
        geom_bar(aes(fill = Var1) , stat = "identity") +
        scale_fill_manual(values = col_s) +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
              # 旋转x轴标签45度
              legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "", fill = NULL)+labs(y = 'Relative proportion (%)') +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
      l1 <- length(unique(z1$Var2))
      p1
      ggsave(p1, file = paste0("whole_shared_",g1,".pdf"),width = l1 / 2 + 3, height = 5)
      ggsave(p1, file = paste0("whole_shared_",g1,".png"),width = l1 / 2 + 3, height = 5)
      
      data4_list[[g1]] <- p1@data
      
      z <- as.data.frame(table(sce3_shared$shared,sce3_shared$cluster,sce3_shared$cloneSize))
      z <- z[z$Var2 %in% c("CD8 Tc","CD8 Tpex","CD8 Tex-term","CD4 Treg"),]
      z1 <- z %>%
        group_by(Var1,Var2) %>%
        mutate(percent = Freq / sum(Freq) * 100) %>%
        ungroup()
      p1 <- ggplot(z1, aes(x = Var1, y = percent)) +
        geom_bar(aes(fill = Var3) , stat = "identity") +
        scale_fill_manual(values = color1) +
        facet_wrap(~ Var2,nrow = 1,scales = "free") +
        
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 0.25),
              # 旋转x轴标签45度
              text = element_text(color = "black"),          # 所有文本
              title = element_text(color = "black"),         # 所有标题
              axis.text = element_text(color = "black"),     # 坐标轴刻度标签
              axis.title = element_text(color = "black"),     # 坐标轴标题
              legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(x = "", fill = NULL)+labs(y = 'Relative proportion (%)') +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
      p1
      l1 <- length(unique(z1$Var2))
      ggsave(p1, file = paste0("whole_shared_",g1,"2.pdf"),width = l1 / 2 + 6, height = 4.5)
      ggsave(p1, file = paste0("whole_shared_",g1,"2.png"),width = l1 / 2 + 6, height = 4.5)
    }
    data4_list$all <- z_whole
    write_xlsx(data4_list, "whole_shared.xlsx")
    
    ################分成四组########################
    z <- sce2_shared@meta.data
    #z <- z[z$surgery == "Post",]
    unique(z$shared2)
    z <- as.data.frame(table(z$shared,z$cluster,z$orig.ident))
    z1 <- z %>%
      group_by(Var2,Var3) %>%
      mutate(percent = Freq / sum(Freq) * 100) %>%
      ungroup()
    
    comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
    #comp1 <- list(c("Pre_R","Post_R"),c("Pre_NR","Post_NR"))
    comp2 <- list(c("Post_R","Post_NR"))
    order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
    color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")

    z1 <- z1[z1$Var1 == "Shared", ]
    z1 <- z1 %>%
      mutate(id = str_extract(Var3, "^[^_]+"))
    z1$group <- meta$group[match(z1$Var3,meta$id)]
    z1$group <- factor(z1$group,levels = order1)
    library(scales)
    p1 <- ggplot(z1, aes(x = group, y = percent)) +
      geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
      scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
      #geom_point(size = 3,color='red') +
      geom_line(aes(group = id), color = 'gray', lwd = 0.5) +  #配对样本间连线
      geom_point() +  #绘制散点
      facet_wrap( ~ Var2,scales = "free_y") +
      scale_y_continuous(
        labels = label_number(accuracy = 1)  # accuracy=1 表示整数
      ) +
      theme(panel.grid = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            panel.background = element_blank(),
            plot.title = element_text(size = 15, hjust = 0.5),
            plot.subtitle = element_text(size = 15, hjust = 0.5),
            axis.text = element_text(size = 15, color = 'black'),
            axis.title = element_text(size = 15, color = 'black')) +
      labs(x = '', y = 'Relative abundance (%)')+
      stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
    p1
    n1 <- length(unique(z1$Var2))
    n1 <- sqrt(n1)
    ggsave(p1, filename = paste0("four_effect_shared.png"),width = n1 * 2 + 4,height = n1 * 1.5)
    ggsave(p1, filename = paste0("four_effect_shared.pdf"),width = n1 * 2 + 4,height = n1 * 1.5)
    
    library(ggforce)
    p1 <- StartracDiversity(sce2, 
                            type = j,
                            #type = "RNA_snn_res.1", 
                            group.by = "orig.ident") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      scale_fill_manual(values = colors_set3_named)
    p1
    data_three <- p1@data
    
    p1 <- ggplot(data_three, aes(x = majorCluster, y = value)) +
      geom_boxplot(aes(fill = majorCluster), outliers = F) +
      facet_wrap(~ variable, scales = "free_y",ncol = 1) +  # 按variable分面，y轴自由
      labs(title = "Boxplot by MajorCluster and Variable",
           x = "Major Cluster",
           y = "Value") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),      # 删除所有网格线
        legend.position = "none"  # 隐藏图例，因为x轴已经显示了
      ) +
      scale_fill_manual(values = colors_set3_named)
    
    l1 <- length(unique(Idents(sce2)))
    ggsave(p1,file = "three_index.pdf",width = l1 /2 + 3, height = 6)
    ggsave(p1,file = "three_index.png",width = l1 /2 + 3, height = 6)
    
    
    diversity1 <- StartracDiversity(sce2, 
                                    type = j, #好像没用？默认Idents(sce2)？
                                    exportTable = T,
                                    group.by = "orig.ident")
    z <- as.data.frame(table(diversity1$group)) #每个样本一个值。
    #https://www.nature.com/articles/s41586-018-0694-x
    #expansion (expa), migration (migr) and transition (tran)

    colnames(diversity1) <- c("id","Cluster","Migration","Transition","Expansion")
    diversity2 <- diversity1 %>%    
      pivot_longer(cols = -c(id,Cluster), names_to = "variable", values_to = "value")
    diversity2 <- na.omit(diversity2)
    meta2 <- meta %>% select(id, surgery, effect)
    diversity3 <- diversity2 %>% left_join(meta2,by = c("id" = "id"))
    #diversity3$group <- paste0(diversity3$surgery,"_",diversity3$effect)
    #diversity3$group <- factor(diversity3$group, levels = order1)
    diversity3$value <- as.numeric(diversity3$value)
    # if (j != "subtype"){
    #   diversity3$Cluster <- paste0("C",diversity3$Cluster)
    # }

    # diversity4 <- diversity3[diversity3$Cluster == "CD8 Tcm",] #C6,7,8,13,16,23都没有想要的结果。
    # p1 <- ggplot(diversity4, aes(x = group, y = value)) +
    #   geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
    #   scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
    #   #geom_point(size = 3,color='red') +
    #   geom_point() +  #绘制散点
    #   facet_wrap( ~ variable,scales = "free_y",ncol = 3) +
    #   theme(panel.grid = element_blank(),
    #         axis.line = element_line(colour = 'black', linewidth = 1),
    #         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    #         panel.background = element_blank(),
    #         plot.title = element_text(size = 15, hjust = 0.5),
    #         plot.subtitle = element_text(size = 15, hjust = 0.5),
    #         axis.text = element_text(size = 15, color = 'black'),
    #         axis.title = element_text(size = 15, color = 'black')) +
    #   labs(x = '', y = 'Relative abundance (%)')+
    #   stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
    # p1
    # n1 <- length(unique(data3$variable))
    # n1 <- sqrt(n1)
    # ggsave(p1, filename = paste0("three_index.png"),width = n1 * 2 + 4,height = n1 * 3)
    # ggsave(p1, filename = paste0("three_index.pdf"),width = n1 * 2 + 4,height = n1 * 3)
    #################subTCELL################################################
    if (j == "subtype"){
      id3 <- c("CD8 Tc","CD8 Tpex","CD4 Treg","CD8 Tex-term")
    }else{
      id3 <- c("CD4 Treg","CD4 Tcm","CD8 Tex")
    }

    scmeta <- sce2@meta.data
    
    for (k in id3) {
      dir.create(k)
      setwd(k)
      diversity4 <- diversity3[diversity3$Cluster == k,] #C6,7,8,13,16,23都没有想要的结果。
      diversity4 <- diversity4 %>%
        mutate(number = str_extract(id, "^[^_]+"))
      
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      diversity4$group <- meta$group[match(diversity4$id,meta$id)]
      diversity4$group <- factor(diversity4$group,levels = order1)
      
      p1 <- ggplot(diversity4, aes(x = group, y = value)) +
        geom_boxplot(aes(color = group), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') +
        geom_line(aes(group = number), color = 'gray', lwd = 0.1)  +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ variable,scales = "free_y",ncol = 3) +
        theme(panel.grid = element_blank(),
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              plot.subtitle = element_text(size = 15, hjust = 0.5),
              axis.text = element_text(size = 15, color = 'black'),
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = '', y = 'Relative abundance (%)')+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1
      n1 <- 8
      n1 <- sqrt(n1)
      ggsave(p1, filename = paste0("three_effect_index.png"),width = n1 * 2 + 4,height = n1 * 1.5)
      ggsave(p1, filename = paste0("three_effect_index.pdf"),width = n1 * 2 + 4,height = n1 * 1.5)
      
      id4 <- rownames(scmeta)[scmeta$cluster == k]
      combined.TCR.sub <- combined.TCR
      for (k2 in 1:length(combined.TCR.sub)) {
        tcr_data <- combined.TCR.sub[[k2]]
        # 找出在 sce 中存在的 barcode
        keep <- tcr_data$barcode %in% id4
        combined.TCR.sub[[k2]] <- tcr_data[keep, ]
      }
      diversity <- clonalQuant(combined.TCR.sub, 
                               cloneCall="strict", 
                               chain = "both", 
                               scale = TRUE)
      diversity <-  diversity$data
      diversity5 <- diversity %>% left_join(meta2,by = c("values" = "id"))
      
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      diversity5$group <- meta$group[match(diversity5$values,meta$id1)]
      diversity5$group <- factor(diversity5$group,levels = order1)
      diversity5 <- diversity5 %>% mutate(number = str_extract(values, "^[^_]+"))
      p1 <- ggplot(diversity5, aes(x = group, y = scaled)) +
        geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_point() +  #绘制散点
        geom_line(aes(group = number), color = 'gray', lwd = 0.5) +  #配对样本间连线
        #facet_wrap( ~ variable,scales = "free_y",ncol = 3) +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        #labs(x = '', y = )+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1
      n1 <- 1
      #n1 <- length(unique(data3$variable))
      n1 <- sqrt(n1)
      ggsave(p1, filename = paste0("diversity1_effect.png"),width = n1 * 2 + 4,height = n1 * 3.5)
      ggsave(p1, filename = paste0("diversity1_effect.pdf"),width = n1 * 2 + 4,height = n1 * 3.5)
      
      clonalHomeostasis <- clonalHomeostasis(combined.TCR.sub, 
                                             cloneCall = "gene")
      clonalHomeostasis <- clonalHomeostasis$data
      clonaldata <- clonalHomeostasis %>% left_join(meta2,by = c("Var1" = "id"))
      
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      clonaldata$group <- meta$group[match(clonaldata$Var1,meta$id1)]
      clonaldata$group <- factor(clonaldata$group,levels = order1)
      clonaldata <- clonaldata %>% mutate(number = str_extract(Var1, "^[^_]+"))
      p1 <- ggplot(clonaldata, aes(x = group, y = value)) +
        geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_point() +  #绘制散点
        geom_line(aes(group = number), color = 'gray', lwd = 0.5) +  #配对样本间连线
        facet_wrap( ~ Var2,scales = "free_y",ncol = 3) +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        #labs(x = '', y = 'Relative abundance')+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1
      
      n1 <- length(unique(clonaldata$Var2))
      n1 <- sqrt(n1)
      ggsave(p1, filename = paste0("clonalHomeostasis_effect.png"),width = n1 * 2 + 4,height = n1 * 3.5)
      ggsave(p1, filename = paste0("clonalHomeostasis_effect.pdf"),width = n1 * 2 + 4,height = n1 * 3.5)
      clonalOverlap(combined.TCR.sub, 
                    cloneCall = "strict", 
                    method = "morisita")
      setwd("../")
    }
    #bulk TCR 结果对应。
    # scmeta1 <- sce2@meta.data[sce2@meta.data$cluster == "C8",]
    # tumor_R_post_pre <- read.delim("../../../../KC2023/TCR/result2/no-pair/tumor_good_post_pre/picture1/nucleic.txt")
    # split_list <- strsplit(scmeta1$CTnt, ";|_")
    # z <- unlist(split_list)
    # sum(tumor_R_post_pre$type %in% z)
    # z1 <- tumor_R_post_pre$type[tumor_R_post_pre$type %in% z]
    # pattern <- paste(z1, collapse = "|")
    # z2 <- scmeta1[grepl(pattern, scmeta1$CTnt),]
    # z3 <- tumor_R_post_pre[tumor_R_post_pre$type %in% z1,]
    # 
    # tumor_post_pre <- read.delim("../../../KC2023/TCR/result2/no-pair/tumor_post_pre/picture1/nucleic.txt")
    # tumor_post_pre <- tumor_post_pre[tumor_post_pre$Post > tumor_post_pre$Pre,]
    # split_list <- strsplit(scmeta1$CTnt, ";|_")
    # z <- unlist(split_list)
    # sum(tumor_post_pre$type %in% z)
    # z1 <- tumor_post_pre$type[tumor_post_pre$type %in% z]
    # pattern <- paste(z1, collapse = "|")
    # z2 <- scmeta1[grepl(pattern, scmeta1$CTnt),]
    # z3 <- tumor_post_pre[tumor_post_pre$type %in% z1,]
    #########grouped#######################
    data2 <- data1_g[[i]]
    data3 <- data2@meta.data
    data3$cluster <- data3[,j]
    # if (j != "subtype") {
    #   data3$cluster <- paste0("C",data3$RNA_snn_res.1.5)
    # }
    z <- table(data3$cluster,data3$orig.ident,data3$cloneSize)
    data4 <- data.frame(z)
    data4 <- data4[data4$Var3 != "None ( < X <= 0)",]
    
    colnames(data4) <- c("cluster","id","class","freq")
    
    data5 <- data4 %>%
      group_by(id,cluster) %>%
      mutate(percent = freq / sum(freq) * 100) %>%
      ungroup()
    # p1 <- ggplot(data5[data5$id == id[11],], aes(x = cluster, y = percent)) +
    #   geom_bar(aes(fill = class) , stat = "identity") +
    #   theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black',linewidth = 1),
    #         # 旋转x轴标签45度
    #         legend.key = element_rect(color = "grey70", linewidth = 0.3),  # 浅灰色细边框
    #         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    #   labs(y = " ", fill = NULL)+labs(y = 'Relative proportion (%)') +
    #   scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)))  # 底部不留空，顶部留5%空间
    # 
    # p1
    # l1 <- length(unique(data4$type))
    # ggsave(p1, file = "whole_type.pdf",width = (l1 / 2) + 4, height = 6)
    
    data6 <- data5 %>% left_join(meta2,by = c("id" = "id"))
    for (k in id3) {
      sce3_shared <-  sce2_shared[,sce2_shared$cluster == k]
      z <- as.data.frame(table(sce3_shared$orig.ident,sce3_shared$shared,
                               sce3_shared$cloneSize))
      z <- z %>%
        group_by(Var1,Var2) %>%
        mutate(percent = Freq / sum(Freq) * 100) %>%
        ungroup()
      z1 <- z[z$Var2 == "Shared",]
      ggplot(z1, aes(x = Var1, y = percent, fill = Var3)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7)
      z$group <- meta$group[match(z$Var1,meta$id)]
      z$number <- meta$number[match(z$Var1,meta$id)]
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      z$group <- factor(z$group,levels = order1)
      
      z1 <- z[z$group %in% c("Post_R","Post_NR"),]
      p1 <- ggplot(z1, aes(x = Var2, y = percent)) +
        geom_boxplot(aes(color = Var1), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_line(aes(group = number), color = 'gray', lwd = 0.1) +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ Var3,scales = "free_y") +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = "", y = "Percent (%)")+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1 
      # ggsave(p1, file = paste0(k,"/shared_un_four.pdf"),width = 8, height = 5)
      # ggsave(p1, file = paste0(k,"/shared_un_four.png"),width = 8, height = 6)
      
      z1 <- z[z$Var2 == "Shared",]
      p1 <- ggplot(z1, aes(x = group, y = percent)) +
        geom_boxplot(aes(color = group), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_line(aes(group = number), color = 'gray', lwd = 0.1) +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ Var3,scales = "free_y") +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = "", y = "Percent (%)")+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1 
      ggsave(p1, file = paste0(k,"/shared_four.pdf"),width = 8, height = 5)
      ggsave(p1, file = paste0(k,"/shared_four.png"),width = 8, height = 6)
      
      z1 <- z[z$Var2 != "Shared",]
      p1 <- ggplot(z1, aes(x = group, y = percent)) +
        geom_boxplot(aes(color = group), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_line(aes(group = number), color = 'gray', lwd = 0.1) +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ Var3,scales = "free_y") +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = "", y = "Percent (%)")+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1 
      ggsave(p1, file = paste0(k,"/shared_un_four.pdf"),width = 8, height = 5)
      ggsave(p1, file = paste0(k,"/shared_un_four.png"),width = 8, height = 6)
      
      
      
      data7 <- data7 %>% mutate(number = str_extract(id, "^[^_]+"))
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      data7$group <- meta$group1[match(data7$id,meta$id)]
      data7$group <- factor(data7$group,levels = order1)
      p1 <- ggplot(data7, aes(x = group, y = percent)) +
        geom_boxplot(aes(color = group), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_line(aes(group = number), color = 'gray', lwd = 0.1) +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ class,scales = "free_y") +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = "", y = "Percent (%)")+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1 
      ggsave(p1, file = paste0(k,"/four_effect.pdf"),width = 8, height = 5)
      ggsave(p1, file = paste0(k,"/four_effect.png"),width = 8, height = 6)
      
      
      
      data7 <- data6[data6$cluster == k,]
      data7 <- data7 %>% mutate(number = str_extract(id, "^[^_]+"))
      comp1 <- list(c("Post_R","Pre_R"),c("Post_R","Post_NR"),c("Post_NR","Pre_NR"),c("Pre_R","Pre_NR"))
      order1 <- c("Pre_NR","Post_NR","Pre_R","Post_R")
      color1 <- c("#4682B4","#87CEFA","#BC8F8F","#D2B48C")
      data7$group <- meta$group[match(data7$id,meta$id)]
      data7$group <- factor(data7$group,levels = order1)
      p1 <- ggplot(data7, aes(x = group, y = percent)) +
        geom_boxplot(aes(color = group), show.legend = F, 
                     width = 0.6,outlier.colour = "black",
                     outlier.size = 0.5) +  #箱线图
        scale_fill_manual(values = c("#1F77B4", "#2CA02C","#D62728", "#FF7F0E")) +  #设置颜色
        #geom_point(size = 3,color='red') + 
        geom_line(aes(group = number), color = 'gray', lwd = 0.1) +  #配对样本间连线
        geom_point(size = 0.5) +  #绘制散点
        facet_wrap( ~ class,scales = "free_y") +
        theme(panel.grid = element_blank(), 
              axis.line = element_line(colour = 'black', linewidth = 1),
              axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
              panel.background = element_blank(), 
              plot.title = element_text(size = 15, hjust = 0.5), 
              plot.subtitle = element_text(size = 15, hjust = 0.5), 
              axis.text = element_text(size = 15, color = 'black'), 
              axis.title = element_text(size = 15, color = 'black')) +
        labs(x = "", y = "Percent (%)")+
        stat_compare_means(method = "wilcox.test",comparisons = comp1,vjust = 1.5)
      p1 
      ggsave(p1, file = paste0(k,"/four_effect.pdf"),width = 8, height = 5)
      ggsave(p1, file = paste0(k,"/four_effect.png"),width = 8, height = 6)
    }
    setwd("../")
  }
  setwd("../")
}
