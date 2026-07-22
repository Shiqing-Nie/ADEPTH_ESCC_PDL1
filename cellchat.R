library(CellChat)
library(patchwork)
library(Seurat)

data.dir <- '13-cellchat'
dir.create(data.dir)
setwd(data.dir)
sce <- readRDS("../../3-auto_annotate/sce.clean.rds")

meta2 <- read.delim("../../meta.txt")

group1 <- unique(meta2$group)
load("../finaltype/sce_meta.Rdata")

library(future)
options(future.seed = TRUE)
options(future.globals.maxSize = 10 * 1024^3)  # 2GB
plan("multisession", workers = 4) # 并行处理

a <- sum(rownames(sce_meta) != rownames(sce@meta.data))
if (a > 0){
  print("something was wrong")
}
meta <- sce_meta
meta <- meta %>% dplyr::select(sample,celltype)
#meta3 <- meta %>% left_join(meta2, by = c("orig.ident" = "id"))
meta$group <- meta2$group[match(meta$sample,meta2$id)]

chat1 <- c("celltype2","celltype3")
sce_meta$celltype3[sce_meta$celltype == "Epithelial"] <- sce_meta$celltype2[sce_meta$celltype == "Epithelial"]

list0 <- c("all","surgery","effect","four")
j <- chat1[2]
for (j in chat1[2]) {
  for (k in list0) {
    dir.create(k)
    setwd(k)
    sce$celltype2 <- sce$celltype
    sce$celltype2 <- sce_meta[,j]
    
    a <- sum(colnames(sce) != rownames(sce_meta))
    if (a > 0){
      print("something was wrong")
    }
    meta0 <- meta
    if (k == "all"){
      meta0$group <- "all"
    }
    if (k == "surgery"){
      meta0$group[meta0$group %in% c("Pre_R","Pre_NR")] <- "Pre"
      meta0$group[meta0$group %in% c("Post_R","Post_NR")] <- "Post"
    }
    if (k == "effect"){
      meta0$group[meta0$group %in% c("Pre_R","Post_R")] <- "R"
      meta0$group[meta0$group %in% c("Pre_NR","Post_NR")] <- "NR"
    }
    group2 <- unique(meta0$group)
    for (i in group2) {
      cell.use = rownames(meta0)[meta0$group == i & meta0$celltype %in% c("Myeloid","Mast","T_NK","Epithelial")] #
      # B. 为 CelChat 分析子集输入数据
      sce.use <- sce[,cell.use]
      sce.use <- sce.use[,sce.use$celltype2 != "contamination"]
      sce.use$samples <- as.factor(sce.use$sample)
      cellchat <- createCellChat(object = sce.use, group.by = "celltype2",assay = "SCT")  #从seruat对象直接构建更方便。
      
      #每组数据需独立运行以下步骤：
      #(1) 配体-受体数据库预处理
      CellChatDB <- CellChatDB.human# 如果分析小鼠数据，请使用 CellChatDB.mouse
      showDatabaseCategory(CellChatDB)
      #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")# 使用分泌信号
      cellchat@DB <- CellChatDB
      #(2) 计算细胞通讯概率
      # 使用 CellChatDB 的一个子集进行细胞间通讯分析
      cellchat <- subsetData(cellchat)  # 过滤低表达基因
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      
      cellchat <- computeCommunProb(cellchat, type = "triMean")# 计算通讯概率
      #(3) 整合通路并过滤
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      #用户可以过滤掉某些细胞群中细胞数量较少的细胞间通讯。默认情况下，每个细胞群中用于细胞间通讯的最小细胞数量为 10。
      ####在信号通路水平推断细胞间通讯####
      cellchat <- computeCommunProbPathway(cellchat)
      #计算聚合的细胞间通讯网络
      cellchat <- aggregateNet(cellchat)
      saveRDS(cellchat, file =  paste0(i,".rds"))
    }
    setwd("../")
  }
}
plan("multisession", workers = 1) # 并行处理
for (l1 in list0[3]) {
  setwd(l1)
  for (j in chat1[2]) {
    #setwd(j)
    cellchat.NR <- readRDS("NR.rds")
    cellchat.R <- readRDS("R.rds")
    object.list <- list(NR = cellchat.NR,
                        R = cellchat.R)
    
    # 对object.list中的每个对象执行（假设object.list已定义）
    for (i in 1:length(object.list)) {
      object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
    }
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    #> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
    cellchat
    data1 <-  subsetCommunication(cellchat) #所有的结果
    save(data1,file = "data1.Rdata")
    
    id0 <- names(data1)
    for (k in id0) {
      data2 <- data1[[k]]
      id1 <- c("Mast-HLA-A","Macro-SPARC")
      dir.create(k)
      setwd(k)
      for (i in id1) {
        mcaro_tpex <- data2[data2$source == i & data2$target == "CD8 Tpex",]
        source <- as.data.frame(table(mcaro_tpex$receptor)) %>% filter(Freq > 1)
        mcaro_tpex <- mcaro_tpex[mcaro_tpex$receptor %in% source$Var1,]
        
        #mcaro_tpex$prob <- log(mcaro_tpex$prob)
        p1 <- ggplot(mcaro_tpex, aes(x = ligand, y = receptor)) +
          geom_tile(aes(fill = prob), color = "gray") +
          scale_fill_gradient(low = "orange", high = "red") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid = element_blank(),      # 删除所有网格线
                panel.background = element_blank()) +
          labs(x = "Sender", y = "Ligand-Receptor Pair", fill = "Regulatory potential",
               title = paste0(i," → CD8 Tpex (",k,")"))
        p1
        ggsave(p1, file = paste0(i,"_tpex.pdf"), height = 3,width = 8)
      } 
      setwd("../")
    }
    
    ##########详细结果########################
    # # 查看细胞群名称（对应行列）
    # rownames(cellchat@net$Pre_NR$prob)  # 源细胞群
    # colnames(cellchat@net$Pre_NR$prob)  # 靶细胞群
    # # 查看配体-受体对名称（第三维）
    # dimnames(cellchat@net$Pre_NR$prob)[[3]][1]  # 或直接 head(cellchat@LR$LRsig)
    # #提取特定细胞群间的通信
    # # 示例：提取细胞群 "C0" -> "Bcell" 的所有配体-受体概率
    # prob_fibro_macro <- cellchat@net$Pre_NR$prob["C0", "C22", ]
    # head(sort(prob_fibro_macro, decreasing = TRUE))  # 按概率排序
    # # 示例：查看配体-受体对 "TGFB1_TGFBR1_TGFBR2" 在所有细胞群中的概率
    # lr_name <- "TGFB1_TGFBR1_TGFBR2"
    # prob_lr <- cellchat@net$Pre_NR$prob[, , lr_name]
    # head(prob_lr)
    # # 将三维数组转换为二维表格（长格式）
    # library(reshape2)
    # df_prob <- melt(cellchat@net$Pre_NR$prob)
    # colnames(df_prob) <- c("source", "target", "LR_pair", "probability")
    # head(df_prob[order(df_prob$probability, decreasing = TRUE), ])
    # 
    # data2 <- data1$Pre_NR
    # sum(data2$source == "Bcell" & data2$target == "Bcell")
    # sum(data2$source == "Bcell" & data2$target == "C0")
    
    ##################################
    dir.create("picture")
    setwd("picture/")
    ############circle###################################
    if (j == "celltype3"){
      levels(object.list$NR@idents)
      weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
      # png(file = paste0("CellChat_number.png"), 
      #     width = 3000, height = 1500, res = 300) # 可根据需要调整尺寸和分辨率:cite[6]

      
      sources_li <- c("Macro-SPARC","Mast-HLA-A")
      for (sources1 in sources_li) {
        pdf(file = paste0("CellChat_number_","_",sources1,".pdf"),
            width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]
        #sources1 <- c("Mast-HLA-A","Macro-SPARC")
        targets1 <- c("CD8 Tc","CD8 Tcm-ANXA1","CD8 Tcm-NFKBIA","CD8 Tem-HSPA1A",
                      "CD8 Tem-IFI6","CD8 Proliferation","CD8 Temra","CD8 Tpex","CD8 Tex-term")
        par(mfrow = c(1,2), xpd=TRUE)
        for (i in c(1,2)) {
          netVisual_circle(object.list[[i]]@net$count, 
                           weight.scale = T, label.edge= T, 
                           # sources.use = c("CD8 Temra","cDC1","Malignant cell"),
                           # targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
                           #                 "CD4 Treg-RTKN2",
                           #                 "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
                           #                 "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
                           sources.use = sources1,
                           targets.use = targets1,
                           # 修改顶点标签字体大小
                           vertex.label.cex = 0.7,  # 增加节点标签字体大小（默认约0.8）
                           # 修改边标签字体大小
                           edge.label.cex = 0.8,    # 边标签字体大小
                           # 修改图例字体大小（如果有图例）
                           #legend.cex = 1.0,# 图例字体大小
                           remove.isolate = T,
                           #color.use = my_colors,
                           edge.weight.max = 100,
                           title.name = paste0("Number of interactions - ", names(object.list)[i]))
        }
        dev.off()
        
        # png(file = paste0("CellChat_stregth.png"),
        #     width = 2000, height = 2000, res = 300) # 可根据需要调整尺寸和分辨率:cite[6]
        pdf(file = paste0("CellChat_strenth_","_",sources1,".pdf"),
            width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]
        par(mfrow = c(1,2), xpd=TRUE)
        for (i in c(1,2)) {
          netVisual_circle(object.list[[i]]@net$weight, 
                           weight.scale = T, label.edge= T, 
                           # sources.use = c("CD8 Temra","cDC1","Malignant cell"),
                           # targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
                           #                 "CD4 Treg-RTKN2",
                           #                 "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
                           #                 "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
                           sources.use = sources1,
                           targets.use = targets1,
                           remove.isolate = T,
                           edge.weight.max = 100,
                           title.name = paste0("Strength of interactions - ", names(object.list)[i]))
        }
        dev.off() 
      }
      # pdf(file = paste0("CellChat_DC_number.pdf"), 
      #     width = 10, height = 6) # 可根据需要调整尺寸和分辨率:cite[6]
      # par(mfrow = c(1,2), xpd=TRUE)
      # for (i in c(3,4)) {
      #   netVisual_circle(object.list[[i]]@net$count, 
      #                    weight.scale = T, label.edge= T, 
      #                    sources.use = c("CD8 Temra","cDC1","Malignant cell"),
      #                    targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
      #                                    "CD4 Treg-RTKN2",
      #                                    "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
      #                                    "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
      #                    remove.isolate = T,
      #                    edge.weight.max = 100,
      #                    title.name = paste0("Number of interactions - ", names(object.list)[i]))
      # }
      # dev.off()    
      # 
      # pdf(file = paste0("CellChat_DC_strength.pdf"), 
      #     width = 10, height = 6) # 可根据需要调整尺寸和分辨率:cite[6]
      # par(mfrow = c(1,2), xpd=TRUE)
      # for (i in c(3,4)) {
      #   netVisual_circle(object.list[[i]]@net$weight, 
      #                    weight.scale = T, label.edge= T, 
      #                    sources.use = c("CD8 Temra","cDC1","Malignant cell"),
      #                    targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
      #                                    "CD4 Treg-RTKN2",
      #                                    "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
      #                                    "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
      #                    remove.isolate = T,
      #                    edge.weight.max = 100,
      #                    title.name = paste0("Strength of interactions - ", names(object.list)[i]))
      # }
      # dev.off()
      
      # png(file = paste0("CellChat_tc_number.png"), 
      #     width = 2000, height = 2000, res = 300) # 可根据需要调整尺寸和分辨率:cite[6]
      pdf(file = paste0("CellChat_tex_number.pdf"), 
          width = 10, height = 10,) # 可根据需要调整尺寸和分辨率:cite[6]
      par(mfrow = c(1,2), xpd=TRUE)
      # 只填充C0与C1-C5之间的通讯（双向）
      for (i in c(1,2)) {
        cell_groups <- c("CD8 Tpex","Macro-SPARC","Malignant cell","Mast-HLA-A")
        custom_net <- matrix(0, 
                             nrow = length(cell_groups), 
                             ncol = length(cell_groups),
                             dimnames = list(cell_groups, cell_groups))
        custom_net[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")] <- object.list[[i]]@net$count[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")]
        custom_net[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")] <- object.list[[i]]@net$count[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")]
        
        
        # 使用自定义矩阵绘图
        netVisual_circle(custom_net,
                         weight.scale = T, label.edge= T, 
                         remove.isolate = T,
                         edge.weight.max = 100,
                         title.name = paste0("Number of interactions - ", names(object.list)[i]))
      }
      dev.off()
      
      # png(file = paste0("CellChat_tc_strength.png"), 
      #     width = 2000, height = 2000, res = 300) # 可根据需要调整尺寸和分辨率:cite[6]
      pdf(file = paste0("CellChat_tex_strength.pdf"), 
          width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]
      par(mfrow = c(1,2), xpd=TRUE)
      # 只填充C0与C1-C5之间的通讯（双向）
      for (i in c(1,2)) {
        cell_groups <- c("CD8 Tpex","Macro-SPARC","Malignant cell","Mast-HLA-A")
        custom_net <- matrix(0, 
                             nrow = length(cell_groups), 
                             ncol = length(cell_groups),
                             dimnames = list(cell_groups, cell_groups))
        custom_net[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")] <- object.list[[i]]@net$weight[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")]
        custom_net[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")] <- object.list[[i]]@net$weight[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")]
        
        
        # 使用自定义矩阵绘图
        netVisual_circle(custom_net,
                         weight.scale = T, label.edge= T, 
                         remove.isolate = T,
                         edge.weight.max = 100,
                         title.name = paste0("Stregth of interactions - ", names(object.list)[i]))
      }
      dev.off()
    }
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    p <- gg1 + gg2
    p
    #ggsave(p, file = "counts_weight.pdf", width = 10, height = 6)
    #ggsave(p, file = "counts_weight.png", width = 10, height = 6)
    #默认是第二组比第一组，红色增加，蓝色降低
    #
    levels(cellchat@meta$datasets)  #查看顺序
    #netVisual_diffInteraction(cellchat, weight.scale = T,edge.width = c(0.01, 0.03))
    #edge.width.max = 0.1
    #netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    #线太粗了，不会调
    #在 CellChat 中比较两组样本（如正常皮肤 NL 和病变皮肤 LS）时，细胞数量差异确实会影响通讯网络的分析结果，但 CellChat 已通过内置的统计方法对细胞数量进行了校正。
    
    #The top colored bar plot represents the sum of each column of the absolute values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions or interaction strength between the two conditions. In the colorbar, red(or blue) represents increased(or decreased) signaling in the second dataset compared to the first one.
    gg1 <- netVisual_heatmap(cellchat,comparison = c("NR", "R")) #sources.use = c(1,4)
    #> Do heatmap based on a merged object
    gg2 <- netVisual_heatmap(cellchat, comparison = c("R", "NR"),measure = "weight")
    #In the colorbar, red(or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
    
    #> Do heatmap based on a merged object
    # pdf("counts_R_weight2.pdf", width = 10, height = 6)
    # gg1 + gg2
    # dev.off()
    #png("counts_weight2.png", width = 10, height = 6, units = "in", res = 300)
    pdf("counts_weight2.pdf", width = 10, height = 6)
    gg1 + gg2
    dev.off()
    
    gg1 <- netVisual_heatmap(cellchat,comparison = c("NR", "R")) #sources.use = c(1,4)
    
    #> Do heatmap based on a merged object
    gg2 <- netVisual_heatmap(cellchat, comparison = c("NR", "R"),measure = "weight")
    #> Do heatmap based on a merged object
    # pdf("counts_NR_weight2.pdf", width = 10, height = 6)
    # gg1 + gg2
    # dev.off()
    #png("counts_NR_weight2.png", width = 13, height = 8, units = "in", res = 300)
    #gg1 + gg2
    #dev.off()
    
    #The above differential network analysis only works for pairwise datasets. If there are more datasets for comparison, CellChat can directly show the number of interactions or interaction strength between any two cell populations in each dataset.
    
    #To better control the node size and edge weights of the inferred networks across different datasets, CellChat computes the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    
    #
    head(cellchat@idents$joint)          # 显示前几个细胞的分类
    levels(cellchat@idents$joint)        # 显示所有分类水平（因子顺序）
    table(cellchat@idents$joint)         # 统计各分类的细胞数量
    #group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))  #对应着12种细胞
    group.cellType <- levels(cellchat@idents$joint) #暂时不分组
    group.cellType <- factor(group.cellType, levels = c(group.cellType))
    object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
    par(mfrow = c(1,3), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    
    par(mfrow = c(1,3), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
    
    
    num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
    weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
    gg <- list()
    
    for (i in 1:length(object.list)) {
      gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
    }
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    # pdf("counts_weight3.pdf", width = 10, height = 5)
    # patchwork::wrap_plots(plots = gg)
    # dev.off()
    #png("counts_weight3.png", width = 10, height = 5, units = "in", res = 300)
    pdf("counts_weight3.pdf", width = 12, height = 5)
    patchwork::wrap_plots(plots = gg)
    dev.off()
    
    ###########取子集############################
    object.list_sub <- object.list
    for (su in names(object.list_sub)) {
      # 如果知道具体细胞索引
      cellchat_sub <- object.list_sub[[su]]
      keep_celltypes <- c("CD8 Tc", "CD8 Tem-HSPA1A", "CD8 Temra", "CD8 Tcm-ANXA1", 
                          "CD8 Tcm-NFKBIA", "CD8 Tex-term", "CD8 Tpex", "CD8 Tem-IFI6",
                          "CD8 Proliferation", "Macro-SPARC")
      object.list_sub[[su]] <- subsetCellChat(cellchat_sub, idents.use = keep_celltypes)
    }
    
    for (i in 1:length(object.list_sub)) {
      gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_sub[[i]], title = names(object.list_sub)[i])
    }
    pdf("counts_weight3.pdf", width = 12, height = 5)
    patchwork::wrap_plots(plots = gg)
    dev.off()
    ########################################
    
    
    #Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS.
    #gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "C0", signaling.exclude = "MIF")
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0
    #> The following `from` values were not present in `x`: 0, -1
    #gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0, 2
    #> The following `from` values were not present in `x`: 0, -1
    #patchwork::wrap_plots(plots = list(gg1,gg2))
    if (j == "celltype2"){
      id_t <- c("CD8 Tc")
      id_mye <- c("Mast")
    }else{
      id_t <- c("CD8 Tpex")
      id_mye <- c("Mast-HLA-A","Macro-SPARC")
    }
    # rankNet(object.list$Post_R, mode = "single", 
    #         measure = "weight", sources.use = id_t, 
    #         targets.use = NULL)
    
    gg1 <- rankNet(cellchat, mode = "comparison", 
                   measure = "weight", sources.use = id_t, 
                   comparison = c(1, 2),
                   targets.use = NULL, stacked = T, do.stat = TRUE)
    data1 <- gg1$data %>% group_by(name) %>%
      mutate(percent = (contribution / sum(contribution)) * 100)
    order_t <- data1 %>% filter(group == "Post_R") %>%
      arrange(desc(percent))
    data1$name <- factor(data1$name, levels = order_t$name)
    ggplot(data1, aes(x = name,y = percent, fill = group)) +  
      geom_bar(position = "stack",stat = "identity") + 
      coord_flip()
    #统计差异是否显著do.stat = TRUE
    gg2 <- rankNet(cellchat, mode = "comparison", 
                   measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
    gg2
    #默认是comparison = c(1, 2),
    ggsave(gg1,filename = "important.png",width = 6, height = 6)
    #ggsave(gg1,filename = "important.pdf",width = 6, height = 6)
    
    library(ComplexHeatmap)
    i = 1
    # combining all the identified signaling pathways from different datasets
    pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10)
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10)
    z <- ht1@matrix
    rowSums(z,na.rm = T)
    # pdf("outgoing.pdf", width = 5, height = 10)
    # ht1 + ht2
    # dev.off()
    png("outgoing.png", width = 10, height = 12, units = "in", res = 300)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    
    dir.create("more1")
    names1 <- names(object.list)
    comp1 <- list(R = c(1,2),NR = c(1,2))
    comp2 <- list(R = 2,NR = 1)
    comp3 <- list(R = "R",NR = "NR")
    source1 <- list(to_tc = c(id_mye),
                    malignant_to = c("Malignant cell"),
                    to_malignant = c(id_mye,id_t),
                    my_tc = c(id_mye),
                    tc_dc = c(id_t)
                    #treg_tc = c("CD4 Treg"),
                    #tc_treg = c("CD8 Tc")
    )
    target1 <- list(to_tc = c(id_t),
                    malignant_to = c(id_mye,id_t),
                    to_malignant = c("Malignant cell"),
                    my_tc = c(id_t),
                    tc_dc =  c(id_mye)
                    #treg_tc = c("CD8 Tc"),
                    #tc_treg = c("CD4 Treg")
    )
    picture1 <- list()
    for(k1 in names(comp1)){
      for (k2 in names(source1)) {
        gg1 <- netVisual_bubble(cellchat, sources.use = source1[[k2]], 
                                targets.use = target1[[k2]],  
                                #comparison = c(3, 4), max.dataset = 4, 
                                comparison = comp1[[k1]], max.dataset = comp2[[k1]], 
                                title.name = paste0("Increased signaling in ",comp3[[k1]]), 
                                angle.x = 45, remove.isolate = T)
        #ggsave(gg1,file = paste0("Tc_tumor_increse.pdf"),height = 5, width = 5)
        z1 <- gg1@data
        z1 <- length(unique(z1$interaction_name))
        names2 <- paste(names1[comp1[[k1]][1]],names1[comp1[[k1]][2]],sep = "_")
        picture1[[k1]][[k2]] <- gg1
        ggsave(gg1,file = paste0("more1/",names2,"_",k2,"_",comp3[[k1]],"_increse.pdf"),height = z1/3, width = 5)
      }
    }
    setwd("../")
  }
  setwd("../")
}

###基于全部数据##################################
setwd("all/")
cellchat <- readRDS("all.rds")

cellchat <- netAnalysis_computeCentrality(cellchat)
data1 <-  subsetCommunication(cellchat) #所有的结果
save(data1,file = "data1.Rdata")

###########取子集############################
object.list_sub <- cellchat
keep_celltypes <- c("CD8 Tc", "CD8 Tem-HSPA1A", "CD8 Temra", "CD8 Tcm-ANXA1", 
                    "CD8 Tcm-NFKBIA", "CD8 Tex-term", "CD8 Tpex", "CD8 Tem-IFI6",
                    "CD8 Proliferation", "Macro-SPARC")
object.list_sub <- subsetCellChat(object.list_sub, idents.use = keep_celltypes)
gg <- netAnalysis_signalingRole_scatter(object.list_sub, title = "All group")

pdf("counts_weight3.pdf", width = 6.5, height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

##########详细结果########################
dir.create("picture")
setwd("picture/")

id1 <- c("Mast-HLA-A","Macro-SPARC")
for (sources1 in id1) {
  
  pdf(file = paste0(sources1,"_CellChat_number.pdf"),
      width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]
  #sources1 <- c("Mast-HLA-A","Macro-SPARC")
  
  targets1 <- c("CD8 Tc","CD8 Tcm-ANXA1","CD8 Tcm-NFKBIA","CD8 Tem-HSPA1A",
                "CD8 Tem-IFI6","CD8 Proliferation","CD8 Temra","CD8 Tpex","CD8 Tex-term")
  netVisual_circle(cellchat@net$count, 
                   weight.scale = T, label.edge= T, 
                   # sources.use = c("CD8 Temra","cDC1","Malignant cell"),
                   # targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
                   #                 "CD4 Treg-RTKN2",
                   #                 "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
                   #                 "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
                   sources.use = sources1,
                   targets.use = targets1,
                   # 修改顶点标签字体大小
                   vertex.label.cex = 0.7,  # 增加节点标签字体大小（默认约0.8）
                   # 修改边标签字体大小
                   edge.label.cex = 0.8,    # 边标签字体大小
                   # 修改图例字体大小（如果有图例）
                   #legend.cex = 1.0,# 图例字体大小
                   remove.isolate = T,
                   #color.use = my_colors,
                   edge.weight.max = 100,
                   title.name = paste0("Number of interactions"))
  dev.off()
  
  # png(file = paste0("CellChat_stregth.png"),
  #     width = 2000, height = 2000, res = 300) # 可根据需要调整尺寸和分辨率:cite[6]
  pdf(file = paste0(sources1,"_CellChat_STRENGTH.pdf"),
      width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]
  
  netVisual_circle(cellchat@net$weight, 
                   weight.scale = T, label.edge= T, 
                   # sources.use = c("CD8 Temra","cDC1","Malignant cell"),
                   # targets.use = c("CD8 Temra","cDC1","Malignant cell","CD8 Tcm-FOS",
                   #                 "CD4 Treg-RTKN2",
                   #                 "CD8 Tcm-ANXA1","CD8 Tem-ME1","Tex-term",
                   #                 "CD8 Proliferation","CD8 Tem-RAPGEF2","Tpex"),
                   sources.use = sources1,
                   targets.use = targets1,
                   remove.isolate = T,
                   edge.weight.max = 100,
                   title.name = paste0("Strength of interactions"))
  dev.off() 
}


pdf(file = paste0("CellChat_tex_number.pdf"), 
    width = 10, height = 10,) # 可根据需要调整尺寸和分辨率:cite[6]
cell_groups <- c("CD8 Tpex","Macro-SPARC","Malignant cell","Mast-HLA-A")
custom_net <- matrix(0, 
                     nrow = length(cell_groups), 
                     ncol = length(cell_groups),
                     dimnames = list(cell_groups, cell_groups))
custom_net[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")] <- cellchat@net$count[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")]
custom_net[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")] <- cellchat@net$count[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")]


# 使用自定义矩阵绘图
netVisual_circle(custom_net,
                 weight.scale = T, label.edge= T, 
                 remove.isolate = T,
                 edge.weight.max = 100,
                 title.name = paste0("Number of interactions"))
dev.off()

pdf(file = paste0("CellChat_tex_strength.pdf"), 
    width = 10, height = 10) # 可根据需要调整尺寸和分辨率:cite[6]

cell_groups <- c("CD8 Tpex","Macro-SPARC","Malignant cell","Mast-HLA-A")
custom_net <- matrix(0, 
                     nrow = length(cell_groups), 
                     ncol = length(cell_groups),
                     dimnames = list(cell_groups, cell_groups))
custom_net[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")] <- cellchat@net$weight[c("CD8 Tpex","Macro-SPARC"), c("Malignant cell","Macro-SPARC","Mast-HLA-A")]
custom_net[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")] <- cellchat@net$weight[ c("Malignant cell","Macro-SPARC","Mast-HLA-A"),c("CD8 Tpex","Macro-SPARC")]


# 使用自定义矩阵绘图
netVisual_circle(custom_net,
                 weight.scale = T, label.edge= T, 
                 remove.isolate = T,
                 edge.weight.max = 100,
                 title.name = paste0("Stregth of interactions"))
dev.off()


gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
ggsave(gg1,file = paste0("signalingRole.pdf"),height = 6, width = 8)


#######具体通路######################
id_t <- c("CD8 Tpex")
id_mye <- c("Mast-HLA-A","Macro-SPARC")


source1 <- list(to_tc = c(id_mye),
                malignant_to = c("Malignant cell"),
                to_malignant = c(id_mye,id_t),
                my_tc = c(id_mye),
                tc_dc = c(id_t)
                #treg_tc = c("CD4 Treg"),
                #tc_treg = c("CD8 Tc")
)
target1 <- list(to_tc = c(id_t),
                malignant_to = c(id_mye,id_t),
                to_malignant = c("Malignant cell"),
                my_tc = c(id_t),
                tc_dc =  c(id_mye)
                #treg_tc = c("CD8 Tc"),
                #tc_treg = c("CD4 Treg")
)
dir.create("more1")
picture1 <- list()
for (k2 in names(source1)) {
  gg1 <- netVisual_bubble(cellchat, sources.use = source1[[k2]], 
                          targets.use = target1[[k2]],  
                          #comparison = c(3, 4), max.dataset = 4, 
                          angle.x = 45, remove.isolate = T)
  #ggsave(gg1,file = paste0("Tc_tumor_increse.pdf"),height = 5, width = 5)
  z1 <- gg1@data
  z1 <- length(unique(z1$interaction_name))
  picture1[[k2]] <- gg1
  ggsave(gg1,file = paste0("more1/",k2,"_increse.pdf"),height = z1/3, width = 5)
}

for (i in id1) {
  mcaro_tpex <- data1[data1$source == i & data1$target == "CD8 Tpex",]
  source <- as.data.frame(table(mcaro_tpex$receptor)) %>% filter(Freq > 1)
  mcaro_tpex <- mcaro_tpex[mcaro_tpex$receptor %in% source$Var1,]
  
  #mcaro_tpex$prob <- log(mcaro_tpex$prob)
  p1 <- ggplot(mcaro_tpex, aes(x = ligand, y = receptor)) +
    geom_tile(aes(fill = prob), color = "gray") +
    scale_fill_gradient(low = "orange", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),      # 删除所有网格线
          panel.background = element_blank()) +
    labs(x = "Sender", y = "Ligand-Receptor Pair",  fill = "Regulatory potential",
         title = "Macro-SPARC → CD8 Tpex")
  p1
  ggsave(p1, file = paste0(i,"_tpex.pdf"), height = 3,width = 8)
}

id2 <- id1
id3 <- rev(id2)

for (i in 1:2) {
  mcaro_tpex <- data1[data1$source == id2[i] & data1$target == id3[i],]
  source <- as.data.frame(table(mcaro_tpex$receptor)) %>% filter(Freq > 1)
  mcaro_tpex <- mcaro_tpex[mcaro_tpex$receptor %in% source$Var1,]
  
  #mcaro_tpex$prob <- log(mcaro_tpex$prob)
  p1 <- ggplot(mcaro_tpex, aes(x = ligand, y = receptor)) +
    geom_tile(aes(fill = prob), color = "gray") +
    scale_fill_gradient(low = "orange", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),      # 删除所有网格线
          panel.background = element_blank()) +
    labs(x = "Sender", y = "Ligand-Receptor Pair",  fill = "Regulatory potential",
         title = paste0(id2[i], " → ", id3[i]))
  p1
  ggsave(p1, file = paste0(id2[i], "-", id3[i],".pdf"), height = 3,width = 8)
}
