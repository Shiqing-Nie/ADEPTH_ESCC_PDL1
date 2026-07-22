rm(list=ls())
#srun --partition=gpu --nodes=1 --ntasks=1 --cpus-per-task=10 --mem=64000 --time=12:00:00 --pty /bin/bash
options(stringsAsFactors = F) 

source('scRNA_scripts/lib.R')
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(future)
#revise
plan(multisession, workers = 1)  
options(future.globals.maxSize = 100000 * 1024^2)

###### step1:导入数据 ######   
dir <- list.dirs(path = "rna_cellranger/",full.names = FALSE, recursive = FALSE)
length(unique(dir)) == length(dir)
fs <- sub("^([^_]+_[^_]+)_.*$", "\\1", dir)
library(tidyverse)
samples <- fs
samples 
sceList = lapply(samples,function(pro){ 
  print(pro)  
  tmp = Read10X(paste0("./rna_cellranger/",pro,"_out")) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 200 )
  # sce =CreateSeuratObject(counts =  ct ,
  #                         project =  pro)
  return(sce)
})
names(sceList) <- samples
#判断是否有丰度为0的gene
geneid <- list()
for (i in names(sceList)) {
  z <- sum(rowSums(sceList[[i]]@assays$RNA$counts) == 0)
  if (z > 0) {
    print("origin counts: have 0 abudnace gene")
  }
  geneid[[i]] <- rownames(sceList[[i]])
}

dir.create("1-double")
setwd("1-double/")

save(geneid,file = "1-double/geneid.Rdata")
#Go to the Windows computer, mapping
load("geneid2.Rdata")
for (i in names(geneid2)) {
  id1 <- geneid2[[i]]
  if (length(id1) != length(rownames(sceList[[i]]))) {
    print("length was wrong")
  }
  a <- sum(id1 == rownames(sceList[[i]]))
  print(a)
  rownames(sceList[[i]]) <- id1
  sceList[[i]] <-  sceList[[i]][!grepl("^ENSG", rownames(sceList[[i]])),]
}

samples
#名称不对，修改,注意检查顺序
samples <- c("HDH_post","HDH_pre","LGC_post","LGC_pre","LGY_post","LGY_pre",
             "LRX_post","LRX_pre","TQH_post","TQH_pre","WDF_post","WDF_pre",
             "XCH_post","XCH_pre","ZCX_post","ZCX_pre","ZLQ_post","ZLQ_pre",
             "ZZS_post","ZZS_pre")
names(sceList) <- samples
for (i in samples) {
  sceList[[i]]$orig.ident <- i
}

zc <- do.call(rbind,lapply(sceList, dim))
rownames(zc) <- samples
colnames(zc) <- c("nfeature","ncell")


source("../scRNA_scripts/double.R")
for (i in 1:length(sceList)) {
  sce <- sceList[[i]]
  doublet <- basic_double(sce)
  #sce2 <- subset(sce, !(rownames(sce@meta.data) %in% z))
  a <- sum(rownames(sce@meta.data) != rownames(doublet))
  if (a > 0 ){
    cat("\033[31m", "doublet order was wrong", "\033[0m\n")
    break
  }else{
    sce$doublet <- doublet[,ncol(doublet)]
    sceList[[i]] <- sce
  }
}


setwd("../")
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  ) 
#names(sce.all@assays$RNA@layers) <- paste0("counts.",samples)

# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
dim(sce.all[["RNA"]]$counts )

#as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
# fivenum(sce.all$nFeature_RNA)
# table(sce.all$nFeature_RNA>800) 
# sce.all=sce.all[,sce.all$nFeature_RNA>800]
# sce.all

library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)

meta <- read.table("meta.txt",header = T, row.names = NULL)

sce.all$sample <- sce.all$orig.ident
sp='human'
# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
samples <- meta$id[match(samples,meta$id)]
a <- sum(meta$id[match(rownames(zc),meta$id)] != samples)
if (a > 0){
  print("something was wrong")
}
rownames(zc) <- samples
t0 <- data.frame(matrix(0,nrow = length(samples), ncol = 2))
rownames(t0) <- samples
colnames(t0) <- c("nfeature-filter","ncell-filter")

for (i in samples){
  print(i)
  t1 <- rownames(sce.all.filt@meta.data)[sce.all.filt@meta.data$sample == i]
  t2 <- rownames(sce.all.filt@assays[["RNA"]]@cells) %in% t1
  t3 <- sce.all.filt@assays[["RNA"]]@layers$counts
  t3 <- t3[,t2,drop = F]
  t3 <- t3[rowSums(t3) > 0,,drop = F]
  t4 <- dim(t3)
  t0[i,] <- t4
  rm(t1,t2,t3,t4)
}
z <- sum(rownames(zc) != rownames(t0))
if (z > 0){
  print("sample order is wrong")
}
zc <- cbind(zc,t0)
write.table(data.frame(ID=rownames(zc),zc),"counts.txt",row.names=F, quote = F, sep = "\t")
id1 <- c("ZCX_post","ZCX_pre")
sce.all.filt <- sce.all.filt[,!c(sce.all.filt$orig.ident %in% id1)]
save(sce.all.filt,file = "filter.Rdata")
setwd('../')
getwd()

###### step3: harmony整合多个单细胞样品 ######
set.seed(10086)  #注意种子，然每次都不一样
table(sce.all.filt$sample)
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  # sce.all.filt <- NormalizeData(sce.all.filt, 
  #                               normalization.method = "LogNormalize",
  #                               scale.factor = 1e4) 
  # sce.all.filt <- FindVariableFeatures(sce.all.filt)
  # sce.all.filt <- ScaleData(sce.all.filt)
  sce.all.filt[["RNA"]] <- split(sce.all.filt[["RNA"]], f = sce.all.filt$orig.ident)
  #should split counts first
  sce.all.filt <- SCTransform(sce.all.filt, 
              vars.to.regress = c("percent_mito", "nCount_RNA"))
  
  # 设置默认assay为SCT
  DefaultAssay(sce.all.filt) <- "SCT"
  sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
  p1 <- ElbowPlot(sce.all.filt,ndims = 50)
  ggsave(plot=p1, filename="ElbowPlot.pdf",width = 10, height = 6)
  #sce.all.filt <- IntegrateLayers(object = sce.all.filt, method = CCAIntegration, normalization.method = "SCT", verbose = F)
  #choice dim0 according to p1
  sce.all.int = run_harmony(input_sce = sce.all.filt, dim0 = 30,assay = "SCT") #50也没啥大问题 #16更好看
  setwd('../')
}
#sce.all.int <- readRDS("2-harmony/sce.all_int.rds")
##########添加污染信息##############
# meta2 <- sceList[[1]]@meta.data
# rownames(meta2) <- paste0(meta2$orig.ident,"_",rownames(meta2))
# for (i in 2:length(sceList)) {
#   meta3 <- sceList[[i]]@meta.data
#   rownames(meta3) <- paste0(meta3$orig.ident,"_",rownames(meta3))
#   if (ncol(meta3) >4){
#     meta2 <- rbind(meta2,meta3)
#   }
# }
# meta3 <- meta2[rownames(sce.all.int@meta.data),]
# sum(is.na(meta3$Contamination))
# sce.all.int$Contamination <- meta3$Contamination
# 
# DimPlot(sce.all.int,group.by = "orig.ident",reduction = "umap")
# FeaturePlot(sce.all.int, features = "Contamination", reduction = "umap")
# dim(sce.all.int)
# low_con_scobj <- sce.all.int[,sce.all.int$Contamination < 0.01] #保留小于等于0.2的细胞
# DimPlot(low_con_scobj,group.by = "orig.ident",reduction = "umap")
# dim(low_con_scobj)
# #############################

#####################################
sce.all.int@meta.data$SCT_snn_res.0.1
colnames(sce.all.int@meta.data)
DimPlot(sce.all.int,group.by = "sample",reduction = "umap")

DimPlot(sce.all.int,group.by = "SCT_snn_res.0.05",reduction = "umap",raster = F,
        label = T,repel = T)
DimPlot(sce.all.int,group.by = "SCT_snn_res.0.1",reduction = "umap",raster = F,
        label = T,repel = T)
DimPlot(sce.all.int,group.by = "SCT_snn_res.1",reduction = "umap",raster = F,
        label = T,repel = T)

# meta1 <- read.delim("../rename.txt")
# input_sce.all@meta.data$orig.ident <- meta1$id[match(input_sce.all@meta.data$orig.ident,meta1$ID)]
# p1 <- DimPlot(input_sce.all, reduction = "umap", group.by = "orig.ident") + 
#   ggtitle("Harmony")
# ggsave(p1, filename="harmony.pdf",width = 8,height = 7)
# ggsave(p1, filename="harmony.png",width = 8,height = 7)
#######下面代码也可以不运行
###### step4:  看标记基因库 ######
# 原则上分辨率是需要自己肉眼判断，取决于个人经验
# 为了省力，我们直接看 0.1 和 0.8 即可
table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$SCT_snn_res.0.05) 
table(sce.all.int$SCT_snn_res.1) 


#p2_tree=clustree(sce.all.int@meta.data, prefix = "SCT_snn_res.")
#ggsave(plot=p2_tree, filename="2-harmony/Tree_diff_resolution.pdf",height = 8, width = 8)
#根据该树判断分为多少个类群合适"2-harmony/Tree_diff_resolution.pdf"
p1 <- DimPlot(sce.all.int,group.by = "doublet")
ggsave(p1,filename = "1-double/whole.pdf", height = 5, width = 6)
#0.1合适分大群
getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "SCT_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
sce.all.int <- PrepSCTFindMarkers(sce.all.int)
cluster.markers <- FindAllMarkers(sce.all.int, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
save(cluster.markers, file = "markers.Rdata")
top_markers <- cluster.markers[cluster.markers$p_val_adj < 1e-10,]

top20_markers <- top_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20)    
# top10_markers <- top10_markers %>%
#   group_by(cluster) %>%
#   summarise(marker_genes = paste(gene, collapse = ", ")) %>%
#   ungroup()
write.table(top20_markers,"top20_marker_FindAllMarkers.txt",row.names = F, sep = "\t", quote = F)

table(sce.all.int@active.ident) 


source('../scRNA_scripts/esophagus-marker.R')
genes_to_check = esophagus_markers_list
dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
genes_to_check = lapply(genes_to_check, function(x) x[!x %in% dup])
p1 <- DotPlot(sce.all.int,features = genes_to_check) +
  theme(axis.text.x=element_text(angle=45,hjust = 1))
w=length( unique(unlist(genes_to_check)) )/5+6;w
z_width <- length(unique(Idents(sce.all.int)))
ggsave(paste('marker.pdf'), height =  z_width / 2 ,width  = w + 4)

p1 <- DimPlot(sce.all.int,label = T)
ggsave(p1,filename = "resolution.pdf", height = 5, width = 6)
setwd('../') 
