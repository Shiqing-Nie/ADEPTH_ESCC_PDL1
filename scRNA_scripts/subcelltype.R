subcelltype <- function(input_sce, dim0 = 30){
  #seuratObj <- RunHarmony(input_sce, "orig.ident")
  seuratObj <- input_sce
  names(seuratObj@reductions)
  
  # p = DimPlot(seuratObj,reduction = "umap",label=T ) 
  # ggsave(filename='umap-by-orig.ident-after-harmony',plot = p)
  # seuratObj <- RunTSNE(seuratObj, dims = 1:dim0, 
  #                      reduction = "harmony") #非常慢
  input_sce=seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:dim0) #16 
  input_sce.all=input_sce
  
  #设置不同的分辨率，观察分群效果(选择哪一个？)
  for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5)) {
    input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                               resolution = res, algorithm = 1)
  }
  colnames(input_sce.all@meta.data)
  apply(input_sce.all@meta.data[,grep("RNA_snn",colnames(input_sce.all@meta.data))],2,table)
  
  p1_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                     ggtitle("louvain_0.01"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                     ggtitle("louvain_0.1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                     ggtitle("louvain_0.2"))
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 16, height = 5)
  
  p1_dim=plot_grid(ncol = 3, DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                     ggtitle("louvain_0.8"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                     ggtitle("louvain_1"), DimPlot(input_sce.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                     ggtitle("louvain_0.3"))
  ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 16, height = 5)
  
  p2_tree=clustree(input_sce.all@meta.data, prefix = "RNA_snn_res.")
  ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 8, width = 8)
  table(input_sce.all@active.ident) 
  saveRDS(input_sce.all, "sce.all_int.rds")
  return(input_sce.all)
  
}
