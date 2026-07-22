Tcell=c("CD2","CD3D","CD3E","CD3G")
Bcell=c("CD19","CD79A","MS4A1","JCHAIN","MZB1")
Myeloid=c("CD68","LYZ","CD14","CD163","IL3RA","LAMP3","CLEC4C","TPSAB1")
Muscle=c("PLN","MYH11","ACTG2","CNN1","MYL9","TAGLN")
Epithelial=c("EPCAM","SFN","KRT5","KRT14")
Endothelial=c("VWF","PECAM1","ENG","CDH5")
Fibroblast=c("FN1","DCN","COL1A1","COL1A2","COL3A1","COL6A1")
Glandular=c("AGR2","KRT23","WFDC2")
NK=c("FCGR3A","NCAM1","KLRD1","FGFBP2","KLRC1","NKG7")
Plasma=c("JCHAIN", "MZB1","SLAMF7", "IGKC")
Mast=c("CPA3","TPSAB1","KIT")
Pericyte=c("RGS5","MCAM","ACTA2")
immune=c("PTPRC")

esophagus_markers_list = list(
  Tcell=Tcell,
  Bcell=Bcell,
  Myeloid=Myeloid,
  Muscle=Muscle,
  Epithelial=Epithelial,
  Endothelial=Endothelial,
  Fibroblast=Fibroblast,
  Glandular=Glandular,
  NK=NK,
  Plasma=Plasma,
  Mast=Mast,
  Pericyte=Pericyte,
  immune=immune
)

# Tcell=c("CD2","CD3D","CD3E","CD3G")
# Bcell=c("CD19","CD79A","MS4A1","JCHAIN","MZB1","SLAMF7", "IGKC")
# Myeloid=c("CD68","LYZ","CD14","IL3RA","LAMP3","CLEC4C","TPSAB1")  
# #Muscle=c("PLN","MYH11","ACTG2","CNN1","MYL9","TAGLN")
# Epithelial=c("EPCAM","SFN","KRT5","KRT14")
# Endothelial=c("VWF","PECAM1","ENG","CDH5")
# Fibroblast=c("FN1","DCN","COL1A1","COL1A2","COL3A1","COL6A1")
# #Glandular=c("AGR2","KRT23","WFCD2")
# NK=c("FCGR3A","NCAM1","KLRD1","FGFBP2","KLRC1")
# #Plasma=c("JCHAIN", "MZB1","SLAMF7", "IGKC") 加到B细胞
# #Mast=c("CPA3","TPSAB1","KIT")
# Pericyte=c("RGS5","MCAM","ACTA2")
# 
# esophagus_markers_list = list(
#   Tcell=Tcell,
#   Bcell=Bcell,
#   Myeloid=Myeloid,
#   #Muscle=Muscle, 
#   Epithelial=Epithelial,
#   Endothelial=Endothelial,  
#   Fibroblast=Fibroblast, 
#   #Glandular=Glandular,
#   NK=NK,
#   #Plasma=Plasma,
#   #Mast=Mast,
#   Pericyte=Pericyte
# ) 

all_strings <- unlist(esophagus_markers_list)
# 2. 统计每个字符串的出现次数
string_counts <- table(all_strings)
# 3. 找出重复的字符串（出现次数大于1）
duplicates <- names(string_counts)[string_counts > 1]

markers_list <- c(
  'esophagus_markers_list'
)

####Tcell####
CD4 <- c("CD4","IL7R")
CD8 <- c("CD8A","CD8B","GZMK")
Tcell_list <- list(CD4 = CD4, CD8 = CD8)

CD4_T <- list(navie = c("CCR7","SELL", "TCF7", "LEF1", "IL7R"),
            Memory = c("LMNA","IL7R"),
            Exhausted = c("CXCL13","CARD16"),
            Treg = c("FOXP3","IL2RA","IKZF2","CD4","CD25"),
            Helper17 = c("KLRB1","CCL20","RORA"),
            Follicular_helper = c("CXCR5","CXCL13","MAF"))
CD8_T <- list(navie = c("TCF","SELL","LEF1","CCR7","CD45RA"),
              Cytotoxic = c("GZMK","GZMN","NKG7","KLRG1","GZMA","GNLY","PRF1","GZMB","IFNG"),
              Memory = c("LMNA","CCR7"),
              Exhausted = c("ENTPD1","PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","SPRY"))
####Bcell####
Bcell_list <- list(Follicular = c("MS4A1"),
                   Naive = c("IGHD", "TCL1A", "FCER2", "CD23"), 
                   memory = c("FCRL4", "NR4A1","TNFRSF13B","TACI"),
                   Active = c("GPR183","TNFRSF13B","CD69","CD83","HLA-DRA"),
                   Proliferative = c("THBA1B"))
#,Plasma = c("IGHG1","JCHAIN")
####Myeloid####
myeloid_list <- list(Monocyte = c("LYZ","FCN1","CCL3"),
                     Macrophage	= c("C1QA","CCL18"),
                     Dendritic = c("CCL17","CCR7","CCL22"))
macrophage_list <- list(m1 = c("CD86"),
                        m2 = c("CD163","C1QA","CD68","CD206"))
Dendritic_list <- list(Activated = c("CD40", "FSCN1", "CCR7"),
                  cDC = c("CD1A","C1DC","CLEC10A"), #Conventional DC(cDC)
                  pDC = c("LILRA4"),
                  tDC = c("IDO1", "FSCN", "LAMP3"))  #Plasmacytoid DC(pDC) Tolerogenic DC(tDC)
cDC_list <- list(cDC1 = c("CLEC9A", "XCR1", "IRF8", "BATF3"),
                 cDC2 = c("FCER1A", "CD1C"))
                 

####Fibroblast#####

Fibroblast_list <- list(iCAF = c("CFD","CXCL12","DCN"),
                     myoCAF = c("ACTA2","TAGLN"))
Fibroblast_list2 <- list(normal = c("IGFBP6","PI16","DPT","CXCL12","GPX3", "DCN"),
                      cancer = c("FAP","MMP1","MMP11","POSTN"))
#Epithelial cells
QP <- c("KRT15", "CXCL14", "IL1R2", "EEF1A1", "CCL2", "SERPINE2", "ZFP36", "ZFP36L2", "DST", "MYC")
#Quiescent progenitor cells(QP) epithelial cells
CY <- c("UBE2C", "TOP2A", "CCNB1", "NUSAP1", "CDKN3", "CCNB2", "BIRC5", "MKI67", "CDC20", "CENPW")
#Cycling epithelial cells(CY)
MD <- c("ANXA1", "KRT13", "GBP6", "LY6D", "KRT4", "CSTA", "RHCG", "S100A8", "TRIM29", "LYPD3")
#Mucosal defense epithelial cells(MD)
TD <- c("SPRR3", "SPRR2E", "SPRR2A", "ECM1", "EMP1", "CNFN", "MAL", "CRNN", "LCN2", "IL1RN")
#Terminal differentiation epithelial cells (TD)
HY <- c("HIF1A", "NDRG1", "NEAT", "HSPA1A", "TMSB10", "PTHLH", "MMP10", "LAMC2", "EGLN3", "SAA1")
#Hypoxia-related stress epithelial cells (HY)
RS <- c("KRT17", "SOD2", "CTAG2", "CCL20", "TMEM45A", "MT1G", "TYMP", "DDIT4", "FABP4", "ANO1", "FADD")
#Reactive oxygen species-related stress epithelial cells (RS)
DO <- c("SPP1", "CES1", "FGFBP2", "GSTM2", "GSTM3", "AKR1C3", "PTGR1", "NQO1", "AKR1C2", "GPX2", "CBR1")
#Deoxidation epithelial cells (DO)
AP <- c("CD74", "HLA-DRA", "HLA-DRB1", "VIM", "HLA-DPB1", "HLA-DQB1", "IFITM1", "HLA-DMA", "IL32", "CFH")
#Antigen presenting epithelial cells (AP)
#cancer-related programs including HY, RS, DO, and AP
Epithelial_list <- list(QP = QP,
                        CY = CY,
                        MD = MD,
                        TD = TD,
                        HY = HY,
                        RS = RS,
                        DO = DO,
                        AP = AP)

