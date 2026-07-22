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
  Pericyte=Pericyte
)


all_strings <- unlist(esophagus_markers_list)
# 2. 统计每个字符串的出现次数
string_counts <- table(all_strings)
# 3. 找出重复的字符串（出现次数大于1）
duplicates <- names(string_counts)[string_counts > 1]


####Tcell####
t1 <- list(cd4 = c("CD4", "IL7R", "CCR7", "SELL"),
           cd8 = c("CD8A","CD8B","GZMB","PRF1","GNLY"),
           gamaT = c("TRDC", "TRDV1","TRDV2","TRGV9"), #Gamma Delta T cells
           NKT = c("ZBTB16","CD56","KLRB1"),  #Natural Killer T cells
           NKT2 =  c("EOMES", "XCL1", "XCL2", "CXCR6", "TIGIT", "LAG3"),
           MAIT = c("TRAV1-2", "SLC4A10", "RORC", "RORA"), #Mucosal-Associated Invariant T cells
           IEL = c("CD160", "KIR2DL4", "TMIGD2", "ITGAE"), #Intraepithelial Lymphocytes
           DNT = c("GZMK")) #Double Negative T cells
t2 <- list(cd4 = c("CD4"),
           cd8 = c("CD8A","CD8B"))

t3 <- list(naive = c("CCR7","SELL","IL7R","TCF7","LEF1","LTB","FOXO1","GZMK"),
           acti = c("SELL","CD40LG","ANXA1","IL2RA","CD69"), #Activation
           memo = c("CCR7", "TCF7", "CD69", "NR4A1", "MYADM", "GATA3", "TBX21"), #Memory
           effe = c("BATF3","IL2RA","IFNG","IRF4","MYC","SLC7A5","SLC7A1","XCL1","GZMB","CCL3","CCL4","IL2","PRF1","NKG7","GNLY"), #Effector
           tole = c("EGR2","IZUMO1R","CD200","DGKZ","BTLA","CTLA4"), #Tolerant
           exha = c("PDCD1","LAG3","HAVCR2","CD244","CD160","IL10R","EOMES","NR4A2","PTGER4","TOX","TOX2","TIGIT","CTLA4","ENTPD1"), #Exhausted
           pro = c("MKI67","TK1","STMN1"), #Proliferation
           toxi = c("GZMK","GZMA","GZMB","NKG7","PRF1","IFNG","GNLY"), #Cytotoxicity
           regu = c("FOXP3","IL2RA","CTLA4","TIGIT","TNFRSF18","IL10","IKZF2","CCR8"), #Regulatory
           earlym =  c("CCR2","CX3CR1","IL1RRAP","ZEB2")) #Early Memory

t4 <- list(naive = c("CCR7","SELL","IL7R","TCF7","LEF1","LTB","FOXO1","GZMK"),
           memo = c("CCR7", "TCF7", "CD69", "NR4A1", "MYADM", "GATA3", "TBX21"),
           effe = c("BATF3","IL2RA","IFNG","IRF4","MYC","SLC7A5","SLC7A1","XCL1","GZMB","CCL3","CCL4","IL2","PRF1","NKG7","GNLY"),
           exha = c("PDCD1","LAG3","HAVCR2","CD244","CD160","IL10R","EOMES","NR4A2","PTGER4","TOX","TOX2","TIGIT","CTLA4","ENTPD1"),
           regu = c("FOXP3","IL2RA","CTLA4","TIGIT","TNFRSF18","IL10","IKZF2","CCR8"))

#Progenitor-like exhausted SPRY1+CD8+ T cells potentiate responsiveness to neoadjuvant PD-1 blockade in esophageal squamous cell carcinoma
t5 <- list(Tem = c("GNLY","FOS","JUN","IL7R"),
           Tex = c("PDCD1", "HAVCR2", "LAG3","ENTPD1", "CXCL13","HLA-DR"),
           Trm = c("ZNF683"))
#Single-cell multi-stage spatial evolutional map of esophageal carcinogenesis
#CD8
t6 <- list(NM = c("TCF7", "CCR7","SELL"), #naı¨ve/memory
           E = c("GZMB", "GZMA", "GZMH", "NKG7", "GNLY"), #effector T 
           EX = c("PDCD1", "CXCL13", "CTLA4", "HAVCR2")) #exhausted T cells 
#CD4
t7 <- list(NM = c("TCF7", "CCR7","SELL"), #naı¨ve/memory
           TFH = c("CXCL13", "CD200"), #T follicular helpers
           Tre = c("FOXP3", "IL2RA", "CTLA4"))  
#以下2个为综述Unraveling the tumor microenvironment of esophageal squamous cell carcinoma through single-cell sequencing: A comprehensive review
#CD4
t8 <- list(navie = c("CCR7","SELL", "TCF7", "LEF1", "IL7R"),
            Memory = c("LMNA","IL7R"),
            Exhausted = c("CXCL13","CARD16"),
            Treg = c("FOXP3","IL2RA","IKZF2","CD4","CD25"),
            Helper17 = c("KLRB1","CCL20","RORA"),
            Follicular_helper = c("CXCR5","CXCL13","MAF"))
#CD8
t9 <- list(navie = c("TCF","SELL","LEF1","CCR7","CD45RA"),
              Cytotoxic = c("GZMK","GZMN","NKG7","KLRG1","GZMA","GNLY","PRF1","GZMB","IFNG"),
              Memory = c("LMNA","CCR7"),
              Exhausted = c("ENTPD1","PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","SPRY"))
tcell_marker <- list(t1 = t1, t2 = t2, t3 = t3,t4 = t4, t5 = t5, t6 = t6,
                     t7 = t7, t8 = t8, t9 = t9)
####Bcell####
#综述Unraveling the tumor microenvironment of esophageal squamous cell carcinoma through single-cell sequencing: A comprehensive review
#CD4
b1 <- list(Follicular = c("MS4A1"),
                   Naive = c("IGHD", "TCL1A", "FCER2", "CD23"), 
                   memory = c("FCRL4", "NR4A1","TNFRSF13B","TACI"),
                   Active = c("GPR183","TNFRSF13B","CD69","CD83","HLA-DRA"),
                   Proliferative = c("THBA1B"),
                   BGCs = c("MEF2B","BCL6","NEIL1"))
bcell_marker <- list(b1 = b1)
#BGCs来自 PMID:38734159


#,Plasma = c("IGHG1","JCHAIN")
####Myeloid####
#综述Unraveling the tumor microenvironment of esophageal squamous cell carcinoma through single-cell sequencing: A comprehensive review
#CD4
myeloid_list <- list(Monocyte = c("LYZ","FCN1","CCL3"), #"LYZ"在巨噬细胞中极高表达。
                     Macrophage	= c("C1QA","CCL18"),
                     Neutrophil = c("FCGR3B","S100A8","CSF3R"),
                     Dendritic = c("CCL17","CCR7","CCL22"))
macrophage_list <- list(m1 = c("CD86","CD80","IL1B", "IL1A", "TNF", "IL10"),
                        m2 = c("CD163","C1QA","CD68","CD206","MRC1"))
#"MRC1"就是CD206
Dendritic_list <- list(Activated = c("CD40", "FSCN1", "CCR7"),
                  cDC = c("CD1A","C1DC","CLEC10A"), #Conventional DC(cDC)
                  pDC = c("LILRA4","CD123"),
                  tDC = c("IDO1", "FSCN", "LAMP3"))  #Plasmacytoid DC(pDC) Tolerogenic DC(tDC)
cDC_list <- list(cDC1 = c("CLEC9A", "XCR1", "IRF8", "BATF3"),
                 cDC2 = c("FCER1A", "CD1C"))


#A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells
#CD68，CD86,是自己加的。
macro <- c("CD68","CD86","INHBA","IL1RN","CCL4", #Macro_INHBA
           "NLRP3","EREG","IL1B", #Macro_NLRP3
           "LYVE1", "PLTP", "SPP1",# #Macro_LYVE1
           "C1QC", "C1QA", "APOE"#Macro_C1QC
)

mono <- c("FCN1", "S100A9", "S100A8", #Mono_CD14
          "FCOR3A", "LST1", "LILRB2" #Mono_CD16
)

DC <- c("LILRA4", "GZMB", "IL3RA", #pDC_LILRA4
        "CLEC9A", "FLT3", "IDO1", #cDC1_CLEC9A
        "CD1C","FICERIA", "HLA-DOA1", #cDC2_CD1C
        "LAMP3", "CCRY", "FSCN1") #cDC3_LAMP3
myeloid1 <- list(macro = macro,mono = mono, DC = DC)


#Integrative analysis of bulk and single-cell gene expression profiles to identify tumor-associated macrophage-derived CCL18 as a therapeutic target of esophageal squamous cell carcinoma
####好像没什么用
mono <- c("BCL2A1", "C15orf48", "OLR1", "IL10", "IL1A", "CCL3", "CXCL8")
macro <- c("C1QC", "CCL18", "FCGR3A", "MMP9", "MSR1", "PLA2G7")
DC <- c("ENTPD1", "CD86", "CCR7", "SERPINB9", "LAMP3", "CCL17", "CCL22")
myeloid2 <- list(macro = macro,mono = mono, DC = DC)

myeloid3 <- list(
  # 单核细胞 (分为亚群)
  Classical_Monocyte = c("CD14", "S100A8", "S100A9", "VCAN"),
  Nonclassical_Monocyte = c("FCGR3A", "MS4A7", "CX3CR1"),
  
  # 巨噬细胞
  Macrophage = c("CD68", "CD163", "MS4A4A", "C1QA", "APOE", "MARCO"), # 注：巨噬细胞标记具有组织特异性
  
  # 树突状细胞 (必须分为亚群)
  cDC1 = c("CLEC9A", "XCR1", "BATF3"),
  cDC2 = c("CD1C", "CLEC10A", "CD1E"),
  pDC = c("IL3RA", "GZMB", "TCF4", "IRF8"),
  
  # 通用或状态标记 (可选)
  #Pan_Myeloid = c("LYZ", "CST3"), # 通用髓系标记
  mregDC = c("LAMP3", "CCR7", "FSCN1") # 成熟/活化DC状态标记
)
#Mature dendritic cells enriched in immunoregulatory molecules (mregDCs): A novel population in the tumour microenvironment and immunotherapy target
#cell parper
#A conserved dendritic-cell regulatory program limits antitumour immunity
#mregDC（mature dendritic cells enriched in immunoregulatory molecules，富含免疫调节分子的成熟树突状细胞）
#DC3 就是 mregDC（富含免疫调节分子的成熟树突状细胞） 的另一种叫法
myeloid_cell <- list(myeloid_list =myeloid_list,
                     macrophage_list =macrophage_list,
                     Dendritic_list = Dendritic_list,
                     cDC_list = cDC_list,
                     myeloid1 = myeloid1,
                     myeloid2 = myeloid2,
                     myeloid3 = myeloid3)           

####Fibroblast#####
#综述Unraveling the tumor microenvironment of esophageal squamous cell carcinoma through single-cell sequencing: A comprehensive review
#CD4
Fibroblast_list <- list(iCAF = c("CFD","CXCL12","DCN"),
                     myoCAF = c("ACTA2","TAGLN"))
#PMID:40068596
Fibroblast_list2 <- list(normal = c("IGFBP6","PI16","DPT","CXCL12","GPX3", "DCN"),
                      cancer = c("FAP","MMP1","MMP11","POSTN"))
fibroblast_cell <- list(fi1 = Fibroblast_list,
                     fi2 = Fibroblast_list2)
####Epithelial cells####
#综述Unraveling the tumor microenvironment of esophageal squamous cell carcinoma through single-cell sequencing: A comprehensive review
#CD4
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
#PMID:40068596
epi2 <- list(basal= c("COL17A1","NOTCH1"),
             proliferative = c("MKI67","JAG1"),
             differentiated = c("ANXA1","ECM1"),
             invasive = c("SOX2","NFE2L2","SPP1","MMP11"))

epi3 <- list(normal = c("LIPF","GKN1", "PGC"), #不太可靠？
             cancer = c("CLDN7","TFF3","CLDN4"),
             cancer1 = c("EPCAM", "TSTA3", "NECTIN4", "SFN", 
                         "KRT5", "KRT16", "KRT18", "KRT17", 
                         "KRT19", "KRT6A", "KRT6B")
             )

epithelial_cell <- list(epi1 = Epithelial_list,
                        epi2 = epi2)
#####endothelial ###########
#PMID:40068596
endo1 <- list(normal = c("NECs","ACKR1","CXCL12","POSTN","ENG","VWF"),
              tumor = c("TECs","PLVAP","SERPINE1", "SPP1"))
endothelial_cell <- list(en1 = endo1)
