Tcell=c("Cd3d", "Cd3e", "Cd3g","Cd2")
Bcell=c("Cd79a", "Ms4a1", "Cd19")
Myeloid <- c("Cd68", "Lyz2", "Cd14", "Cd163", "Il3ra", "Lamp3", "Clec4c", "Tpsab1")
Muscle <- c("Pln", "Myh11", "Actg2", "Cnn1", "Myl9", "Tagln")
Epithelial <- c("Epcam", "Sfn", "Krt5", "Krt14")
Endothelial=c("Pecam1", "Cdh5", "Vwf", "Kdr", "Esam")
Fibroblast=c("Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra", "Pdgfrb", "Acta2")
Glandular <- c("Agr2", "Krt23", "Wfdc2")
NK=c("Nkg7", "Klrd1", "Klrk1", "Prf1", "Gzmb")
Plasma=c("Jchain", "Sdc1", "Xbp1", "Ighm", "Igha", "Mzb1")
Mast <- c("Cpa3", "Tpsab1", "Kit")
Pericyte=c("Rgs5", "Des", "Cspg4", "Pdgfrb")
Immune <- c("Ptprc")
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
  immune=Immune
)

z_width <- length(levels(Idents(sce.all.int)))
genes_to_check = lapply(get(x), str_to_title)
dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
genes_to_check = lapply(genes_to_check, function(x) x[!x %in% dup])
DotPlot(sce.all.int , features = genes_to_check )  + 
 # coord_flip() + 
 theme(axis.text.x=element_text(angle=45,hjust = 1))
w=length( unique(unlist(genes_to_check)) )/5+6;w
ggsave(paste('marker.pdf'), height =  z_width / 2 ,width  = w + 4)
 
  