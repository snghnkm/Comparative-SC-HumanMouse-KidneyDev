# Figure_3.R

# Fig 3A.
ggsave(DimPlot(human_Pax8, label=T, repel=T), filename="human_Pax8.pdf")
ggsave(DimPlot(mouse_Pax8, label=T, repel=T), filename="mouse_Pax8.pdf")

# Fig 3B.
ggsave(UMAPPlot(mouse_Pax8.podocyte_PEC.podocyte.rm7111217, label = TRUE, raster=T, cols=colour_ramp(c("#B2A032", "#1ca187"))(seq(0, 1, length = 9))) + NoLegend(), filename="075_mouse_Pax8.podocyte_PEC.podocyte.rm7111217.UMAP.2.pdf") 

Idents(mouse_Pax8.podocyte_PEC.podocyte.rm7111217) <- factor(Idents(mouse_Pax8.podocyte_PEC.podocyte.rm7111217), levels= c(6,7,4,8,1,2,5,3,0))

# Fig 3D.
selected_genes.mouse <- c("Lhx1", "Slc16a1", "Pax8", "Nphs1", "Magi2", "Podxl", "Anxa1", "Col4a4", "Vegfa")
p2 <-plot_pseudotime_heatmap(mouse_Pax8.podocyte_PEC.podocyte.monocle2[selected_genes.mouse,], return_heatmap = T, show_rownames=T, num_clusters=1, cluster_rows = F, cores = 32)
ggsave(p2, filename="mouse_Pax8.podocyte_PEC.podocyte.trajectory.selected.pdf", height=round(selected_genes  %>% length() / 10), width=3, limitsize = F)

selected_genes <- c("Lhx1", "Slc16a1", "Pax8", "Olfm3", "Cdh11", "Pcdh9", "Nphs1", "Magi2", "Podxl", "Anxa1", "Col4a4", "Vegfa", "Pla2r1", "Dcn", "Mafb")
p2 <-plot_pseudotime_heatmap(human_Pax8.podocyte_PEC.podocyte.monocle2[toupper(selected_genes),], return_heatmap = T, show_rownames=T, num_clusters=1, cluster_rows = F, cores = 32)
ggsave(p2, filename="human_Pax8.podocyte_PEC.podocyte.trajectory.selected.pdf", height=round(selected_genes  %>% length() / 10), width=3, limitsize = F)

# Fig 3E.
selected_genes <- c("Lhx1", "Slc16a1", "Pax8", "Olfm3", "Cdh11", "Pcdh9", "Nphs1", "Magi2", "Podxl", "Anxa1", "Col4a4", "Vegfa", "Pla2r1", "Dcn", "Mafb", "mt-Cytb", "Gm42418", "mt-Co2")

ggsave(FeaturePlot(mouse_Pax8.podocyte_PEC.podocyte.rm7111217, feature=c("Olfm3", "Pcdh9", "Pla2r1", "Dcn"),  label = F, min.cutoff="q10", max.cutoff="q90", raster=T, pt.size=1.5), width=8*2, height=6*2, filename="075_mouse_Pax8.podocyte_PEC.podocyte.rm7111217.featureplot.4.2.pdf") 

ggsave(FeaturePlot(human_Pax8.podocyte_PEC, features = toupper(selected_genes), label = TRUE, repel = TRUE, min.cutoff = "q10", max.cutoff = "q90"), filename = "002_human_Pax8.podocyte_PEC.featureplots.pdf", width = 7*4, height = 7*3)

selected_genes_human_specific <- c("Olfm3", "Pcdh9", "Cdh11", "Pla2r1", "Dcn")
ggsave(FeaturePlot(mouse_Pax8.podocyte_PEC.podocyte, features=selected_genes_human_specific, label=T, repel=T, order=T, min.cutoff = "q10", max.cutoff = "q90", ncol=3), filename="004_2_mouse_Pax8.podocyte_PEC.podocyte.featureplots.pdf", width=7*3, height=7*2)

# Fig 3F-I.
ggsave(UMAPPlot(mouse_Pax8.proximal_tubule, label = TRUE, cols=colour_ramp(c("#a1cfc0", "#73b9d2"))(seq(0, 1, length = 15)), raster=T) + NoLegend(), filename="001_mouse_Pax8.proximal_tubule.res.1.2.2.pdf", width=10)

ggsave(FeaturePlot(mouse_Pax8.proximal_tubule, features=c("Wnt4", "Hnf4a", "Osr2", "Wnt16", "Slc39a5", "Dab2", "Ly6a","Hey1", "Sox11", "Lef1", "Aldob", "Slc34a1", "Spp2", "Tmem252", "Sat1", "C1s1", "Bmp3", "Slc5a12", "Gatm", "Kap", "Slc17a1", "Lhx1", "Ifit3", "Slc17a1", "Umod", "Cldn6", "Serpina6", "Slc7a12", "Slc7a13", "Pitx2", "Aadat", "Slc22a13", "Aqp1", "Cdh6", "Epha7", "Slc5a2", "Cyp2e1", "Slc22a2","Slc17a3", "Cyp2j5", "Acox2", "Slc13a3", "Slc5a11", "Gdf15", "Angpt2", "Satb2", "Corin", "Hnf4g", "Akap12", "Cldn1", "Pou3f3", "Irx1", "Umod", "Cdh4", "Maf", "Hnf1b", "Hnf1a", "Hnf4g", "Thrb"), label=T, repel=T, order=T, min.cutoff = "q10", max.cutoff = "q90"), filename="002_mouse_Pax8.proximal_tubule.featureplots.pdf", width=7*4, height=7*8, limitsize = F)

ggsave(FeaturePlot(human_Pax8.proximal_tubule, features=toupper(c("Wnt4", "Hnf4a", "Osr2", "Wnt16", "Slc39a5", "Dab2", "Ly6a","Hey1", "Sox11", "Lef1", "Aldob", "Slc34a1", "Spp2", "Tmem252", "Sat1", "C1s1", "Bmp3", "Slc5a12", "Gatm", "Kap", "Slc17a1", "Lhx1", "Ifit3", "Slc17a1", "Umod", "Cldn6", "Serpina6", "Slc7a12", "Slc7a13", "Pitx2", "Aadat", "Slc22a13", "Aqp1", "Cdh6", "Epha7", "Slc5a2", "Cyp2e1", "Slc22a2","Slc17a3", "Cyp2j5", "Acox2", "Slc13a3", "Slc5a11", "Gdf15", "Angpt2", "Satb2", "Corin", "Hnf4g", "Akap12", "Cldn1", "Pou3f3", "Irx1", "Umod", "Cdh4")), label=T, repel=T, order=T, min.cutoff = "q10", max.cutoff = "q90"), filename="002_human_Pax8.proximal_tubule.featureplots.pdf", width=7*4, height=7*8, limitsize = F)


# Fig 3J-M.
ggsave(UMAPPlot(mouse_Pax8.distal_nephron, label = F, cols=colour_ramp(c("#73b9d2", "#465584"))(seq(0, 1, length = 17))) + NoLegend(), filename="075_mouse_Pax8.distal_nephron.UMAP.2.tiff", width=10) 

ggsave(UMAPPlot(mouse_Pax8.distal_nephron, label = TRUE, raster=T, cols=colour_ramp(c("#73b9d2", "#465584"))(seq(0, 1, length = 17))) + NoLegend(), filename="075_mouse_Pax8.distal_nephron.UMAP.2.pdf", width=10) 
ggsave(UMAPPlot(mouse_Pax8.distal_nephron, label = F, cols=colour_ramp(c("#73b9d2", "#465584"))(seq(0, 1, length = 17))) + NoLegend(), filename="075_mouse_Pax8.distal_nephron.UMAP.2.tiff", width=10) 

ggsave(FeaturePlot(mouse_Pax8.distal_nephron, features=c("Ptger3", "Cldn10", "Slc12a1", "Foxq1", "Cldn16", "Pappa2", "Nos1", "Slc12a3", "Calb1", "Trpv5", "Bst1", "Wnt4", "Fgf8", "Bmp2", "Lhx1", "Tfap2a", "Tfap2b", "Tmem52b", "Sall3", "Sox9", "Mecom", "Emx1", "Gata3", "Rprm", "Ifit3", "Satb2", "Pou3f3", "Lgr5", "Slc6a12", "Bcl6", "Sim2", "Umod", "Irx1", "Esrrb", "Esrrg", "Zfp503", "Wnt7b", "Scnn1b", "Aqp2", "Trpm6", "Gdf15", "Clcnka", "Aqp1", "Slc14a2", "Pitx2", "Corin"), label=T, repel=T, order=T, min.cutoff = "q10", max.cutoff = "q90"), filename="002_mouse_Pax8.distal_nephron.featureplots.pdf", width=7*4, height=7*5)
