# Figure 2.

# Fig 2B.
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223 <- readRDS("Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Rds")

DefaultAssay(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223) <- 'integrated'
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223 <- FindClusters(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, resolution = 1.0)
DefaultAssay(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223) <- "RNA"
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223 <- NormalizeData(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223)

saveRDS(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, "015_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.res1.0.RNA_norm.Rds", compress=F) # final Mouse N-UE object

Idents(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223) <- factor(Idents(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223), levels= c(
    0, 2, 18, 13, 19, 11, 1, 16, 9, 20, 6, 23, 4, 8, 3, 10, 21, 25, 7, 5, 15, 12, 22, 17, 14, 26, 24))


Mp0_NUE_color <- unique(c(
colour_ramp(c("#cf5742", "#B2A032"))(seq(0, 1, length = 2+7)),
colour_ramp(c("#B2A032", "#1ca187"))(seq(0, 1, length = 2+3)),
colour_ramp(c("#1ca187", "#a1cfc0"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#a1cfc0", "#73b9d2"))(seq(0, 1, length = 2+2)),
colour_ramp(c("#73b9d2", "#8590b0"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#8590b0", "#465584"))(seq(0, 1, length = 2+3)),
colour_ramp(c("#465584", "#7C5B7A"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#7C5B7A", "#958abc"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#958abc", "#6a2449"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#6a2449", "#bd3760"))(seq(0, 1, length = 2+0)),
colour_ramp(c("#bd3760", "#dd88a1"))(seq(0, 1, length = 2+0))
))
length(Mp0_NUE_color)

p <- DotPlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features=rev(c(
  "Cited1", "Six2", "Bmper", "Wnt4", "Wt1", "Mafb",  "Nphs2", "Tagln", "Osr2", "Hnf4a", "Lrp2", "Slc34a1", "Fgf8", "Tfap2b", "Slc12a1", "Irx1", "Umod", "Calb1", "Ret", "Rprm", "Aqp2","Upk1b", "Foxi1", "Atp6v1g3" )), cols = c("lightgrey", "red"), col.min = 0) + coord_flip()

ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features=c("Wnt4","Eya1", "Ret", "Scnn1g"), label = F, min.cutoff="q10", max.cutoff="q90", order=T, pt.size=1.5), filename="016_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.UMAP.2.featureplot.2.pdf", width=9*2, height=7*2) 
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features=c("Hoxd10", "Hoxd11"), label = F, min.cutoff="q10", max.cutoff="q90", order=T, pt.size=1.5), filename="016_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.UMAP.2.featureplot.Hoxd10.Hoxd11.pdf", width=9*2, height=7) 


ggsave(p, filename="016_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.dotplot.markers.pdf", width=12, height=6)
ggsave(UMAPPlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, label=T, repel=T, raster=T, pt.size=0.1, cols=Mp0_NUE_color) + NoLegend(), filename="011_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.UMAP.res1.0.pdf")
ggsave(UMAPPlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, label=F, repel=F, raster=10000000, pt.size=0.1, cols=Mp0_NUE_color) + NoLegend(), filename="011_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.UMAP.res1.0.tiff")
ggsave(UMAPPlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, group.by = "replicate", raster=T, pt.size=0.1) + NoLegend(), filename="017_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.by.replicate.pdf")

t <- Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223@meta.data %>% dplyr::select(c(integrated_snn_res.1, replicate))
t$integrated_snn_res.1 <- as.character(t$integrated_snn_res.1)
t$integrated_snn_res.1 <- factor(t$integrated_snn_res.1, levels=c(0, 2, 18, 13, 19, 11, 1, 9, 16, 20, 6, 23, 4, 8, 3, 7, 10, 21, 5, 15, 12, 25,  22, 17, 14, 26, 24))

ggsave(
    ggplot(t, aes(x=integrated_snn_res.1, fill=replicate)) + geom_bar(position = "fill") + ggpubr::theme_pubr(),
    filename="018_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.geom_bar.replicate.pdf",
    width=12,
    height=6,
    limitsize=F
)

as.data.frame(table(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223@active.ident)) %>% readr::write_excel_csv("019_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.cellcountsbycluster.csv")

ggsave(
    Seurat::VlnPlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, c(
        "nCount_RNA", "nFeature_RNA", "percent.mt"),
    pt.size=0, ncol=1),
    filename="020_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.QC.vlnplot.pdf",
    width=12,
    height=3*3,
    limitsize = F
)

ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Nphs2", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Nphs2.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Slc34a1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Slc34a1.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Slc12a1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Slc12a1.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Aqp2", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Aqp2.featureplot.pdf")


ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Akap12", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Akap12.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Cldn1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Cldn1.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Wnt16", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Wnt16.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Umod", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Umod.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Foxi1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Foxi1.featureplot.pdf")
ggsave(FeaturePlot(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, features="Upk1b", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 4, order=T), filename="Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Upk1b.featureplot.pdf")


# Fig 2C
DefaultAssay(hm.integrated.rm8) <- "novIntPeaks"

# re-compute the UMAP using corrected LSI embeddings
hm.integrated.rm8 <- RunUMAP(hm.integrated.rm8, dims = 2:30, reduction = 'harmony')

p5 <- DimPlot(hm.integrated.rm8, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
ggsave(p5, filename="090_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.UMAP.pdf")

hm.integrated.rm8 <- FindNeighbors(hm.integrated.rm8, reduction = 'harmony', dims = 2:30)


hm.integrated.rm8 <- hm.integrated.rm8 %>%
    FindClusters(verbose = T, algorithm = 3, resolution = 2.0)
ggsave(UMAPPlot(hm.integrated.rm8, label=T, repel=T), filename="091_1_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.res2.0.UMAP.pdf")

hm.integrated.rm8 <- hm.integrated.rm8 %>%
    FindClusters(verbose = T, algorithm = 3, resolution = 1.6)
ggsave(UMAPPlot(hm.integrated.rm8, label=T, repel=T, raster=T)+NoLegend(), filename="091_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.res1.6.UMAP.pdf")
ggsave(UMAPPlot(hm.integrated.rm8, label=F, repel=F, raster=F, cols=Hufetal_NUE_snATAC_color)+NoLegend(), filename="091_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.res1.6.UMAP.tiff")

Hufetal_NUE_snATAC_color <- unique(c(
    colour_ramp(c("#cf5742", "#B2A032"))(seq(0, 1, length = 2+3)),
    colour_ramp(c("#B2A032", "#1ca187"))(seq(0, 1, length = 2+1)),
    colour_ramp(c("#1ca187", "#a1cfc0"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#a1cfc0", "#73b9d2"))(seq(0, 1, length = 2+3)),
    colour_ramp(c("#73b9d2", "#8590b0"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#8590b0", "#465584"))(seq(0, 1, length = 2+5)),
    colour_ramp(c("#465584", "#7C5B7A"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#7C5B7A", "#958abc"))(seq(0, 1, length = 2+1)),
    colour_ramp(c("#958abc", "#bd3760"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#bd3760", "#dd88a1"))(seq(0, 1, length = 2+0))
))
length(Hufetal_NUE_snATAC_color)

# Fig 2C

Hufetal_NUE_color <- unique(c(
    colour_ramp(c("#cf5742", "#B2A032"))(seq(0, 1, length = 2+3)),
    colour_ramp(c("#B2A032", "#1ca187"))(seq(0, 1, length = 2+2)),
    colour_ramp(c("#1ca187", "#a1cfc0"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#a1cfc0", "#73b9d2"))(seq(0, 1, length = 2+2)),
    colour_ramp(c("#73b9d2", "#8590b0"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#8590b0", "#465584"))(seq(0, 1, length = 2+3)),
    colour_ramp(c("#465584", "#7C5B7A"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#7C5B7A", "#958abc"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#958abc", "#6a2449"))(seq(0, 1, length = 2+1)),
    colour_ramp(c("#6a2449", "#bd3760"))(seq(0, 1, length = 2+0)),
    colour_ramp(c("#bd3760", "#dd88a1"))(seq(0, 1, length = 2+0))
))
length(Hufetal_NUE_color)


ggsave(UMAPPlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, label = TRUE, raster=T) + NoLegend(), filename="075_integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.UMAP.2.pdf") 
ggsave(UMAPPlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, label = F, cols=Hufetal_NUE_color) + NoLegend(), filename="075_integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.UMAP.2.tiff") 


# Fig 2D
t <- Hp0_10_SCT.combined.rm022.rm3134.rm303132@meta.data %>% dplyr::select(c(seurat_clusters, orig.ident, replicate, sex))
t$seurat_clusters <- as.character(t$seurat_clusters)

t <- t %>% mutate(stage = case_when(
  replicate %in% c("BGI_M1", "BGI_M2", "BGI_M3", "BGI_F1") ~ "1_Early",
  replicate %in% c("BGI_F2", "BGI_F3") ~ "2_Middle",
  replicate %in% c("NOV_M4", "NOV_M5", "NOV_F5", "NOV_F6") ~ "3_Late"
))

t$seurat_clusters <- factor(t$seurat_clusters, levels=c(4,2,14,22,16,21,17,0,8,12,7,20,23,9,24,5,25,6,1,26,30,13,31,29,18,10,3,28,19,15,11,27))

ggsave(
    ggplot(t, aes(x=seurat_clusters, fill=stage)) + geom_bar(position = "fill"),
    filename="080_Hp0_10_SCT.combined.rm022.rm3134.rm303132.stage.bar.pdf", width=12, height=6, limitsize=F
)

Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- AddMetaData(
  object = Hp0_10_SCT.combined.rm022.rm3134.rm303132,
  metadata = t$stage,
  col.name = 'stage'
)

ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, split.by = "stage", label = FALSE) + NoLegend(), filename="076_integrated.N.RNA_norm.UMAP.stage.pdf", width=21, height=7)

table(t$seurat_clusters, t$stage)
prop.table(table(t$seurat_clusters, t$stage))

ggsave(
    ggplot(t, aes(x=seurat_clusters, fill=replicate)) + geom_bar(position = "fill") + ggpubr::theme_pubr(),
    filename="080_integrated.N.RNA_norm.replicate.bar.pdf", width=12, height=6, limitsize=F
)

ggsave(ggpubr::ggbarplot(prop.table(table(t$seurat_clusters, t$stage),2) %>% as.data.frame(), "Var1", "Freq", fill="Var2", label=FALSE, position = position_dodge(0.9)),
  filename="080_integrated.N.RNA_norm.replicate.bar.by.cluster.pdf", width=12, height=4)

# Fig 2E

dotplot_features = c("EYA1", "BMPER", "TOP2A", "WNT4", "LHX1", "MAFB", "NPHS2", "AKAP12", "LRP2", "SLC34A1", "SATB2", "TFAP2B", "SLC12A1", "UMOD", "TMEM52B", "RET", "AQP2", "UPK3A", "FOXI1")
p <- DotPlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features=rev(dotplot_features), col.min=0) + coord_flip()
ggsave(p, filename="010_integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.dotplot.markers.eps", width=7.7, height=7.0)

# Fig 2G

ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="NPHS2", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.NPHS2.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="SLC34A1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.SLC34A1.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="SLC12A1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.SLC12A1.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="AQP2", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.AQP2.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="CALB1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.CALB1.featureplot.pdf", width=15, height=15)

ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="AKAP12", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.AKAP12.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="CLDN1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.CLDN1.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="UMOD", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.UMOD.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="FOXI1", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.FOXI1.featureplot.pdf", width=15, height=15)
ggsave(FeaturePlot(integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm, features="UPK1A", min.cutoff="q10", max.cutoff="q90", raster=F, pt.size = 1.5, order=T), filename="integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.UPK1A.featureplot.pdf", width=15, height=15)


# Fig 2M
hm.integrated.rm8.GA.averages <- AverageExpression(hm.integrated.rm8, assays =  "RNA", return.seurat = TRUE)
p <- DoHeatmap(hm.integrated.rm8.GA.averages, features = rownames(hm.integrated.rm8.GA.averages[['RNA']]), size = 3, draw.lines = F)

mat <- p$data %>% dplyr::select(Feature, Cell, Expression) %>% pivot_wider(names_from=Feature, values_from=Expression) %>% tibble::column_to_rownames("Cell")

hm.integrated.rm8.GeneActivity.AllMarkers <- readRDS("095_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.GeneActivity.AllMarkers.Rds")
GA_top10_per_cluster <- hm.integrated.rm8.GeneActivity.AllMarkers %>% group_by(cluster) %>% top_n(wt=-p_val_adj, n=3) %>% top_n(wt=avg_logFC, n=3*7) %>% ungroup() %>% select(c(`gene`)) %>% unique()
GA_top10_per_cluster # total 206 GA

t_mat <- t(mat) %>% as.data.frame() %>% rownames_to_column(var="gene")

t_mat3 <- t_mat %>% filter(gene %in% pull(GA_top10_per_cluster, gene))

t_mat3 %>% write_excel_csv("014_hm.integrated.rm8.GA.averages.pheatmap.csv")

GA_order <- c(c("HMCN1", "CAPZB", "NKAIN3", "FAT3", "BMPER", "LAMA4", "EYA1", "DCHS2", "GABBR2", "ITGA8", "COLGALT2", "GFRA1", "ATP2B2", "KIF26B", "ELAVL2", "DCC", "AUTS2", "HOXC4", "CDH11", "AFF3", "HAS2", "TRPM5", "ALK", "SH3PXD2A", "SLIT1", "HOXC6", "XXYLT1", "KIAA2012.1", "RP11-834C11.12", "SORCS2", "KIAA2012", "LYPD1", "CNTFR", "SH3RF3", "DAPL1", "CUX2", "COL13A1", "LTBP1", "MEIS1", "SLC6A5", "BASP1", "TMEM132D", "SLIT3", "WBSCR17", "GLI2", "WT1", "CCBE1", "UNC5C", "CD6", "CRYGN", "OMA1", "MGAT5B", "FBLN2", "OLFM3", "PCDH9", "NPHS2", "PODXL", "PTPRO", "MCC", "SYNPO", "NPHS1", "ARHGAP28", "GRK5", "KLHL29", "ZBTB7C", "FRY", "CBLB", "CLIC5", "IQCJ", "CDC14A", "TJP1", "TGFBR3", "PLCE1", "IQGAP2", "MAGI2", "THSD7A", "PLA2R1", "ST3GAL6", "MME", "CRIM1", "ITGB3", "MEGF11", "PITPNC1", "FUT9", "CDH6", "PDZK1IP1", "HNF4G", "TTBK1", "SLC27A2", "C1QTNF3", "CHRNA4", "CTD-2501B8.1", "SLC6A19", "ACMSD", "ADAMTS12", "HNF4A", "SMIM24", "GLYATL1", "SLC5A8", "ANPEP", "ALPL", "DAB2", "PCK1", "LRP2", "SLC7A7", "C9", "PRLR", "SLC34A3", "SLC34A1", "RAB11FIP3", "NIT2", "TNFSF10", "CUBN", "SLC13A3", "AK4", "ALDOB", "GLYAT", "B3GNT3", "SATB2", "IRX1", "ARFGEF3", "CASR", "GRID1", "CDH4", "DAB1", "SLC39A8", "DNAH5", "PAPPA2", "NTRK3", "ERBB4", "ESRRG", "STK32B", "CLCN5", "ACPP", "GABRB3", "SLC12A1", "LNX1", "AGPAT4", "KCNJ1", "BTBD11", "TBC1D4", "PTGER3", "TPST2", "THSD4", "SIM2", "TFCP2L1", "PRKCQ", "SLC4A4", "MECOM", "COBLL1", "LIMCH1", "IQCK", "ADGRL2", "EHF", "WDR72", "GRHL2", "PRDM16", "SALL3", "TMEM52B", "FAR2", "PAH", "TBC1D9", "KCTD8", "CCDC192", "GALNT13", "GATA3", "APCDD1L", "AHI1", "SCNN1B", "ALDH3B2", "SCNN1G", "PCDH7", "PZP", "SDR42E1", "CAMK4", "NCALD", "TWIST1", "KCNN4", "MYEOV", "COL5A2", "KANK4", "LGI1", "MMP2", "BCAT1", "WNT9B", "ST8SIA6", "CSMD3", "KALRN", "TBX3", "SLC9A4", "EYA4", "SYT8", "TMEM45B", "PPL", "SCIN", "ELF5", "WNT7B", "NECTIN4", "PDE7B", "AIM1", "LINC00675", "PPARG", "RNF223", "KRT15", "IL1RN", "S100P"))

p2 <-pheatmap::pheatmap((t_mat3 %>% column_to_rownames("gene"))[GA_order,], clustering_method="ward.D2", display_numbers = F, fontsize_row = 10, cluster_cols=F, cluster_rows=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100))
ggsave(p2, filename="014_hm.integrated.rm8.GA.averages.pheatmap.final.eps", width=30, height=40, limitsize = FALSE)

# Fig 2N
hm.integrated.rm8.differential.motifs <- FindAllMarkers(object = hm.integrated.rm8,  assay="novIntPeaks_JASPAR_chromVAR",
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_common_ATAC' #important?
)
saveRDS(hm.integrated.rm8.differential.motifs, "hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.Rds", compress=F)
hm.integrated.rm8.differential.motifs %>% write_excel_csv("hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.csv")
hm.integrated.rm8.differential.motifs %>% dplyr::select(cluster, gene) %>% group_by(cluster) %>% mutate(id=row_number()) %>% pivot_wider(names_from = cluster, values_from = gene) %>% as_tibble() %>% write_excel_csv("hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.pivot_wider.csv")


DefaultAssay(hm.integrated.rm8) <- "novIntPeaks_JASPAR_chromVAR"
hm.integrated.rm8.JASPAR.averages <- AverageExpression(hm.integrated.rm8, assays =  "novIntPeaks_JASPAR_chromVAR", return.seurat = TRUE)


pfm.JASPAR2020 <- readRDS("new_signac_Mp0/JASPAR2020.9606.10090.pfm.Rds")


p <- DoHeatmap(hm.integrated.rm8.JASPAR.averages, features = unique(ID(pfm.JASPAR2020)), size = 3, draw.lines = F)

mat <- p$data %>% dplyr::select(Feature, Cell, Expression) %>% pivot_wider(names_from=Feature, values_from=Expression) %>% tibble::column_to_rownames("Cell")

pfm.ID.name <- tibble(id=ID(pfm.JASPAR2020), name=name(pfm.JASPAR2020))

t_mat <- left_join(t(mat) %>% as.data.frame()%>% rownames_to_column(var="motif"), pfm.ID.name, by=c("motif" = "id"))
t_mat <- t_mat %>% unite("id_name", c("motif","name"))
t_mat <- t_mat %>% unique()
rownames(t_mat) <- t_mat$id_name
t_mat %>% rownames_to_column("motif") %>% write_excel_csv("hm.integrated.rm8.differential.motifs.JASPAR_chromVAR.averages.pheatmap.csv")
t_mat <- t_mat %>% dplyr::select(!c(id_name))
t_mat2 <- t_mat %>% mutate(maxma = do.call(pmax, .[-1])) %>% dplyr::filter(maxma >2) %>% dplyr::select(!c("maxma"))
t_mat2 %>% dim()
# t_mat2 %>% arrange()
# t_mat2 %>% map_df(sort, decreasing = TRUE)

p2 <-pheatmap::pheatmap(t_mat, clustering_method="ward.D2", display_numbers = F, fontsize_row = 10, cluster_cols=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100))
ggsave(p2, filename="014_Mp0_N_UE_atac.newUMAP.rm01091113.rm0825.rm24.chromVAR.JASPAR_chromVAR.averages.pheatmap.pdf", width=30, height=40, limitsize = FALSE)

p2 <-pheatmap::pheatmap(t_mat2, clustering_method="ward.D2", display_numbers = F, fontsize_row = 10, cluster_cols=T, cluster_rows=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100))
ggsave(p2, filename="014_Mp0_N_UE_atac.newUMAP.rm01091113.rm0825.rm24.chromVAR.JASPAR_chromVAR.averages.pheatmap.test.pdf", width=30, height=40, limitsize = FALSE)


hm.integrated.rm8.differential.motifs <- readRDS("new_signac_Hp0/hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.Rds")

hm.integrated.rm8.differential.motifs.id <- left_join(hm.integrated.rm8.differential.motifs, pfm.ID.name, by=c("gene" = "id"))
hm.integrated.rm8.differential.motifs.id %>% write_excel_csv("hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.gene.id.csv")

hm.integrated.rm8.differential.motifs.id <- read_csv("hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.gene.id.csv")
hm.integrated.rm8.differential.motifs.id %>% dplyr::select(cluster, name) %>% group_by(cluster) %>% mutate(id=row_number()) %>% pivot_wider(names_from = cluster, values_from = name) %>% as_tibble() %>% write_excel_csv("hm.integrated.rm8.differential.motifs.chromVARall.FindAllMarkers.gene.id.pivot_wider.csv")

motif_top10_per_cluster <- hm.integrated.rm8.differential.motifs.id %>% group_by(cluster) %>% top_n(wt=-p_val_adj, n=10) %>% ungroup() %>% select(c(`gene`,`name`)) %>% unique()
motif_top10_per_cluster # total 319 Motif

t_mat <- t(mat) %>% as.data.frame() %>% rownames_to_column(var="motif")
t_mat %>% head()

t_mat3 <- t_mat %>% filter(motif %in% pull(motif_top10_per_cluster, gene))

t_mat4 <- left_join(t_mat3, pfm.ID.name, by=c("motif" = "id")) %>% unite("id_name", c("motif","name")) %>% unique() %>% tibble()
t_mat4 %>% write_excel_csv("hm.integrated.rm8.differential.motifs.chromVAR.JASPAR_chromVAR.averages.pheatmap.csv")

motif_order <- c("MA0783.1_PKNOX2", "MA1572.1_TGIF2LY", "MA1571.1_TGIF2LX", "MA1638.1_HAND2", "MA0796.1_TGIF1", "MA1109.1_NEUROD1", "MA0669.1_NEUROG2", "MA0797.1_TGIF2", "MA1570.1_TFAP4(var.2)", "MA1618.1_Ptf1a", "MA1655.1_ZNF341", "MA1123.2_TWIST1", "MA1485.1_FERD3L", "MA1635.1_BHLHE22(var.2)", "MA0698.1_ZBTB18", "MA0500.2_MYOG", "MA0048.2_NHLH1", "MA0521.1_Tcf12", "MA0816.1_Ascl2", "MA1472.1_BHLHA15(var.2)", "MA0832.1_Tcf21", "MA0665.1_MSC", "MA1100.2_ASCL1", "MA0091.1_TAL1::TCF3", "MA0691.1_TFAP4", "MA0667.1_MYF6", "MA1619.1_Ptf1a(var.2)", "MA1641.1_MYF5", "MA1468.1_ATOH7", "MA1467.1_ATOH1(var.2)", "MA0909.2_HOXD13", "MA0687.1_SPIC", "MA0600.2_RFX2", "MA0510.2_RFX5", "MA1529.1_NHLH2", "MA1118.1_SIX1", "MA1119.1_SIX2", "MA1503.1_HOXB9", "MA0911.1_Hoxa11", "MA0631.1_Six3", "MA0485.2_HOXC9", "MA1473.1_CDX4", "MA1506.1_HOXD10", "MA0633.1_Twist2", "MA0465.2_CDX2", "MA0775.1_MEIS3", "MA0873.1_HOXD12", "MA0907.1_HOXC13", "MA0651.1_HOXC11", "MA0908.1_HOXD11", "MA0906.1_HOXC12", "MA0483.1_Gfi1b", "MA0878.2_CDX1", "MA0650.2_HOXA13", "MA0913.2_HOXD9", "MA1113.2_PBX2", "MA1640.1_MEIS2(var.2)", "MA0899.1_HOXA10", "MA1418.1_IRF3", "MA0629.1_Rhox11", "MA0905.1_HOXC10", "MA1639.1_MEIS1(var.2)", "MA1623.1_Stat2", "MA1419.1_IRF4", "MA0149.1_EWSR1-FLI1", "MA1652.1_ZKSCAN5", "MA0497.1_MEF2C", "MA0080.5_SPI1", "MA0081.2_SPIB", "MA0695.1_ZBTB7C", "MA0774.1_MEIS2", "MA0051.1_IRF2", "MA1155.1_ZSCAN4", "MA0508.3_PRDM1", "MA1524.1_MSGN1", "MA0528.2_ZNF263", "MA0050.2_IRF1", "MA1642.1_NEUROG2(var.2)", "MA1647.1_PRDM4", "MA0038.2_GFI1", "MA0056.2_MZF1", "MA0057.1_MZF1(var.2)", "MA0901.2_HOXB13", "MA1487.1_FOXE1", "MA1606.1_Foxf1", "MA0613.1_FOXG1", "MA0030.1_FOXF2", "MA0846.1_FOXC2", "MA0851.1_Foxj3", "MA0614.1_Foxj2", "MA0481.3_FOXP1", "MA0845.1_FOXB1", "MA1607.1_Foxl2", "MA0498.2_MEIS1", "MA0593.1_FOXP2", "MA0070.1_PBX1", "MA0847.2_FOXD2", "MA1103.2_FOXK2", "MA0852.2_FOXK1", "MA0032.2_FOXC1", "MA0052.4_MEF2A", "MA1630.1_Znf281", "MA0084.1_SRY", "MA0768.1_LEF1", "MA1497.1_HOXA6", "MA1421.1_TCF7L1", "MA1627.1_Wt1", "MA1522.1_MAZ", "MA0069.1_PAX6", "MA1116.1_RBPJ", "MA0781.1_PAX9", "MA0523.1_TCF7L2", "MA0769.2_TCF7", "MA0779.1_PAX1", "MA1505.1_HOXC8", "MA0634.1_ALX3", "MA1481.1_DRGX", "MA0674.1_NKX6-1", "MA1519.1_LHX5", "MA0912.2_HOXD3", "MA0675.1_NKX6-2", "MA0681.2_PHOX2B", "MA0880.1_Dlx3", "MA1500.1_HOXB6", "MA1501.1_HOXB7", "MA0709.1_Msx3", "MA0708.1_MSX2", "MA0782.2_PKNOX1", "MA0471.2_E2F6", "MA0758.1_E2F7", "MA0632.2_TCFL5", "MA0078.1_Sox17", "MA1530.1_NKX6-3", "MA0865.1_E2F8", "MA0151.1_Arid3a", "MA1653.1_ZNF148", "MA0442.2_SOX10", "MA0087.1_Sox5", "MA0515.1_Sox6", "MA1125.1_ZNF384", "MA0867.2_SOX4", "MA0514.1_Sox3", "MA1152.1_SOX15", "MA0517.1_STAT1::STAT2", "MA0868.2_SOX8", "MA0854.1_Alx1", "MA0725.1_VSX1", "MA0702.2_LMX1A", "MA0892.1_GSX1", "MA0703.2_LMX1B", "MA0661.1_MEOX1", "MA0707.1_MNX1", "MA1463.1_ARGFX", "MA0601.1_Arid3b", "MA0793.1_POU6F2", "MA0630.1_SHOX", "MA0135.1_Lhx3", "MA0075.3_PRRX2", "MA0704.1_Lhx4", "MA0662.1_MIXL1", "MA0723.1_VAX2", "MA0706.1_MEOX2", "MA0718.1_RAX", "MA0700.2_LHX2", "MA0031.1_FOXD1", "MA0480.1_Foxo1", "MA0850.1_FOXP3", "MA0848.1_FOXO4", "MA0033.2_FOXL1", "MA0042.2_FOXI1", "MA1601.1_ZNF75D", "MA0157.2_FOXO3", "MA0732.1_EGR3", "MA0849.1_FOXO6", "MA0495.3_MAFF", "MA1521.1_MAFA", "MA0659.2_MAFG", "MA0842.2_NRL", "MA0117.2_Mafb", "MA1520.1_MAF", "MA0461.2_Atoh1", "MA0018.4_CREB1", "MA0656.1_JDP2(var.2)", "MA0488.1_JUN", "MA1593.1_ZNF317", "MA1121.1_TEAD2", "MA0090.3_TEAD1", "MA0808.1_TEAD3", "MA0809.2_TEAD4", "MA0156.2_FEV", "MA0028.2_ELK1", "MA0759.1_ELK3", "MA0763.1_ETV3", "MA0474.2_ERG", "MA0475.2_FLI1", "MA0760.1_ERF", "MA1143.1_FOSL1::JUND(var.2)", "MA0098.3_ETS1", "MA0765.2_ETV5", "MA0645.1_ETV6", "MA0076.2_ELK4", "MA1483.1_ELF2", "MA0609.2_CREM", "MA0762.1_ETV2", "MA0641.1_ELF4", "MA0750.2_ZBTB7A", "MA1133.1_JUN::JUNB(var.2)", "MA0065.2_Pparg::Rxra", "MA0859.1_Rarg", "MA1574.1_THRB", "MA0855.1_RXRB", "MA0728.1_Nr2f6(var.2)", "MA0115.1_NR1H2::RXRA", "MA1537.1_NR2F1(var.2)", "MA0856.1_RXRG", "MA0017.2_NR2F1", "MA1550.1_PPARD", "MA0677.1_Nr2f6", "MA1148.1_PPARA::RXRA", "MA0857.1_Rarb", "MA0093.3_USF1", "MA0871.2_TFEC", "MA0512.2_Rxra", "MA0663.1_MLX", "MA0526.3_USF2", "MA1494.1_HNF4A(var.2)", "MA0664.1_MLXIPL", "MA0620.3_MITF", "MA0114.4_HNF4A", "MA0692.1_TFEB", "MA0504.1_NR2C2", "MA0484.2_HNF4G", "MA1112.2_NR4A1", "MA0831.2_TFE3", "MA1541.1_NR6A1", "MA1534.1_NR1I3", "MA1111.1_NR2F2", "MA0837.1_CEBPE", "MA0838.1_CEBPG", "MA0160.1_NR4A2", "MA0505.1_Nr5a2", "MA0829.2_SREBF1(var.2)", "MA0466.2_CEBPB", "MA0046.2_HNF1A", "MA0153.2_HNF1B", "MA0858.1_Rarb(var.2)", "MA0496.3_MAFK", "MA0591.1_Bach1::Mafk", "MA0501.1_MAF::NFE2", "MA0592.3_ESRRA", "MA0141.3_ESRRB", "MA0643.1_Esrrg", "MA0089.2_NFE2L1", "MA1633.1_BACH1", "MA0841.1_NFE2", "MA0491.2_JUND", "MA0655.1_JDP2", "MA1134.1_FOS::JUNB", "MA1144.1_FOSL2::JUND", "MA1130.1_FOSL2::JUN", "MA1137.1_FOSL1::JUNB", "MA1128.1_FOSL1::JUN", "MA0490.2_JUNB", "MA0478.1_FOSL2", "MA0476.1_FOS", "MA0835.2_BATF3", "MA1135.1_FOSB::JUNB", "MA1141.1_FOS::JUND", "MA0477.2_FOSL1", "MA1622.1_Smad2::Smad3", "MA1138.1_FOSL2::JUNB", "MA1634.1_BATF", "MA0462.2_BATF::JUN", "MA0099.3_FOS::JUN", "MA1101.2_BACH2", "MA0489.1_JUN(var.2)", "MA1132.1_JUN::JUNB", "MA1142.1_FOSL1::JUND", "MA1540.1_NR5A1", "MA0683.1_POU4F2", "MA1499.1_HOXB4", "MA0784.1_POU1F1", "MA0788.1_POU3F3", "MA0627.2_POU2F3", "MA0790.1_POU4F1", "MA0786.1_POU3F1", "MA0787.1_POU3F2", "MA1115.1_POU5F1", "MA0792.1_POU5F1B", "MA0785.1_POU2F1", "MA1471.1_BARX2", "MA0628.1_POU6F1", "MA0071.1_RORA", "MA0740.1_KLF14", "MA0685.1_SP4", "MA1517.1_KLF6", "MA1516.1_KLF3", "MA0079.4_SP1", "MA0145.3_TFCP2", "MA1105.2_GRHL2", "MA0789.1_POU3F4", "MA0647.1_GRHL1", "MA1131.1_FOSL2::JUN(var.2)", "MA0003.4_TFAP2A", "MA0745.2_SNAI2", "MA1620.1_Ptf1a(var.3)", "MA0499.2_MYOD1", "MA1621.1_Rbpjl", "MA1631.1_ASCL1(var.2)", "MA1558.1_SNAI1", "MA1559.1_SNAI3", "MA0830.2_TCF4", "MA0522.3_TCF3", "MA0744.2_SCRT2", "MA1648.1_TCF12(var.2)", "MA0820.1_FIGLA", "MA0103.3_ZEB1", "MA0743.2_SCRT1", "MA0742.1_Klf12", "MA0037.3_GATA3", "MA0766.2_GATA5", "MA1104.2_GATA6", "MA0036.3_GATA2", "MA0482.2_GATA4", "MA1528.1_NFIX(var.2)", "MA0119.1_NFIC::TLX1", "MA1643.1_NFIB", "MA1527.1_NFIC(var.2)", "MA0598.3_EHF", "MA0473.3_ELF1", "MA0640.2_ELF3", "MA0493.1_Klf1", "MA0136.2_ELF5", "MA0062.3_GABPA", "MA0761.2_ETV1", "MA1508.1_IKZF1", "MA0764.2_ETV4", "MA0525.2_TP63", "MA0861.1_TP73", "MA0106.3_TP53")


t_mat4 %>% head()
p2 <-pheatmap::pheatmap((t_mat4 %>% column_to_rownames("id_name"))[motif_order,], clustering_method="ward.D2", display_numbers = F, fontsize_row = 10, cluster_cols=F, cluster_rows=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100))
ggsave(p2, filename="014_hm.integrated.rm8.differential.motifs.chromVAR.JASPAR_chromVAR.averages.pheatmap.final.eps", width=30, height=40, limitsize = FALSE)


# Fig S2A.


# Human and Mouse N-UE correlation plot

geneTable <- read_csv("new_seurat_Mp0_Hp0_integrated/geneTrans.txt")
Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm <- readRDS("new_seurat_Mp0/015_Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.res1.0.RNA_norm.Rds")

mouseToHumanGeneName <- left_join(
    x=Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm@assays$RNA@counts@Dimnames[[1]] %>% as_tibble(),
    y=geneTable %>% dplyr::select(c(Human.Symbol, Mouse.Symbol)), by=c("value" = "Mouse.Symbol")
) %>% mutate(newnames = ifelse(is.na(Human.Symbol), value, Human.Symbol)) %>% pull(newnames)

mouseToHumanGeneName.integrated <- left_join(
    x=Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm@assays$integrated@data@Dimnames[[1]] %>% as_tibble(),
    y=geneTable %>% dplyr::select(c(Human.Symbol, Mouse.Symbol)), by=c("value" = "Mouse.Symbol")
) %>% mutate(newnames = ifelse(is.na(Human.Symbol), value, Human.Symbol)) %>% pull(newnames)


# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)
#_____________________________________

Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm <- RenameGenesSeurat(obj = Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, newnames = mouseToHumanGeneName)

Mp0_RNA.av.exp <- AverageExpression(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm)$RNA
Mp0_RNA.av.exp %>% rownames() %>% head()
Mp0_RNA.av.exp %>% colnames() %>% head()

Mp0_RNA.av.exp.variablefeature <- log1p(Mp0_RNA.av.exp[mouseToHumanGeneName.integrated,])
Mp0_RNA.av.exp.all <- log1p(Mp0_RNA.av.exp[mouseToHumanGeneName,])


# Hp0_10_SCT_sex.combined.res1.0.rm152931.rm33 <- readRDS("new_seurat_Hp0/Hp0_10_SCT_sex.combined.res1.0.rm152931.rm33.RNA_Norm.Rds")
Hp0_10_SCT_sex.combined.res1.0.rm152931.rm33 <- readRDS("newest_seurat_Hp0/105_integrated.RNA_norm.res1.6.rm1927.res1.0.23.clusters.RNA_norm.Rds")
Hp0_RNA.av.exp <- AverageExpression(Hp0_10_SCT_sex.combined.res1.0.rm152931.rm33)$RNA

Hp0_RNA.av.exp.variablefeature <- log1p(Hp0_RNA.av.exp[VariableFeatures(Hp0_10_SCT_sex.combined.res1.0.rm152931.rm33[['integrated']]),])
Hp0_RNA.av.exp.all <- log1p(Hp0_RNA.av.exp)


Mp0_Hp0_integrated_variable_feature_intersect <- intersect(rownames(Mp0_RNA.av.exp.variablefeature), rownames(Hp0_RNA.av.exp.variablefeature))
Mp0_Hp0_integrated_variable_feature_intersect %>% length()
Mp0_Hp0_RNA_intersect <- intersect(rownames(Mp0_RNA.av.exp), rownames(Hp0_RNA.av.exp))
Mp0_Hp0_RNA_intersect %>% length()

Hp0_RNA.av.exp.variablefeature.selected <- Hp0_RNA.av.exp.variablefeature[Mp0_Hp0_integrated_variable_feature_intersect,]
Mp0_RNA.av.exp.variablefeature.selected <- Mp0_RNA.av.exp.variablefeature[Mp0_Hp0_integrated_variable_feature_intersect,]


# Nephron_ureteric
df1 <- Mp0_RNA.av.exp.variablefeature.selected %>% dplyr::select(c('0','2','1','18','13','11','19','16','9','7','5','20','6','23','4','8','3','10','21','25','7','15','12','22','17','14','26','24'))

#  %>% dplyr::select(c('6','9','13','25','19','18','10','11','12','28','15','22','24','31','34'))
df1 %>% dim()
# df2 <- Hp0_RNA.av.exp.variablefeature.selected %>% dplyr::select(c('0','5','17','18','13','1','11','23','20','9','29','7','21','6','3','24','8','2','22'))


df2 <- Hp0_RNA.av.exp.variablefeature.selected %>% dplyr::select(c('2','3','11','17','10','9','8','0','13','12','18','14','21','6','7','15','5','19','20','4','1','16','22'))

#  %>% dplyr::select(c('4','2','14','22','16','17','0','21','8','12','7','23','9','24','5','20','25','6','1','26'))
df2 %>% dim()

#all
res <- data.frame()
for (i in 1:ncol(df1)) {
    for (j in 1:ncol(df2)) {
        res[i,j] <- cor(df1[,i], df2[,j], method = "spearman")
    }
}
rownames(res) <- colnames(df1)
colnames(res) <- colnames(df2)
ggsave(pheatmap(res, cluster_cols=T, cluster_rows=T),filename="correlation.heatmap.nephue.only.test.2.pdf")
# ggsave(pheatmap(res, cluster_cols=F, cluster_rows=F),filename="correlation.heatmap.nephue.test.no.cluster.pdf")

ggsave(pheatmap::pheatmap(res, display_numbers = F, cluster_cols=F, cluster_rows=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100)),filename="correlation.heatmap.nephue.only.test.no.cluster.2.pdf", width=7, height=7)



