# Figure 1.

# Fig 1B.
UMAPPlot(object=Mp0_13i_3x3_SCT.remove.10.15.21, label=T, repel=T)

# Fig 1D.
DefaultAssay(Mp0_13i_3x3_SCT.remove.10.15.21)  <- "RNA"
Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm <- NormalizeData(Mp0_13i_3x3_SCT.remove.10.15.21, verbose = T)
Idents(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm) <- factor(Idents(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm), levels= c(
# NPC
6, 9, 25, 19,
# iNPC
13, 11, 15, 12, 18,
# Pod
10,
# PT
28,
# DT/LOH
22,
# UE
24, 31, 34,
# Vasc
0, 17, 26, 23,
# Immune
21, 27, 29,
# Interst
20, 8, 30, 2, 4, 3, 5, 16, 1, 14, 7, 32,
# Neuronal
33
))

p <- VlnPlot(object=Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, features=c(
#NPC
"Cited1", "Six2",
#iNPC
"Bmper", "Wnt4", "Lhx1",
#Podocyte
"Mafb", "Nphs2", #"Nphs1", #"Zbtb7c",
#PT
#"Slc7a13", #"Slc27a2", #"Slc22a6",
"Hnf4a", #"Spp2", #"Aldob", "Slc34a1", #"Lrp2",
#DT/loH
"Slc12a1", #"Tfap2b", #"Tmem213",
#PC
"Calb1", "Aqp2", #"Fxyd4",
#IC
"Foxi1", #"Atp6v1g3", #"Atp6v0d2",
#DMC
"Upk1b", #"Fxyd3", #"Sfn",
#Vascular
"Aplnr", "Plvap", "Aqp1", "Sox17",
#Imm
"C1qc", "Ccr2", "S100a9",
#Interstitial
"Eln", "Smoc2", "Clca3a1", "Fibin", #"Igfbp5",
"Ren1", "Gata3", "Myh11", "Actg2",
#Neural
"Sox10", "Foxd3"
), ncol = 1, pt.size=0)
ggsave(p, filename="Mp0_13i_3x3_SCT.remove.10.15.21.final_v9.split.markers.pdf", width=10, height=60, limitsize = FALSE)

# Fig 1C.
mouse_final_vlnplot_features <-c("Eya1", "Bmper", "Mafb", "Nphs2", "Hnf4a", "Cubn", "Slc34a1", "Pou3f3", "Slc12a1", "Umod", "Tmem52b", "Elf5", "Krt19", "Foxq1", "Egfl7", "Plvap", "Aplnr", "Ikzf1", "Cd36", "Pdgfrb", "Pdgfra", "Dcn", "Postn", "Ntn1", "Dkk2", "Tbx18")

ggsave(
    Seurat::VlnPlot(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, 
    features=mouse_final_vlnplot_features,
    pt.size=0, ncol=1),
    filename="023_Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm.mouse_final_vlnplot_features.vlnplot.pdf",
    width=12,
    height=3*length(mouse_final_vlnplot_features),
    limitsize = F
)

# Fig 1H.
atac.newUMAP.rm01091113.rm0825.rm24 <- readRDS("atac.newUMAP.rm01091113.rm0825.rm24.Rds")

Idents(atac.newUMAP.rm01091113.rm0825.rm24) <- factor(Idents(atac.newUMAP.rm01091113.rm0825.rm24),
 levels= c(0,4,2,6,10,9,24,5,13,3,1,19,18,8,15,12,7,17,22,16,14,11,21,20,23))

DefaultAssay(atac.newUMAP.rm01091113.rm0825.rm24) <- "RNA"
p <- DotPlot(atac.newUMAP.rm01091113.rm0825.rm24, features=c("Eya1", "Six2", "Spock2", "Meox2", "Zbtb7c", "Wt1", "Podxl", "R3hdml", "Nphs2", "Srgap1", "Magi2", "Col4a2", "Akap12", "Osr2", "Jag1", "Hnf4a", "Aldob", "Hnf4g", "Slc5a12", "Slc34a1", "Spp2", "Serpina6", "Aadat", "Aqp1", "Greb1", "Lef1", "Pax2", "Tmem52b", "Mecom", "Tfap2a", "Trpv5", "Umod", "Irx1", "Irx2", "Wnt7b", "Pou3f3", "Tfap2b", "Slc12a1", "Kcnj1", "Atp6v1b1", "Slc12a3", "Nos1", "Wnt11", "Gfra1", "Wnt9b", "Gata3", "Scnn1b", "Scnn1g", "Aqp2", "Fxyd4", "Aqp4", "Foxq1", "Foxi1", "Atp6v1g3", "Upk1b", "Upk3a"), col.min = 0) + coord_flip()

ggsave(p, filename="atac.newUMAP.rm01091113.rm0825.rm24.GA.dotplots.pdf", width=20, height=16, limitsize=F)

# Fig 1G.
hm.integrated.rm21.rm2526 <- readRDS("063_Hp0_BGI_NOV_filtered_hm_integrated.rm21.rm2526.novIntPeaks.res1.6.Rds")
ggsave(UMAPPlot(hm.integrated.rm21.rm2526, label=F, repel=F, raster=100000000)+NoLegend(), filename="060_Hp0_BGI_NOV_filtered_hm_integrated.rm21.rm2526.novIntPeaks.res1.6.UMAP.tiff")
ggsave(UMAPPlot(hm.integrated.rm21.rm2526, label=T, repel=T)+NoLegend(), filename="060_Hp0_BGI_NOV_filtered_hm_integrated.rm21.rm2526.novIntPeaks.res1.6.UMAP.eps")

# Fig 1I.
hm.integrated.rm21.rm2526 <- readRDS("063_Hp0_BGI_NOV_filtered_hm_integrated.rm21.rm2526.novIntPeaks.res1.6.Rds")
DefaultAssay(hm.integrated.rm21.rm2526) <- "RNA"
hm.integrated.rm21.rm2526 <- NormalizeData(hm.integrated.rm21.rm2526)

human_vlnplot_final_features <- c("EYA1", "BMPER", "MAFB", "NPHS2", "HNF4A", "CUBN", "SLC34A1", "POU3F3", "SLC12A1", "UMOD", "TMEM52B", "GATA3", "KRT19", "FOXQ1", "EGFL7", "PLVAP", "APLNR", "IKZF1", "CD36", "PDGFRB", "PDGFRA", "DCN", "POSTN", "NTN1", "DKK2", "TBX18")

Idents(hm.integrated.rm21.rm2526) <- factor(Idents(hm.integrated.rm21.rm2526), levels= c(
    0, 2, 16, 10, 15, 6, 20, 9, 19, 7, 17, 11, 4, 1, 3, 5, 21, 25, 14, 23, 24, 26, 22, 8, 13, 12, 
    18 
))

DefaultAssay(hm.integrated.rm21.rm2526) <- "RNA"
p <- DotPlot(hm.integrated.rm21.rm2526, features=rev(human_vlnplot_final_features),
  cols = c("lightgrey", "red"), dot.min=0.25, col.min=0.5) + coord_flip()
ggsave(p, filename="073_2_hm.integrated.rm21.rm2526.final.dotplot.markers.pdf", width=12, height=6)

# Fig 1J.
library(TFBSTools)
library(readr)
library(pheatmap)
library(RColorBrewer)

atac.newUMAP.rm01091113.rm0825.rm24 <- readRDS("atac.newUMAP.rm01091113.rm0825.rm24.Rds")

atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_allJASPAR_chromVAR.averages <- AverageExpression(atac.newUMAP.rm01091113.rm0825.rm24, assays =  "common_ATAC_allJASPAR_chromVAR", return.seurat = TRUE)

JASPAR2020.10090 <- readRDS("JASPAR2020.10090.pfm.Rds")
JASPAR2020.10090.pfm.ID.name <- tibble(id=ID(JASPAR2020.10090), name=name(JASPAR2020.10090))
JASPAR2020.9606 <- readRDS("JASPAR2020.9606.pfm.Rds")
JASPAR2020.9606.pfm.ID.name <- tibble(id=ID(JASPAR2020.9606), name=name(JASPAR2020.9606))
JASPAR2020.10090.9606.pfm.ID.name <- dplyr::bind_rows(JASPAR2020.10090.pfm.ID.name, JASPAR2020.9606.pfm.ID.name)
motifs_features <- atac.newUMAP.rm01091113.rm0825.rm24.differential.motifs %>% pull(gene) %>% unique()
motifs_features <- intersect(motifs_features, JASPAR2020.10090.9606.pfm.ID.name %>% pull(id))

p <- DoHeatmap(atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_allJASPAR_chromVAR.averages, features = motifs_features, size = 3, draw.lines = F)

ggsave(p, filename="atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_9606_10090_JASPAR_chromVAR.averages.heatmap.pdf", width=10, height=40, limitsize = FALSE)
mat <- p$data %>% select(Feature, Cell, Expression) %>% pivot_wider(names_from=Feature, values_from=Expression) %>% column_to_rownames("Cell")

pfm <- readRDS("JASPAR2020.all.pfm.Rds")
pfm.ID.name <- tibble(id=ID(pfm), name=name(pfm))

t_mat <- left_join(t(mat) %>% as.data.frame()%>% rownames_to_column(var="motif"), pfm.ID.name, by=c("motif" = "id"))
t_mat <- t_mat %>% unite("id_name", c("motif","name"))
rownames(t_mat) <- t_mat$id_name
t_mat <- t_mat %>% select(!c(id_name))

t_mat %>% rownames_to_column("Motif") %>% write_excel_csv("atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_9606_10090_JASPAR_chromVAR.averages.pheatmap.csv")

t_mat <- read_csv("atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_9606_10090_JASPAR_chromVAR.averages.pheatmap.csv")
t_mat <- t_mat %>% as.data.frame() %>% tibble::column_to_rownames("Motif")

t_mat_Hoxs <- t_mat %>% filter(Motif %in% Hoxs) %>% as.data.frame() %>% tibble::column_to_rownames("Motif")
t_mat_Foxs <- t_mat %>% filter(Motif %in% Foxs) %>% as.data.frame() %>% tibble::column_to_rownames("Motif")

p1 <-pheatmap(t_mat, clustering_method="ward.D2", fontsize_row = 8, cluster_cols=T, color=colorRamp(c("yellow", "green", "blue"))(100)) 
p1 <-pheatmap(t_mat, clustering_method="ward.D2", fontsize_row = 8, cluster_cols=T, color=colorRampPalette(c("yellow", "green", "blue"))(100)) 

p1 <-pheatmap(t_mat, clustering_method="ward.D2", fontsize_row = 1, cluster_cols=T, color=colorRampPalette(brewer.pal(9,"YlGnBu"))(100)) 
ggsave(p1, filename="atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_9606_10090_JASPAR_chromVAR.averages.pheatmap.color.pdf", width=16, height=24, limitsize = FALSE)

p2 <-pheatmap(t_mat, clustering_method="ward.D2", fontsize_row = 8, cluster_cols=T, display_numbers=T, number_format="%.2f")
ggsave(p2, filename="atac.newUMAP.rm01091113.rm0825.rm24.common_ATAC_9606_10090_JASPAR_chromVAR.averages.pheatmap.pdf", width=16, height=80, limitsize = FALSE)


# Fig 1K.
DefaultAssay(hm.integrated.rm21.rm2526) <- "novIntPeaks_JASPAR_chromVAR"
hm.integrated.rm21.rm2526.JASPAR.averages <- AverageExpression(hm.integrated.rm21.rm2526, assays =  "novIntPeaks_JASPAR_chromVAR", return.seurat = TRUE)

pfm.JASPAR2020 <- readRDS("JASPAR2020.9606.10090.pfm.Rds")

p <- DoHeatmap(hm.integrated.rm21.rm2526.JASPAR.averages, features = unique(ID(pfm.JASPAR2020)), size = 3, draw.lines = F)

mat <- p$data %>% dplyr::select(Feature, Cell, Expression) %>% pivot_wider(names_from=Feature, values_from=Expression) %>% column_to_rownames("Cell")

pfm.ID.name <- tibble(id=ID(pfm.JASPAR2020), name=name(pfm.JASPAR2020))

t_mat <- left_join(t(mat) %>% as.data.frame()%>% rownames_to_column(var="motif"), pfm.ID.name, by=c("motif" = "id"))
t_mat <- t_mat %>% unite("id_name", c("motif","name"))
t_mat <- t_mat %>% unique()
rownames(t_mat) <- t_mat$id_name
t_mat <- t_mat %>% dplyr::select(!c(id_name))

p2 <-pheatmap(t_mat, clustering_method="ward.D2", display_numbers = F, fontsize_row = 10, cluster_cols=F, color=colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"))(100))

ggsave(p2, filename="069_hm.integrated.rm21.rm2526.chromVAR.JASPAR_chromVAR.averages.pheatmap.pdf", width=30, height=40, limitsize = FALSE)

# Figure_S1.R

# Fig S1A.
ggsave(
    VlnPlot(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, c(
        "nCount_RNA", "nFeature_RNA", "percent.mt"),
    pt.size=0, ncol=1),
    filename="004_Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm.QC.vlnplot.pdf",
    width=12,
    height=3*3,
    limitsize = F
)

#plot as percent of cluster/proportion
ggsave(
    ggplot(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm@meta.data, aes(x=seurat_clusters, fill=replicate)) + geom_bar(position = "fill") + ggpubr::theme_pubr(),
    filename="006_Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm.geom_bar.replicate.pdf",
    width=12,
    height=6,
    limitsize=F
)

t <- Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm@meta.data %>% dplyr::select(c(seurat_clusters, replicate))
t$seurat_clusters <- as.character(t$seurat_clusters)
t$seurat_clusters <- factor(t$seurat_clusters, levels=c(6, 9, 25, 19, 13, 11, 18, 10, 12, 28, 15, 22, 31, 24, 34, 0, 17, 23, 26, 27, 29, 21, 20, 8, 30, 2, 4, 3, 5, 16, 1, 14, 7, 32, 33))

ggsave(
    ggplot(t, aes(x=seurat_clusters, fill=replicate)) + geom_bar(position = "fill") + ggpubr::theme_pubr(),
    filename="006_Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm.geom_bar.replicate.pdf",
    width=12,
    height=6,
    limitsize=F
)

# Fig S1B.


# Fig S1D.
ggsave(
    VlnPlot(hm.integrated.rm21.rm2526, c(
        "nCount_ATAC", "nFeature_ATAC", "nCount_common_ATAC", "nFeature_common_ATAC", "duplicate", "mitochondrial", "TSS_fragments", "peak_region_fragments", "on_target_fragments", "enhancer_region_fragments", "blacklist_region_fragments", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_percentile", "nCount_RNA", "nFeature_RNA"),
    pt.size=0),
    filename="062_Hp0_BGI_NOV_filtered_hm_integrated.rm21.rm2526.novIntPeaks.QC.pdf",
    width=20,
    height=16,
    limitsize = F
)

t <- hm.integrated.rm21.rm2526@meta.data %>% dplyr::select(c(seurat_clusters, dataset))
t$seurat_clusters <- as.character(t$seurat_clusters)
t$seurat_clusters <- factor(t$seurat_clusters, levels=c(0, 2, 16, 10, 15, 6, 20, 9, 19, 7, 17, 11, 1, 4, 3, 5, 21, 25, 14, 23, 24, 26, 18, 12, 13, 8, 22))

ggsave(
    ggplot(t, aes(x=seurat_clusters, fill=dataset)) + geom_bar(position = "fill") + ggpubr::theme_pubr(),
    filename="072_hm.integrated.rm21.rm2526.geom_bar.replicate.pdf",
    width=12,
    height=6,
    limitsize=F
)

ggsave(
    VlnPlot(hm.integrated.rm21.rm2526, c("nCount_ATAC", 
              "nFeature_ATAC",
              "peak_region_fragments",
              "blacklist_ratio"),
    pt.size=0, ncol=1),
    filename="073_hm.integrated.rm21.rm2526.QC.vlnplot.pdf",
    width=12,
    height=3*4,
    limitsize = F
)
