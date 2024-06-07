# 014_developmental_trajectory_analysis.r

# Note: This script includes the podocyte trajectory analysis.
# The trajectory analysis of other nephrogenic lineages is done similarly,
# with only differences in subclustering of proximal tubule cells and distal nephron cells.

# Load necessary libraries
library(monocle)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(magrittr)

# Load mouse data
mouse_Pax8 <- readRDS("Mp0_N_UE_nephrogenic_Pax8_scRNA-seq.Rds")

# Load DEGs
mouse_Pax8.FindAllMarkers <- read_csv("Mp0_N_UE_nephrogenic_Pax8_scRNA-seq.DEG.FindAllMarkers.csv")

# Pivot and save DEGs
mouse_Pax8.FindAllMarkers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  as_tibble() %>%
  write_excel_csv("Mp0_N_UE_nephrogenic_Pax8_scRNA-seq.DEG.FindAllMarkers.pivot_wider.csv")

# Set working directory and load human data
setwd("newest_seurat_Hp0")
human_Pax8 <- integrated.RNA_norm.N.rm17.RNA_norm.N.rm02316.RNA_norm.rm61016.rm23.RNA_norm

setwd("final_trajectory_analysis")

# Function to find most correlated genes
GetMostCorrelated <- function(Srt, gene) {  
  comparable.genes <- t(Srt@scale.data[apply(Srt@scale.data, 1, sd) > 0, ])
  ans <- cor(Srt@scale.data[gene, ], comparable.genes, method = "pearson")
  names(ans) <- colnames(comparable.genes)
  sort(ans, decreasing = TRUE)
}

# Function to run monocle analysis with specified markers
monocle_running <- function(seurat_object.filtered, seurat_object.filtered.AllMarkers, trajectory_plot_filename, state_plot_filename, pseudotime_plot_filename) {
  ordering_genes.table <- seurat_object.filtered.AllMarkers %>%
    group_by(cluster) %>%
    arrange(p_val_adj) %>%
    top_frac(n = 0.1, wt = -p_val_adj)
  
  ordering_genes <- ordering_genes.table %>%
    pull(gene) %>%
    unique()
  
  print(length(ordering_genes))

  seurat_object.filtered.monocle2 <- as.CellDataSet(seurat_object.filtered)
  
  seurat_object.filtered.monocle2 %<>%
    monocle::setOrderingFilter(ordering_genes) %>%
    estimateSizeFactors() %>%
    estimateDispersions() %>%
    reduceDimension(method = 'DDRTree') %>%
    orderCells()

  p1 <- plot_cell_trajectory(seurat_object.filtered.monocle2, color_by = "seurat_clusters")
  ggsave(p1, filename = trajectory_plot_filename, width = 6, height = 4)
  ggsave(plot_cell_trajectory(seurat_object.filtered.monocle2, color_by = "State"), filename = state_plot_filename)
  
  p1 <- plot_pseudotime_heatmap(seurat_object.filtered.monocle2[ordering_genes, ], return_heatmap = TRUE, show_rownames = TRUE, num_clusters = 6)
  ggsave(p1, filename = pseudotime_plot_filename, height = round(length(ordering_genes) / 10), width = 6, limitsize = FALSE)

  return(seurat_object.filtered.monocle2)
}

# Function to run monocle analysis with variable features
monocle_running2 <- function(seurat_object.filtered, trajectory_plot_filename, state_plot_filename, pseudotime_plot_filename) {
  DefaultAssay(seurat_object.filtered) <- "RNA"
  seurat_object.filtered <- NormalizeData(seurat_object.filtered)
  ordering_genes <- VariableFeatures(seurat_object.filtered)
  
  seurat_object.filtered.monocle2 <- as.CellDataSet(seurat_object.filtered)
  
  seurat_object.filtered.monocle2 %<>%
    monocle::setOrderingFilter(ordering_genes) %>%
    estimateSizeFactors() %>%
    estimateDispersions() %>%
    reduceDimension(method = 'DDRTree') %>%
    orderCells()

  p1 <- plot_cell_trajectory(seurat_object.filtered.monocle2, color_by = "seurat_clusters")
  ggsave(p1, filename = trajectory_plot_filename, width = 6, height = 4)
  ggsave(plot_cell_trajectory(seurat_object.filtered.monocle2, color_by = "State"), filename = state_plot_filename)
  
  p1 <- plot_pseudotime_heatmap(seurat_object.filtered.monocle2[ordering_genes, ], return_heatmap = TRUE, show_rownames = TRUE, num_clusters = 6)
  ggsave(p1, filename = pseudotime_plot_filename, height = round(length(ordering_genes) / 10), width = 6, limitsize = FALSE)

  return(seurat_object.filtered.monocle2)
}

# Podocyte trajectory analysis (human/mouse)
# Mouse podocyte
mouse_Pax8.podocyte_PEC <- subset(mouse_Pax8, idents = c(1, 8, 9, 2, 10, 0, 4, 14))

# Preprocess and run PCA/UMAP
DefaultAssay(mouse_Pax8.podocyte_PEC) <- "integrated"
mouse_Pax8.podocyte_PEC %<>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 15) %>%
  FindClusters(resolution = 1.6) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

# Normalize data and run feature plots
DefaultAssay(mouse_Pax8.podocyte_PEC) <- "RNA"
mouse_Pax8.podocyte_PEC <- NormalizeData(mouse_Pax8.podocyte_PEC)

# Further subsetting and preprocessing
mouse_Pax8.podocyte_PEC.podocyte <- subset(mouse_Pax8.podocyte_PEC, idents = c(8, 15, 5, 14, 6, 10, 1, 3, 4, 0, 2, 17, 9, 16, 18))

DefaultAssay(mouse_Pax8.podocyte_PEC.podocyte) <- "integrated"
mouse_Pax8.podocyte_PEC.podocyte %<>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 15) %>%
  FindClusters(resolution = 1.6) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

# Run monocle analysis
mouse_Pax8.podocyte_PEC.podocyte.monocle2 <- monocle_running(
  mouse_Pax8.podocyte_PEC.podocyte,
  mouse_Pax8.podocyte_PEC.podocyte.AllMarkers,
  "mouse_Pax8.podocyte_PEC.podocyte.trajectory.pdf",
  "mouse_Pax8.podocyte_PEC.podocyte.state.pdf",
  "mouse_Pax8.podocyte_PEC.podocyte.pseudotime.pdf"
)

# Remove specific clusters and preprocess
mouse_Pax8.podocyte_PEC.podocyte.rm7111217 <- subset(mouse_Pax8.podocyte_PEC.podocyte, idents = c(7, 11, 12, 17), invert = TRUE)

DefaultAssay(mouse_Pax8.podocyte_PEC.podocyte.rm7111217) <- "integrated"
mouse_Pax8.podocyte_PEC.podocyte.rm7111217 %<>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 15) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

# Normalize data and find markers
DefaultAssay(mouse_Pax8.podocyte_PEC.podocyte.rm7111217) <- "RNA"
mouse_Pax8.podocyte_PEC.podocyte <- NormalizeData(mouse_Pax8.podocyte_PEC.podocyte.rm7111217)

mouse_Pax8.podocyte_PEC.podocyte.rm7111217.AllMarkers <- FindAllMarkers(mouse_Pax8.podocyte_PEC.podocyte.rm7111217, only.pos = TRUE)
mouse_Pax8.podocyte_PEC.podocyte.rm7111217.AllMarkers %>%
  write_excel_csv("mouse_Pax8.podocyte_PEC.podocyte.rm7111217.AllMarkers.csv")

mouse_Pax8.podocyte_PEC.podocyte.rm7111217.AllMarkers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  as_tibble() %>%
  write_excel_csv("mouse_Pax8.podocyte_PEC.podocyte.rm7111217.AllMarkers.pivot_wider.csv")

# Get most correlated genes
GetMostCorrelated(mouse_Pax8.podocyte_PEC.podocyte.rm7111217@assays$integrated, "Anxa1") %>%
  write.csv("mouse_Pax8.podocyte_PEC.podocyte.rm7111217.Anxa1.corr.csv")

# Differential gene test
m.diff_test_res <- differentialGeneTest(mouse_Pax8.podocyte_PEC.podocyte.monocle2, fullModelFormulaStr = "~sm.ns(Pseudotime)")
m.diff_test_res.out <- m.diff_test_res[, c("gene_short_name", "pval", "qval")]
m.diff_test_res.out %>% write_excel_csv("mouse_Pax8.podocyte_PEC.podocyte.monocle2.csv")

# Human podocyte
human_Pax8.podocyte_PEC <- subset(human_Pax8, idents = c(25, 6, 24, 11, 19, 5, 2, 26, 23, 0))

DefaultAssay(human_Pax8.podocyte_PEC) <- "integrated"
human_Pax8.podocyte_PEC %<>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 15) %>%
  FindClusters(resolution = 1.6) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

DefaultAssay(human_Pax8.podocyte_PEC) <- "RNA"
human_Pax8.podocyte_PEC <- NormalizeData(human_Pax8.podocyte_PEC)
human_Pax8.podocyte_PEC <- FindVariableFeatures(human_Pax8.podocyte_PEC)

# Find all markers and filter by clusters
human_Pax8.podocyte_PEC.AllMarkers <- FindAllMarkers(human_Pax8.podocyte_PEC, only.pos = TRUE)
human_Pax8.podocyte_PEC.AllMarkers %>%
  filter(cluster %in% c(21, 9)) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = -p_val_adj)

# Further subsetting and preprocessing
human_Pax8.podocyte_PEC.podocyte <- subset(human_Pax8.podocyte_PEC, idents = c(20, 13, 8, 14, 18, 10, 3, 15, 16, 12, 2, 0, 1, 5, 6, 17))

DefaultAssay(human_Pax8.podocyte_PEC.podocyte) <- "integrated"
human_Pax8.podocyte_PEC.podocyte %<>%
  ScaleData(verbose = TRUE) %>%
  RunPCA(npcs = 30, verbose = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 15) %>%
  FindClusters(resolution = 1.6) %>%
  RunUMAP(reduction = "pca", dims = 1:20)

# Normalize data and run feature plots
DefaultAssay(human_Pax8.podocyte_PEC.podocyte) <- "RNA"
human_Pax8.podocyte_PEC.podocyte <- NormalizeData(human_Pax8.podocyte_PEC.podocyte)
human_Pax8.podocyte_PEC.podocyte <- ScaleData(human_Pax8.podocyte_PEC.podocyte)
ggsave(FeaturePlot(human_Pax8.podocyte_PEC.podocyte, features = toupper(selected_genes), label = TRUE, repel = TRUE, order = TRUE, min.cutoff = "q10", max.cutoff = "q90"), filename = "004_human_Pax8.podocyte_PEC.podocyte.featureplots.pdf", width = 7*3, height = 7*4)

# Find all markers and save
human_Pax8.podocyte_PEC.podocyte.AllMarkers <- FindAllMarkers(human_Pax8.podocyte_PEC.podocyte, only.pos = TRUE)
human_Pax8.podocyte_PEC.podocyte.AllMarkers %>%
  write_excel_csv("human_Pax8.podocyte_PEC.podocyte.AllMarkers.csv")

# Run monocle analysis
human_Pax8.podocyte_PEC.podocyte.monocle2 <- monocle_running(
  human_Pax8.podocyte_PEC.podocyte,
  human_Pax8.podocyte_PEC.podocyte.AllMarkers,
  "human_Pax8.podocyte_PEC.podocyte.trajectory.pdf",
  "human_Pax8.podocyte_PEC.podocyte.state.pdf",
  "human_Pax8.podocyte_PEC.podocyte.pseudotime.pdf"
)

# Differential gene test
diff_test_res <- differentialGeneTest(human_Pax8.podocyte_PEC.podocyte.monocle2, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res.out <- diff_test_res[, c("gene_short_name", "pval", "qval")]
diff_test_res.out %>% write_excel_csv("human_Pax8.podocyte_PEC.podocyte.monocle2.csv")

# Get most correlated genes for human podocyte
OLFM3_corr <- GetMostCorrelated(human_Pax8.podocyte_PEC.podocyte@assays$RNA, "OLFM3")
PCDH9_corr <- GetMostCorrelated(human_Pax8.podocyte_PEC.podocyte@assays$RNA, "PCDH9")
CDH11_corr <- GetMostCorrelated(human_Pax8.podocyte_PEC.podocyte@assays$RNA, "CDH11")

# Identify more transient genes
more_transient <- union(
  union(
    OLFM3_corr %>% as.data.frame() %>% filter(`.` > 0.25) %>% rownames(),
    PCDH9_corr %>% as.data.frame() %>% filter(`.` > 0.25) %>% rownames()
  ),
  CDH11_corr %>% as.data.frame() %>% filter(`.` > 0.25) %>% rownames()
)
more_transient
