# 003_Human_fetal_RNA.r

# Define paths for BGI datasets
Hp0_BGI_M1_path <- "M_10_6_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_BGI_M2_path <- "M_11_1_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_BGI_M3_path <- "M_11_6_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_BGI_F1_path <- "F_11_3_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_BGI_F2_path <- "F_15_1_fresh_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_BGI_F3_path <- "F_15_2_frozen_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"

# Load and create Seurat objects for each BGI dataset
Hp0_BGI_M1 <- Read10X(data.dir = Hp0_BGI_M1_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")
Hp0_BGI_M2 <- Read10X(data.dir = Hp0_BGI_M2_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")
Hp0_BGI_M3 <- Read10X(data.dir = Hp0_BGI_M3_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")
Hp0_BGI_F1 <- Read10X(data.dir = Hp0_BGI_F1_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")
Hp0_BGI_F2 <- Read10X(data.dir = Hp0_BGI_F2_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")
Hp0_BGI_F3 <- Read10X(data.dir = Hp0_BGI_F3_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_BGI")

# Merge BGI datasets into a single Seurat object
Hp0_BGI_six <- merge(x = Hp0_BGI_M1, y = c(Hp0_BGI_M2, Hp0_BGI_M3, Hp0_BGI_F1, Hp0_BGI_F2, Hp0_BGI_F3),
                     add.cell.ids = c("BGI_M1", "BGI_M2", "BGI_M3", "BGI_F1", "BGI_F2", "BGI_F3"))

# Add metadata: replicate, sex, and mitochondrial percentage
replicate <- substring(rownames(Hp0_BGI_six@meta.data), 1, 6)
names(replicate) <- rownames(Hp0_BGI_six@meta.data)
Hp0_BGI_six <- AddMetaData(Hp0_BGI_six, metadata = replicate, col.name = "replicate")

sex <- substring(rownames(Hp0_BGI_six@meta.data), 5, 5)
names(sex) <- rownames(Hp0_BGI_six@meta.data)
Hp0_BGI_six <- AddMetaData(Hp0_BGI_six, metadata = sex, col.name = "sex")

Hp0_BGI_six[["percent.mt"]] <- PercentageFeatureSet(Hp0_BGI_six, pattern = "^MT-")

# Subset based on RNA features and mitochondrial percentage
Hp0_BGI_six <- subset(Hp0_BGI_six, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 35)
Hp0_BGI_six <- subset(Hp0_BGI_six, subset = nCount_RNA > 500 & nCount_RNA < 30000)

# Perform SCTransform normalization
Hp0_BGI_six_SCT <- SCTransform(Hp0_BGI_six, vars.to.regress = "nCount_RNA", verbose = TRUE)

# Run PCA and UMAP, and find neighbors and clusters
Hp0_BGI_six_SCT <- RunPCA(Hp0_BGI_six_SCT, npcs = 30, verbose = TRUE)
Hp0_BGI_six_SCT <- RunUMAP(Hp0_BGI_six_SCT, dims = 1:30, verbose = TRUE)
Hp0_BGI_six_SCT <- FindNeighbors(Hp0_BGI_six_SCT, reduction = "pca", dims = 1:20, k.param = 15, verbose = TRUE)
Hp0_BGI_six_SCT <- FindClusters(Hp0_BGI_six_SCT, resolution = 0.8, verbose = TRUE)

# Save UMAP plots
ggsave(UMAPPlot(Hp0_BGI_six_SCT, label = TRUE) + NoLegend(), filename = "009_Hp0_BGI_six_SCT.UMAP.pdf")
ggsave(UMAPPlot(Hp0_BGI_six_SCT, group.by = 'sex', label = TRUE) + NoLegend(), filename = "009_Hp0_BGI_six_SCT_sex.UMAP.sex.pdf")

# Define paths for NOV datasets
Hp0_NOV_F4_path <- "F_11_2_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_NOV_M4_path <- "M_17_1_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_NOV_M5_path <- "M_17_5_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_NOV_F5_path <- "F_17_1_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"
Hp0_NOV_F6_path <- "F_17_6_scRNA_5_0_1_include_introns/outs/filtered_feature_bc_matrix/"

# Load and create Seurat objects for each NOV dataset
Hp0_NOV_F4 <- Read10X(data.dir = Hp0_NOV_F4_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_NOV")
Hp0_NOV_M4 <- Read10X(data.dir = Hp0_NOV_M4_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_NOV")
Hp0_NOV_M5 <- Read10X(data.dir = Hp0_NOV_M5_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_NOV")
Hp0_NOV_F5 <- Read10X(data.dir = Hp0_NOV_F5_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_NOV")
Hp0_NOV_F6 <- Read10X(data.dir = Hp0_NOV_F6_path) %>% CreateSeuratObject(min.cells = 3, min.features = 300, project = "Hp0_NOV")

# Merge NOV datasets into a single Seurat object
Hp0_NOV_five <- merge(x = Hp0_NOV_F4, y = c(Hp0_NOV_M4, Hp0_NOV_M5, Hp0_NOV_F5, Hp0_NOV_F6),
                      add.cell.ids = c("NOV_F4", "NOV_M4", "NOV_M5", "NOV_F5", "NOV_F6"))

# Add metadata: replicate, sex, and mitochondrial percentage
replicate <- substring(rownames(Hp0_NOV_five@meta.data), 1, 6)
names(replicate) <- rownames(Hp0_NOV_five@meta.data)
Hp0_NOV_five <- AddMetaData(Hp0_NOV_five, metadata = replicate, col.name = "replicate")

sex <- substring(rownames(Hp0_NOV_five@meta.data), 5, 5)
names(sex) <- rownames(Hp0_NOV_five@meta.data)
Hp0_NOV_five <- AddMetaData(Hp0_NOV_five, metadata = sex, col.name = "sex")

Hp0_NOV_five[["percent.mt"]] <- PercentageFeatureSet(Hp0_NOV_five, pattern = "^MT-")

# Subset based on RNA features and mitochondrial percentage
Hp0_NOV_five <- subset(Hp0_NOV_five, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 35)
Hp0_NOV_five <- subset(Hp0_NOV_five, subset = nCount_RNA > 500 & nCount_RNA < 30000)

# Further subset excluding NOV_F4
Hp0_NOV_five_subset_F4 <- subset(Hp0_NOV_five, subset = replicate == "NOV_F4", invert = TRUE)

# Perform SCTransform normalization
Hp0_NOV_five_subset_F4_SCT <- SCTransform(Hp0_NOV_five_subset_F4, vars.to.regress = "nCount_RNA", verbose = TRUE)
Hp0_NOV_five_subset_F4_SCT <- RunPCA(Hp0_NOV_five_subset_F4_SCT, npcs = 30, verbose = TRUE)
Hp0_NOV_five_subset_F4_SCT <- RunUMAP(Hp0_NOV_five_subset_F4_SCT, dims = 1:30, verbose = TRUE)
Hp0_NOV_five_subset_F4_SCT <- FindNeighbors(Hp0_NOV_five_subset_F4_SCT, reduction = "pca", dims = 1:20, k.param = 15, verbose = TRUE)
Hp0_NOV_five_subset_F4_SCT <- FindClusters(Hp0_NOV_five_subset_F4_SCT, resolution = 0.8, verbose = TRUE)

# Save UMAP plots
ggsave(UMAPPlot(Hp0_NOV_five_subset_F4_SCT, label = TRUE) + NoLegend(), filename = "009_Hp0_NOV_five_subset_F4_SCT.UMAP.pdf")
ggsave(UMAPPlot(Hp0_NOV_five_subset_F4_SCT, group.by = 'sex', label = TRUE) + NoLegend(), filename = "009_Hp0_NOV_five_subset_F4_SCT_sex.UMAP.sex.pdf")

# Integration of BGI and NOV datasets
# Prepare for integration
Hp0_10_SCT.features <- SelectIntegrationFeatures(object.list = list(Hp0_BGI_six_SCT, Hp0_NOV_five_subset_F4_SCT), nfeatures = 3000)

# Prepare objects for SCT integration
Hp0_10_SCT.list <- PrepSCTIntegration(object.list = list(Hp0_BGI_six_SCT, Hp0_NOV_five_subset_F4_SCT), anchor.features = Hp0_10_SCT.features, verbose = TRUE)

# Find integration anchors (This step may take a long time)
Hp0_10_SCT.anchors <- FindIntegrationAnchors(object.list = Hp0_10_SCT.list, anchor.features = Hp0_10_SCT.features, normalization.method = "SCT", dims = 1:20)

# Integrate data
Hp0_10_SCT.combined <- IntegrateData(anchorset = Hp0_10_SCT.anchors, normalization.method = "SCT", dims = 1:20)

# Run PCA and UMAP, and find neighbors and clusters for the integrated data
DefaultAssay(Hp0_10_SCT.combined) <- 'integrated'
Hp0_10_SCT.combined <- RunPCA(Hp0_10_SCT.combined, npcs = 30, verbose = FALSE)
Hp0_10_SCT.combined <- RunUMAP(Hp0_10_SCT.combined, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined <- FindNeighbors(Hp0_10_SCT.combined, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined <- FindClusters(Hp0_10_SCT.combined, resolution = 1.0)

# Save UMAP plots for the integrated data
ggsave(UMAPPlot(Hp0_10_SCT.combined, label = TRUE) + NoLegend(), filename = "075_Hp0_10_SCT.combined.UMAP.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined, group.by = "orig.ident", label = TRUE) + NoLegend(), filename = "076_Hp0_10_SCT.combined.UMAP.orig.ident.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined, group.by = "replicate", label = TRUE) + NoLegend(), filename = "077_Hp0_10_SCT.combined.UMAP.replicate.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined, group.by = "sex", label = TRUE) + NoLegend(), filename = "078_Hp0_10_SCT.combined.UMAP.sex.pdf")

# Save UMAP plots split by clusters
ggsave(UMAPPlot(Hp0_10_SCT.combined, split.by = "seurat_clusters", label = TRUE) + NoLegend(), filename = "145_Hp0_10_SCT.combined.UMAP.split.clusters.pdf", limitsize = FALSE, width = 7 * 35, height = 7)

# Save violin plots for QC metrics
p <- VlnPlot(object = Hp0_10_SCT.combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
ggsave(p, filename = "144_Hp0_10_SCT.combined.UMAP.Vlnplots.QC.pdf", width = 7 * 3, height = 7)

# Set default assay to RNA and normalize data
DefaultAssay(Hp0_10_SCT.combined) <- "RNA"
Hp0_10_SCT.combined <- NormalizeData(Hp0_10_SCT.combined)

# Remove clusters 0 and 22 and re-process the data
Hp0_10_SCT.combined.rm022 <- subset(Hp0_10_SCT.combined, idents = c(0, 22), invert = TRUE)
DefaultAssay(Hp0_10_SCT.combined.rm022) <- 'integrated'
Hp0_10_SCT.combined.rm022 <- RunPCA(Hp0_10_SCT.combined.rm022, npcs = 30, verbose = FALSE)
Hp0_10_SCT.combined.rm022 <- RunUMAP(Hp0_10_SCT.combined.rm022, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022 <- FindNeighbors(Hp0_10_SCT.combined.rm022, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022 <- FindClusters(Hp0_10_SCT.combined.rm022, resolution = 1.0)

# Save UMAP plots for the subset without clusters 0 and 22
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022, label = TRUE) + NoLegend(), filename = "075_Hp0_10_SCT.combined.rm022.UMAP.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022, group.by = "orig.ident", label = TRUE) + NoLegend(), filename = "140_Hp0_10_SCT.combined.rm022.UMAP.orig.ident.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022, group.by = "replicate", label = TRUE) + NoLegend(), filename = "141_Hp0_10_SCT.combined.rm022.UMAP.replicate.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022, group.by = "sex", label = TRUE) + NoLegend(), filename = "142_Hp0_10_SCT.combined.rm022.UMAP.sex.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022, split.by = "seurat_clusters", label = TRUE) + NoLegend(), filename = "145_Hp0_10_SCT.combined.rm022.UMAP.split.clusters.pdf", limitsize = FALSE, width = 7 * 29, height = 7)

# Save violin plots for QC metrics
p <- VlnPlot(object = Hp0_10_SCT.combined.rm022, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
ggsave(p, filename = "144_Hp0_10_SCT.combined.rm022.Vlnplots.QC.pdf", width = 7 * 3, height = 7)

# Remove clusters 31 and 34 and re-process the data
Hp0_10_SCT.combined.rm022.rm3134 <- subset(Hp0_10_SCT.combined.rm022, idents = c(31, 34), invert = TRUE)
DefaultAssay(Hp0_10_SCT.combined.rm022.rm3134) <- 'integrated'
Hp0_10_SCT.combined.rm022.rm3134 <- RunPCA(Hp0_10_SCT.combined.rm022.rm3134, npcs = 30, verbose = FALSE)
Hp0_10_SCT.combined.rm022.rm3134 <- RunUMAP(Hp0_10_SCT.combined.rm022.rm3134, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022.rm3134 <- FindNeighbors(Hp0_10_SCT.combined.rm022.rm3134, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022.rm3134 <- FindClusters(Hp0_10_SCT.combined.rm022.rm3134, resolution = 1.0)

# Save UMAP plots for the subset without clusters 31 and 34
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134, label = TRUE) + NoLegend(), filename = "075_Hp0_10_SCT.combined.rm022.rm3134.UMAP.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134, group.by = "orig.ident", label = TRUE) + NoLegend(), filename = "140_Hp0_10_SCT.combined.rm022.rm3134.UMAP.orig.ident.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134, group.by = "replicate", label = TRUE) + NoLegend(), filename = "141_Hp0_10_SCT.combined.rm022.rm3134.UMAP.replicate.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134, group.by = "sex", label = TRUE) + NoLegend(), filename = "142_Hp0_10_SCT.combined.rm022.rm3134.UMAP.sex.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134, split.by = "seurat_clusters", label = TRUE) + NoLegend(), filename = "145_Hp0_10_SCT.combined.rm022.rm3134.UMAP.split.clusters.pdf", limitsize = FALSE, width = 7 * 29, height = 7)

# Save violin plots for QC metrics
p <- VlnPlot(object = Hp0_10_SCT.combined.rm022.rm3134, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
ggsave(p, filename = "144_Hp0_10_SCT.combined.rm022.rm3134.Vlnplots.QC.pdf", width = 7 * 3, height = 7)

# Remove clusters 30, 31, and 32 and re-process the data
Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- subset(Hp0_10_SCT.combined.rm022.rm3134, idents = c(30, 31, 32), invert = TRUE)
DefaultAssay(Hp0_10_SCT.combined.rm022.rm3134.rm303132) <- 'integrated'
Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- RunPCA(Hp0_10_SCT.combined.rm022.rm3134.rm303132, npcs = 30, verbose = FALSE)
Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- RunUMAP(Hp0_10_SCT.combined.rm022.rm3134.rm303132, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- FindNeighbors(Hp0_10_SCT.combined.rm022.rm3134.rm303132, reduction = "pca", dims = 1:30)
Hp0_10_SCT.combined.rm022.rm3134.rm303132 <- FindClusters(Hp0_10_SCT.combined.rm022.rm3134.rm303132, resolution = 1.0)

# Save UMAP plots for the subset without clusters 30, 31, and 32
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, label = TRUE) + NoLegend(), filename = "075_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.pdf")
Idents(Hp0_10_SCT.combined.rm022.rm3134.rm303132) <- factor(Idents(Hp0_10_SCT.combined.rm022.rm3134.rm303132), levels = c(
  4, 2, 14, 22, 16, 21, 17, 0, 8, 12, 7, 20, 23, 9, 24, 5, 25, 6, 1, 26, 30, 13, 31, 29, 18, 10, 3, 28, 19, 15, 11, 27
))

ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, label = TRUE, raster = 10000000) + NoLegend(), filename = "075_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.eps")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, group.by = "orig.ident", label = TRUE) + NoLegend(), filename = "140_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.orig.ident.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, group.by = "replicate", label = TRUE) + NoLegend(), filename = "141_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.replicate.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, group.by = "sex", label = TRUE) + NoLegend(), filename = "142_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.sex.pdf")
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, split.by = "seurat_clusters", label = TRUE) + NoLegend(), filename = "145_Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.split.clusters.pdf", limitsize = FALSE, width = 7 * 29, height = 7)

# Save violin plots for QC metrics
p <- VlnPlot(object = Hp0_10_SCT.combined.rm022.rm3134.rm303132, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
ggsave(p, filename = "144_Hp0_10_SCT.combined.rm022.rm3134.rm303132.Vlnplots.QC.pdf", width = 7 * 3, height = 7)

# N-UE subclustering

# Define clusters for N-UE subclustering
N_UE_clusters <- c(4, 2, 14, 22, 16, 21, 17, 0, 8, 12, 7, 20, 23, 9, 24, 5, 25, 6, 1, 26)

# Identify cell barcodes for BGI and NOV datasets in N-UE clusters
N_UE.BGI.cell.barcodes <- WhichCells(subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, idents = N_UE_clusters, subset = orig.ident == "Hp0_BGI"))
N_UE.NOV.cell.barcodes <- WhichCells(subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, idents = N_UE_clusters, subset = orig.ident == "Hp0_NOV"))

# Subset BGI and NOV datasets based on N-UE cell barcodes
Hp0_BGI_six_SCT.N_UE <- subset(Hp0_BGI_six_SCT, cells = N_UE.BGI.cell.barcodes)
Hp0_NOV_five_subset_F4.N_UE <- subset(Hp0_NOV_five_subset_F4, cells = N_UE.NOV.cell.barcodes)

# Set default assay to RNA and normalize data
DefaultAssay(Hp0_BGI_six_SCT.N_UE) <- "RNA"
DefaultAssay(Hp0_NOV_five_subset_F4.N_UE) <- "RNA"
Hp0_BGI_six_SCT.N_UE.RNA_norm <- NormalizeData(Hp0_BGI_six_SCT.N_UE, verbose = TRUE)
Hp0_NOV_five_subset_F4.N_UE.RNA_norm <- NormalizeData(Hp0_NOV_five_subset_F4.N_UE, verbose = TRUE)

# Find variable features
Hp0_BGI_six_SCT.N_UE.RNA_norm <- FindVariableFeatures(Hp0_BGI_six_SCT.N_UE.RNA_norm, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
Hp0_NOV_five_subset_F4.N_UE.RNA_norm <- FindVariableFeatures(Hp0_NOV_five_subset_F4.N_UE.RNA_norm, selection.method = "vst", nfeatures = 3000, verbose = TRUE)

# Find integration anchors and integrate data
integrated.anchors <- FindIntegrationAnchors(object.list = list(Hp0_BGI_six_SCT.N_UE.RNA_norm, Hp0_NOV_five_subset_F4.N_UE.RNA_norm), dims = 1:20)
integrated <- IntegrateData(anchorset = integrated.anchors, dims = 1:20)

# Scale data, run PCA, find neighbors, find clusters, and run UMAP
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = TRUE)
integrated <- RunPCA(integrated, npcs = 30, verbose = TRUE)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20, k.param = 15)
integrated <- FindClusters(integrated, resolution = 1.6)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)

# Normalize RNA data in the integrated object
DefaultAssay(integrated) <- "RNA"
integrated.RNA_norm <- NormalizeData(integrated, verbose = TRUE)

# Save UMAP plots
ggsave(UMAPPlot(integrated.RNA_norm, label = TRUE) + NoLegend(), filename = "106_Hp0_10.combined.N_UE.RNA_norm.UMAP.pdf")
ggsave(UMAPPlot(integrated.RNA_norm, label = TRUE) + NoLegend(), filename = "075_integrated.RNA_norm.UMAP.pdf")
ggsave(UMAPPlot(integrated.RNA_norm, group.by = "orig.ident", label = TRUE) + NoLegend(), filename = "076_integrated.RNA_norm.UMAP.orig.ident.pdf")
ggsave(UMAPPlot(integrated.RNA_norm, group.by = "replicate", label = TRUE) + NoLegend(), filename = "077_integrated.RNA_norm.UMAP.replicate.pdf")
ggsave(UMAPPlot(integrated.RNA_norm, group.by = "sex", label = TRUE) + NoLegend(), filename = "078_integrated.RNA_norm.UMAP.sex.pdf")
ggsave(UMAPPlot(integrated.RNA_norm, split.by = "seurat_clusters", label = TRUE) + NoLegend(), filename = "079_integrated.RNA_norm.UMAP.split.clusters.pdf", limitsize = FALSE, width = 7 * 29, height = 7)

# Save violin plots for QC metrics
p <- VlnPlot(object = integrated.RNA_norm, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
ggsave(p, filename = "080_integrated.RNA_norm.Vlnplots.QC.pdf", width = 7 * 3, height = 7)

# Remove clusters 19 and 27 and re-process the data
integrated.RNA_norm.res1.6.rm1927 <- subset(integrated.RNA_norm, idents = c(19, 27), invert = TRUE)
DefaultAssay(integrated.RNA_norm.res1.6.rm1927) <- "integrated"
integrated.RNA_norm.res1.6.rm1927 <- ScaleData(integrated.RNA_norm.res1.6.rm1927, verbose = TRUE)
integrated.RNA_norm.res1.6.rm1927 <- RunPCA(integrated.RNA_norm.res1.6.rm1927, npcs = 30, verbose = TRUE)
integrated.RNA_norm.res1.6.rm1927 <- FindClusters(integrated.RNA_norm.res1.6.rm1927, resolution = 1.0)
integrated.RNA_norm.res1.6.rm1927 <- RunUMAP(integrated.RNA_norm.res1.6.rm1927, reduction = "pca", dims = 1:20)

# Normalize RNA data in the new subset
DefaultAssay(integrated.RNA_norm.res1.6.rm1927) <- "RNA"
integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm <- NormalizeData(integrated.RNA_norm.res1.6.rm1927, verbose = TRUE)
