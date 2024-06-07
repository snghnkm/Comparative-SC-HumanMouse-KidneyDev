# Mp0_RNA_13i_3x3_SCT
# Load necessary libraries for the analysis

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(uwot)

# Load 2018 Mp0 Samples using the 10X Genomics v2 chemistry, Reference Genome v3, Cell Ranger v3
# Read 10X data for each sample and create Seurat objects

# Load old Mp0 3x3 data
Mp0_F1 <- Read10X(data.dir ="old_Mp0_3x3/F1_3x3/filtered_feature_bc_matrix")
Mp0_F2 <- Read10X(data.dir ="old_Mp0_3x3/F2_3x3/filtered_feature_bc_matrix")
Mp0_F3 <- Read10X(data.dir ="old_Mp0_3x3/F3_3x3/filtered_feature_bc_matrix")

Mp0_M5 <- Read10X(data.dir ="old_Mp0_3x3/M5_3x3/filtered_feature_bc_matrix")
Mp0_M6 <- Read10X(data.dir ="old_Mp0_3x3/M6_3x3/filtered_feature_bc_matrix")
Mp0_M7 <- Read10X(data.dir ="old_Mp0_3x3/M7_3x3/filtered_feature_bc_matrix")
Mp0_M4 <- Read10X(data.dir ="old_Mp0_3x3/M4_3x3/filtered_feature_bc_matrix")

# Create Seurat objects for each sample with minimum cell and feature thresholds
F1_obj <- CreateSeuratObject(Mp0_F1, min.cells = 3, min.features = 300, project = "p0_18")
F2_obj <- CreateSeuratObject(Mp0_F2, min.cells = 3, min.features = 300, project = "p0_18")
F3_obj <- CreateSeuratObject(Mp0_F3, min.cells = 3, min.features = 300, project = "p0_18")

M4_obj <- CreateSeuratObject(Mp0_M4, min.cells = 3, min.features = 300, project = "p0_18")
M5_obj <- CreateSeuratObject(Mp0_M5, min.cells = 3, min.features = 300, project = "p0_18")
M6_obj <- CreateSeuratObject(Mp0_M6, min.cells = 3, min.features = 300, project = "p0_18")
M7_obj <- CreateSeuratObject(Mp0_M7, min.cells = 3, min.features = 300, project = "p0_18")

# Merge individual Seurat objects into one combined object
Mp0_sev <- merge(x = F1_obj, y = c(F2_obj,F3_obj,M4_obj,M5_obj,M6_obj,M7_obj), add.cell.ids = c("F1","F2","F3","M4","M5","M6","M7"))

# Add metadata for replicate, sex, and mitochondrial gene percentage to the combined Seurat object
replicate <- substring(rownames(Mp0_sev@meta.data), 1,3)
names(replicate) <- rownames(Mp0_sev@meta.data)
Mp0_sev <- AddMetaData(Mp0_sev, metadata = replicate, col.name = "replicate") 

sex <- substring(rownames(Mp0_sev@meta.data), 1,1)
names(sex) <- rownames(Mp0_sev@meta.data)
Mp0_sev <- AddMetaData(Mp0_sev, metadata = sex, col.name = "sex")

Mp0_sev [["percent.mt"]] <- PercentageFeatureSet(Mp0_sev , pattern = "^mt-")

Mp0_sev  <-subset(Mp0_sev , subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 35)
Mp0_sev  <-subset(Mp0_sev , subset = nCount_RNA > 500 & nCount_RNA < 30000)

# Load and process new 2019 Mp0 samples using 10X Genomics v3 chemistry, Reference Genome v3, and Cell Ranger v3 (sequencing runs 1 and 2)
# Load new Mp0 r1+r2 data
Mp0_F4 <- Read10X(data.dir ="new Mp0 r1+r2/F4/filtered_feature_bc_matrix")
Mp0_F5 <- Read10X(data.dir ="new Mp0 r1+r2/F5/filtered_feature_bc_matrix")
Mp0_F6 <- Read10X(data.dir ="new Mp0 r1+r2/F6/filtered_feature_bc_matrix")
Mp0_M8 <- Read10X(data.dir ="new Mp0 r1+r2/M8/filtered_feature_bc_matrix")
Mp0_M9 <- Read10X(data.dir ="new Mp0 r1+r2/M9/filtered_feature_bc_matrix")
Mp0_M10 <- Read10X(data.dir ="new Mp0 r1+r2/M10/filtered_feature_bc_matrix")

F4_obj <- CreateSeuratObject(Mp0_F4, min.cells = 3, min.features = 300, project = "p0_19")
F5_obj <- CreateSeuratObject(Mp0_F5, min.cells = 3, min.features = 300, project = "p0_19")
F6_obj <- CreateSeuratObject(Mp0_F6, min.cells = 3, min.features = 300, project = "p0_19")
M8_obj <- CreateSeuratObject(Mp0_M8, min.cells = 3, min.features = 300, project = "p0_19")
M9_obj <- CreateSeuratObject(Mp0_M9, min.cells = 3, min.features = 300, project = "p0_19")
M10_obj <- CreateSeuratObject(Mp0_M10, min.cells = 3, min.features = 300, project = "p0_19")

Mp0_six <- merge(x = F4_obj, y = c(F5_obj,F6_obj,M8_obj,M9_obj,M10_obj), add.cell.ids = c("F4","F5","F6","M8","M9","M10"))

# Calculate the percentage of mitochondrial genes and add it to the metadata
replicate <- substring(rownames(Mp0_six@meta.data), 1,3)
names(replicate) <- rownames(Mp0_six@meta.data)
Mp0_six <- AddMetaData(Mp0_six, metadata = replicate, col.name = "replicate") 

sex <- substring(rownames(Mp0_six@meta.data), 1,1)
names(sex) <- rownames(Mp0_six@meta.data)
Mp0_six <- AddMetaData(Mp0_six, metadata = sex, col.name = "sex")

Mp0_six [["percent.mt"]] <- PercentageFeatureSet(Mp0_six , pattern = "^mt-")

Mp0_six  <-subset(Mp0_six , subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 35)
Mp0_six  <-subset(Mp0_six , subset = nCount_RNA > 500 & nCount_RNA < 30000)

# SCTransform normalization and data integration
# Run SCTransform to normalize the Mp0_six dataset

Mp0_sev_SCT  <- SCTransform(Mp0_sev, vars.to.regress = "nCount_RNA", verbose = FALSE)
Mp0_six_SCT  <- SCTransform(Mp0_six, vars.to.regress = "nCount_RNA", verbose = FALSE)


# Prepare for integration of two datasets (Mp0_sev_SCT and Mp0_six_SCT) using Seurat's integration functions
# Select integration features from both datasets, specifying the number of features to use
Mp0_13_SCT.features <- SelectIntegrationFeatures(object.list = list(Mp0_sev_SCT, Mp0_six_SCT), nfeatures = 3000)

# Prepare the datasets for integration by scaling and normalizing the data
# Find integration anchors between the two datasets
# Run the data integration using the identified anchors
Mp0_13_SCT.list <- PrepSCTIntegration(object.list = list(Mp0_sev_SCT, Mp0_six_SCT), anchor.features = Mp0_13_SCT.features, verbose = FALSE)
Mp0_13_SCT.anchors <- FindIntegrationAnchors(object.list = Mp0_13_SCT.list, anchor.features = Mp0_13_SCT.features, normalization.method = "SCT", dims = 1:20)
Mp0_13_SCT.combined <- IntegrateData(anchorset = Mp0_13_SCT.anchors, normalization.method = "SCT", dims = 1:20)

# Rename the integrated Seurat object
Mp0_13i_3x3_SCT <- Mp0_13_SCT.combined

# Set the default assay to "integrated" for downstream analysis
DefaultAssay(Mp0_13i_3x3_SCT) <- "integrated"

# Switch to the RNA assay for visualization purposes, as it contains the original uncorrected values
DefaultAssay(Mp0_13i_3x3_SCT) <- "RNA"
Mp0_13i_3x3_SCT_RNA_norm <- NormalizeData(Mp0_13i_3x3_SCT, verbose = FALSE)

Mp0_13i_3x3_SCT <- Mp0_13_SCT.combined
DefaultAssay(Mp0_13i_3x3_SCT) <- "integrated"

Mp0_13i_3x3_SCT <- RunPCA(Mp0_13i_3x3_SCT, npcs = 30, verbose = FALSE)
Mp0_13i_3x3_SCT <- RunTSNE(Mp0_13i_3x3_SCT, reduction = "pca", dims = 1:30)
Mp0_13i_3x3_SCT <- FindNeighbors(Mp0_13i_3x3_SCT, reduction = "pca", dims = 1:30)
Mp0_13i_3x3_SCT <- FindClusters(Mp0_13i_3x3_SCT, resolution = 1.0)

# Remove weak clusters
# Filter cells based on mitochondrial content and feature counts
Mp0_subset <- subset(Mp0_13i_3x3_SCT, idents = c(10, 15, 21), invert = TRUE)
Mp0_subset <- subset(Mp0_subset, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 30)

# Create a subset for Neph_UE (Nephron Upper Epithelial) cells
Mp0_subset_copy <- Mp0_subset
Mp0_NUE <- subset(Mp0_subset_copy, idents = c(5, 8, 22, 29, 13, 11, 18, 32, 12, 19, 9, 27, 31, 24, 35, 39))


# New Integrated Robj Created________________
#name shorten
Mp0_13i_3x3_SCT<- Mp0_13_SCT.combined

#After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note:  not that useful for gene expression plots, as only 2000 features
DefaultAssay(Mp0_13i_3x3_SCT)  <- "integrated"

#_Correct_mislabeled sexes
corrected_replicate <- Mp0_13i_3x3_SCT@meta.data %>% pull(replicate) %>%
    str_replace("F1_", "F1") %>%
    str_replace("F2_", "F2") %>%    
    str_replace("F3_", "F3") %>%    
    str_replace("M4_", "M4") %>%    
    str_replace("M5_", "M5") %>%    
    str_replace("M6_", "M6") %>%    
    str_replace("M7_", "M7") %>% 
    # Changed ------------------
    str_replace("M8_", "F4") %>%    
    str_replace("M9_", "F5") %>%    
    str_replace("M10", "F6") %>% 
    str_replace("F4_", "M8") %>%    
    str_replace("F5_", "M9") %>%    
    str_replace("F6_", "M10")    
    
corrected_sex <- corrected_replicate %>% str_sub(end=1)

# Xist___________________________
VlnPlot(object= Mp0_13i_3x3_SCT, features = c("Xist"), split.by = "sex", pt.size=0, split.plot = TRUE)
VlnPlot(object= Mp0_13i_3x3_SCT, features = c("Xist"), split.by = "corrected_sex", pt.size=0, split.plot = TRUE)
DoHeatmap(object= Mp0_13i_3x3_SCT, assay='SCT', features=c("Xist"), group.by="sex", slot = "data")
DoHeatmap(object= Mp0_13i_3x3_SCT, assay='SCT', features=c("Xist"), group.by="corrected_sex", slot = "data")

# Ddx3y_________________________
VlnPlot(object= Mp0_13i_3x3_SCT, features = c("Ddx3y"), split.by = "sex", pt.size=0, split.plot = TRUE)
VlnPlot(object= Mp0_13i_3x3_SCT, features = c("Ddx3y"), split.by = "corrected_sex", pt.size=0, split.plot = TRUE)
DoHeatmap(object= Mp0_13i_3x3_SCT, assay='SCT', features=c("Ddx3y"), group.by="sex", slot = "data")
DoHeatmap(object= Mp0_13i_3x3_SCT, assay='SCT', features=c("Ddx3y"), group.by="corrected_sex", slot = "data")

# Remove three weak cluster and reclustering with the same 2000 features

Mp0_13i_3x3_SCT.remove.10.15.21 <- subset(Mp0_13i_3x3_SCT, idents = c("10", "15", "21"), invert = TRUE)

DefaultAssay(Mp0_13i_3x3_SCT.remove.10.15.21) <- 'integrated'
Mp0_13i_3x3_SCT.remove.10.15.21 <- RunPCA(Mp0_13i_3x3_SCT.remove.10.15.21, npcs = 30, verbose = FALSE)
Mp0_13i_3x3_SCT.remove.10.15.21 <- RunUMAP(Mp0_13i_3x3_SCT.remove.10.15.21, reduction = "pca", dims = 1:30)
Mp0_13i_3x3_SCT.remove.10.15.21 <- FindNeighbors(Mp0_13i_3x3_SCT.remove.10.15.21, reduction = "pca", dims = 1:30)
Mp0_13i_3x3_SCT.remove.10.15.21 <- FindClusters(Mp0_13i_3x3_SCT.remove.10.15.21, resolution = 1.0)

# Remove three weak clusters (10, 15, 21) and re-cluster with the same 2000 features
Mp0_13i_3x3_SCT.remove.10.15.21 <- subset(Mp0_13i_3x3_SCT, idents = c("10", "15", "21"), invert = TRUE)

# Set the default assay to 'integrated'
DefaultAssay(Mp0_13i_3x3_SCT.remove.10.15.21) <- 'integrated'

# Perform PCA on the data
Mp0_13i_3x3_SCT.remove.10.15.21 <- RunPCA(Mp0_13i_3x3_SCT.remove.10.15.21, npcs = 30, verbose = FALSE)

# Run UMAP for dimensionality reduction
Mp0_13i_3x3_SCT.remove.10.15.21 <- RunUMAP(Mp0_13i_3x3_SCT.remove.10.15.21, reduction = "pca", dims = 1:30)

# Find neighbors in the PCA space
Mp0_13i_3x3_SCT.remove.10.15.21 <- FindNeighbors(Mp0_13i_3x3_SCT.remove.10.15.21, reduction = "pca", dims = 1:30)

# Find clusters in the data with resolution 1.0
Mp0_13i_3x3_SCT.remove.10.15.21 <- FindClusters(Mp0_13i_3x3_SCT.remove.10.15.21, resolution = 1.0)

# Subclustering analysis
# Load required data objects
load("./Mp0_sev_SCT.Robj")
load("./Mp0_six_SCT.Robj")

# Set the default assay to 'RNA' for both datasets
DefaultAssay(Mp0_sev_SCT) <- "RNA"
DefaultAssay(Mp0_six_SCT) <- "RNA"

# Normalize the RNA data for both datasets
Mp0_sev_SCT.RNA_norm <- NormalizeData(Mp0_sev_SCT, verbose = TRUE)
Mp0_six_SCT.RNA_norm <- NormalizeData(Mp0_six_SCT, verbose = TRUE)

# Load pre-processed RNA normalized data
Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm <- readRDS("Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm.Rds")
Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm

# Define cell barcodes for different clusters and conditions
# N-UE clusters
N_UE.p0_18.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(6, 9, 13, 25, 19, 11, 18, 10, 12, 15, 28, 22, 24, 34, 31), subset = orig.ident == "p0_18"))
N_UE.p0_19.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(6, 9, 13, 25, 19, 11, 18, 10, 12, 15, 28, 22, 24, 34, 31), subset = orig.ident == "p0_19"))

Mp0_sev_SCT.RNA_norm.N_UE <- subset(Mp0_sev_SCT.RNA_norm, cells = N_UE.p0_18.cell.barcodes)
Mp0_six_SCT.RNA_norm.N_UE <- subset(Mp0_six_SCT.RNA_norm, cells = N_UE.p0_19.cell.barcodes)

# Vascular clusters
Vascular.p0_18.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(0, 17, 26, 23), subset = orig.ident == "p0_18"))
Vascular.p0_19.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(0, 17, 26, 23), subset = orig.ident == "p0_19"))

Mp0_sev_SCT.RNA_norm.Vascular <- subset(Mp0_sev_SCT.RNA_norm, cells = Vascular.p0_18.cell.barcodes)
Mp0_six_SCT.RNA_norm.Vascular <- subset(Mp0_six_SCT.RNA_norm, cells = Vascular.p0_19.cell.barcodes)

# Immune clusters
Immune.p0_18.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(21, 27, 29), subset = orig.ident == "p0_18"))
Immune.p0_19.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(21, 27, 29), subset = orig.ident == "p0_19"))

Mp0_sev_SCT.RNA_norm.Immune <- subset(Mp0_sev_SCT.RNA_norm, cells = Immune.p0_18.cell.barcodes)
Mp0_six_SCT.RNA_norm.Immune <- subset(Mp0_six_SCT.RNA_norm, cells = Immune.p0_19.cell.barcodes)

# Interstitial clusters
Interstitial.p0_18.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(1, 2, 3, 4, 8, 5, 16, 14, 30, 20, 7, 32), subset = orig.ident == "p0_18"))
Interstitial.p0_19.cell.barcodes <- WhichCells(subset(Mp0_13i_3x3_SCT.remove.10.15.21.RNA_norm, idents = c(1, 2, 3, 4, 8, 5, 16, 14, 30, 20, 7, 32), subset = orig.ident == "p0_19"))

Mp0_sev_SCT.RNA_norm.Interstitial <- subset(Mp0_sev_SCT.RNA_norm, cells = Interstitial.p0_18.cell.barcodes)
Mp0_six_SCT.RNA_norm.Interstitial <- subset(Mp0_six_SCT.RNA_norm, cells = Interstitial.p0_19.cell.barcodes)

# Integration helper function
integration <- function (p0_18, p0_19) {
    DefaultAssay(p0_18) <- "RNA"
    DefaultAssay(p0_19) <- "RNA"
    
    p0_18 <- FindVariableFeatures(p0_18, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
    p0_19 <- FindVariableFeatures(p0_19, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
    
    integrated.anchors <- FindIntegrationAnchors(object.list = c(p0_18, p0_19), dims = 1:20)
    integrated <- IntegrateData(anchorset = integrated.anchors, dims = 1:20)
    
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated, verbose = TRUE)
    integrated <- RunPCA(integrated, npcs = 30, verbose = TRUE)
    integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20, k.param = 15)
    integrated <- FindClusters(integrated, resolution = 0.8)
    integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
    
    DefaultAssay(integrated) <- "RNA"
    integrated.RNA_norm <- NormalizeData(integrated, verbose = TRUE)
    return (integrated.RNA_norm)
}

# Integrate and normalize data for each cluster type
Mp0_13i_3x3_Norm.N_UE.RNA_norm         <- integration(Mp0_sev_SCT.RNA_norm.N_UE, Mp0_six_SCT.RNA_norm.N_UE)
Mp0_13i_3x3_Norm.Vascular.RNA_norm     <- integration(Mp0_sev_SCT.RNA_norm.Vascular, Mp0_six_SCT.RNA_norm.Vascular)
Mp0_13i_3x3_Norm.Interstitial.RNA_norm <- integration(Mp0_sev_SCT.RNA_norm.Interstitial, Mp0_six_SCT.RNA_norm.Interstitial)
Mp0_13i_3x3_Norm.Immune.RNA_norm       <- integration(Mp0_sev_SCT.RNA_norm.Immune, Mp0_six_SCT.RNA_norm.Immune)