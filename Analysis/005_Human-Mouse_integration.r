# 005_Human-Mouse_integration.r

### Human-mouse transcriptomics integration: Cross-species sc/snRNA-seq 

# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)

# Step 1: Load datasets and prepare for integration

# Load human and mouse datasets
Hp0_10_SCT_sex.combined.N_UE.RNA_norm <- readRDS("newest_seurat_Hp0/105_integrated.RNA_norm.res1.6.rm1927.res0.8.RNA_norm.Rds")
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223 <- readRDS("newest_seurat_Mp0_Hp0_integrated/Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.Rds")

# Annotate species in metadata
Hp0_10_SCT_sex.combined.N_UE.RNA_norm$species <- "human"
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223$species <- "mouse"

# Load gene translation table
geneTable <- read_csv("new_seurat_Mp0_Hp0_integrated/geneTrans.txt")

# Create a mapping of mouse to human gene names
mouseToHumanGeneName <- left_join(
  x = Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223@assays$RNA@counts@Dimnames[[1]] %>% as_tibble(),
  y = geneTable %>% dplyr::select(c(Human.Symbol, Mouse.Symbol)),
  by = c("value" = "Mouse.Symbol")
) %>% mutate(newnames = ifelse(is.na(Human.Symbol), value, Human.Symbol)) %>% pull(newnames)

# Function to rename genes in a Seurat object
RenameGenesSeurat <- function(obj, newnames) {
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data, and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
  } else {
    print("Unequal gene sets: nrow(RNA) != nrow(newnames)")
  }
  
  obj@assays$RNA <- RNA
  return(obj)
}

# Rename genes in the mouse dataset
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223 <- RenameGenesSeurat(obj = Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, newnames = mouseToHumanGeneName)

# Subset datasets by origin
Hp0_10_SCT_sex.combined.N_UE.RNA_norm.BGI <- subset(Hp0_10_SCT_sex.combined.N_UE.RNA_norm, subset = orig.ident == "Hp0_BGI")
Hp0_10_SCT_sex.combined.N_UE.RNA_norm.NOV <- subset(Hp0_10_SCT_sex.combined.N_UE.RNA_norm, subset = orig.ident == "Hp0_NOV")
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_18 <- subset(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, subset = orig.ident == "p0_18")
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_19 <- subset(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223, subset = orig.ident == "p0_19")

# Normalize data using SCTransform
Hp0_10_SCT_sex.combined.N_UE.RNA_norm.BGI <- SCTransform(Hp0_10_SCT_sex.combined.N_UE.RNA_norm.BGI, vars.to.regress = "nCount_RNA", verbose = TRUE)
Hp0_10_SCT_sex.combined.N_UE.RNA_norm.NOV <- SCTransform(Hp0_10_SCT_sex.combined.N_UE.RNA_norm.NOV, vars.to.regress = "nCount_RNA", verbose = TRUE)
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_18 <- SCTransform(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_18, vars.to.regress = "nCount_RNA", verbose = TRUE)
Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_19 <- SCTransform(Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_19, vars.to.regress = "nCount_RNA", verbose = TRUE)

# Prepare list for integration
integration_list <- list(Hp0_10_SCT_sex.combined.N_UE.RNA_norm.BGI, Hp0_10_SCT_sex.combined.N_UE.RNA_norm.NOV, Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_18, Mp0_13i_3x3_Norm.N_UE.RNA_norm.rm2223.p0_19)

# Select integration features and prepare for integration
transcriptome.features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 3000)
transcriptome.list <- PrepSCTIntegration(object.list = integration_list, anchor.features = transcriptome.features)

# Find integration anchors and integrate data
transcriptome.anchors <- FindIntegrationAnchors(object.list = transcriptome.list, normalization.method = "SCT", anchor.features = transcriptome.features, verbose = TRUE, dims = 1:20)
transcriptome.integrated <- IntegrateData(anchorset = transcriptome.anchors, normalization.method = "SCT", dims = 1:20)

# Run PCA, UMAP, and clustering
DefaultAssay(transcriptome.integrated) <- 'integrated'
transcriptome.integrated <- RunPCA(transcriptome.integrated, verbose = FALSE)
transcriptome.integrated <- RunUMAP(transcriptome.integrated, dims = 1:30)
transcriptome.integrated <- FindNeighbors(transcriptome.integrated, reduction = "pca", dims = 1:20, k.param = 15, verbose = TRUE)
transcriptome.integrated <- FindClusters(transcriptome.integrated, resolution = 0.8, verbose = TRUE)

# Save plots
ggsave(DimPlot(transcriptome.integrated, group.by = "species"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.pdf")
ggsave(DimPlot(transcriptome.integrated, group.by = "replicate"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.by.replicate.pdf")
ggsave(DimPlot(transcriptome.integrated, group.by = "sex"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.by.sex.pdf")
ggsave(UMAPPlot(transcriptome.integrated, label = TRUE, repel = TRUE), filename = "010_Hp0_Mp0_integrated_trial5.UMAP.pdf")
ggsave(UMAPPlot(transcriptome.integrated, label = TRUE, repel = TRUE) + NoLegend(), filename = "010_Hp0_Mp0_integrated_trial5.UMAP.NoLegend.pdf")

# Re-cluster at higher resolution
transcriptome.integrated <- FindClusters(transcriptome.integrated, resolution = 1.0, verbose = TRUE)

# Save additional plots
ggsave(DimPlot(transcriptome.integrated, group.by = "species"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.by.species.1.0.pdf")
ggsave(DimPlot(transcriptome.integrated, group.by = "replicate"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.by.replicate.1.0.pdf")
ggsave(DimPlot(transcriptome.integrated, group.by = "sex"), filename = "009_Mp0_Hp0_RNAseq_integrated.2k.pca30.by.sex.1.0.pdf")
ggsave(UMAPPlot(transcriptome.integrated, label = TRUE, repel = TRUE), filename = "010_Hp0_Mp0_integrated_trial5.UMAP.1.0.pdf")
ggsave(UMAPPlot(transcriptome.integrated, label = TRUE, repel = TRUE) + NoLegend(), filename = "010_Hp0_Mp0_integrated_trial5.UMAP.NoLegend.1.0.pdf")

# Subset integrated object by species
transcriptome.integrated.mouse <- subset(transcriptome.integrated, subset = species == "mouse")
DefaultAssay(transcriptome.integrated.mouse) <- "RNA"
transcriptome.integrated.mouse <- NormalizeData(transcriptome.integrated.mouse)

transcriptome.integrated.human <- subset(transcriptome.integrated, subset = species == "human")
DefaultAssay(transcriptome.integrated.human) <- "RNA"
transcriptome.integrated.human <- NormalizeData(transcriptome.integrated.human)

# The integrated object is now prepared and can be further analyzed or saved as needed.

### All-inclusive analysis: snATAC-seq intergration on top of Cross-species sc/snRNA-seq

## Human snRNA-seq and snATAC-seq integration

# Load integrated human RNA data
transcriptome.integrated.human <- readRDS("newest_seurat_Mp0_Hp0_integrated/010_Hp0_Mp0_integrated_trial5.human.RNA_norm.Rds")
DefaultAssay(transcriptome.integrated.human) <- "integrated"
table(transcriptome.integrated.human@active.ident)

# Load human ATAC data
hp0_N_UE.ATAC <- readRDS("new_signac_Hp0/092_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.Rds")
DefaultAssay(hp0_N_UE.ATAC) <- "RNA"

# Find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = transcriptome.integrated.human,
  query = hp0_N_UE.ATAC,
  reference.assay = "integrated",
  query.assay = "RNA",
  reduction = 'cca',
  dims = 1:30
)

# Transfer data to get cell type predictions
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = transcriptome.integrated.human@active.ident,
  weight.reduction = hp0_N_UE.ATAC[["harmony"]],
  dims = 2:30
)

# Add cell type predictions to ATAC metadata
hp0_N_UE.ATAC <- AddMetaData(hp0_N_UE.ATAC, metadata = celltype.predictions)

# Update cluster identities and generate UMAP plots
Idents(object = hp0_N_UE.ATAC) <- 'novIntPeaks_snn_res.1.6'
p1 <- UMAPPlot(hp0_N_UE.ATAC, label = TRUE, repel = TRUE) + NoLegend()
Idents(object = hp0_N_UE.ATAC) <- 'predicted.id'
Idents(object = hp0_N_UE.ATAC) <- factor(x = Idents(hp0_N_UE.ATAC), levels = stringr::str_sort(levels(hp0_N_UE.ATAC), numeric = TRUE))
p2 <- UMAPPlot(hp0_N_UE.ATAC, label = TRUE, repel = TRUE) + NoLegend()
p3 <- UMAPPlot(transcriptome.integrated.human, label = TRUE, repel = TRUE) + NoLegend()

# Save combined UMAP plot
ggsave(p1 + p2 + p3, filename = "018_hp0_N_UE.predicted.id.pdf", width = 18, height = 6)

# Identify variable features
genes.use <- VariableFeatures(transcriptome.integrated.human)

# Transfer imputation data
refdata.integrated <- GetAssayData(transcriptome.integrated.human, assay = "integrated", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata.integrated, weight.reduction = hp0_N_UE.ATAC[["harmony"]])
hp0_N_UE.ATAC[["integrated"]] <- imputation

# Update metadata for RNA and ATAC data
hp0_N_UE.ATAC@meta.data$seurat_clusters <- NULL
transcriptome.integrated.human@meta.data$seurat_clusters <- transcriptome.integrated.human@active.ident
transcriptome.integrated.human$tech <- "RNA"
hp0_N_UE.ATAC$tech <- "ATAC"

# Merge RNA and ATAC data
coembed.subset <- merge(x = transcriptome.integrated.human, y = hp0_N_UE.ATAC)

# Re-cluster and run PCA and UMAP
DefaultAssay(coembed.subset) <- "integrated"
coembed.subset <- ScaleData(coembed.subset, features = genes.use, do.scale = TRUE)
coembed.subset <- RunPCA(coembed.subset, features = genes.use, verbose = TRUE)
coembed.subset <- RunUMAP(coembed.subset, dims = 1:20)

# Update cluster identities
coembed.subset$seurat_clusters <- ifelse(!is.na(coembed.subset$seurat_clusters), coembed.subset$seurat_clusters, coembed.subset$predicted.id)

# Plotting UMAP results
p1 <- DimPlot(coembed.subset, group.by = "tech")
p4 <- DimPlot(coembed.subset, group.by = "orig.ident")
ggsave(p1 + p4, filename = "022_coembed.subset.hp0_N_UE.RNA.refdata.integrated.Rds.tech.orig.ident.pdf", width = 12, height = 6)

p2 <- DimPlot(coembed.subset, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p3 <- DimPlot(coembed.subset, split.by = "tech", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
ggsave(p2 + p3 + plot_layout(widths = c(1, 2)), filename = "023_coembed.subset.hp0_N_UE.RNA.refdata.integrated.Rds.clusters_by_tech.pdf", width = 18, height = 6)

## Mouse scRNA-seq and snATAC-seq integration

# Set working directories
setwd("newest_all_inclusive_analysis")
setwd("newest_all_inclusive_analysis")

# Load required libraries
lapply(c("Signac", "Seurat", "ggplot2", "patchwork", "readr", "tidyr", "dplyr", "tibble", "rtracklayer", "stringr"), library, character.only = TRUE)

# Load integrated mouse RNA data
transcriptome.integrated.mouse <- readRDS("newest_seurat_Mp0_Hp0_integrated/010_Hp0_Mp0_integrated_trial5.mouse.RNA_norm.Rds")
DefaultAssay(transcriptome.integrated.mouse) <- "integrated"

# Load mouse ATAC data
mp0_N_UE.ATAC <- readRDS("new_signac/atac.newUMAP.rm01091113.rm0825.rm24.Rds")
DefaultAssay(mp0_N_UE.ATAC) <- "RNA"

# Load gene translation table and create a mapping from mouse to human gene names
geneTable <- read_csv("new_seurat_Mp0_Hp0_integrated/geneTrans.txt")
mouseToHumanGeneName <- left_join(
    x = mp0_N_UE.ATAC@assays$RNA@counts@Dimnames[[1]] %>% as_tibble(),
    y = geneTable %>% dplyr::select(c(Human.Symbol, Mouse.Symbol)), by = c("value" = "Mouse.Symbol")
) %>% mutate(newnames = ifelse(is.na(Human.Symbol), value, Human.Symbol)) %>% pull(newnames)

# Function to rename genes in a Seurat object
RenameGenesSeurat <- function(obj, newnames) {
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  
  obj@assays$RNA <- RNA
  return(obj)
}

# Rename genes in the ATAC data
mp0_N_UE.ATAC <- RenameGenesSeurat(obj = mp0_N_UE.ATAC, newnames = mouseToHumanGeneName)

# Find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = transcriptome.integrated.mouse,
  query = mp0_N_UE.ATAC,
  reference.assay = "integrated", 
  query.assay = "RNA",
  reduction = 'cca',
  dims = 1:30
)

# Transfer data to get cell type predictions
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = transcriptome.integrated.mouse@active.ident,
  weight.reduction = mp0_N_UE.ATAC[["lsi"]],
  dims = 2:30
)

# Add cell type predictions to ATAC metadata
mp0_N_UE.ATAC <- AddMetaData(mp0_N_UE.ATAC, metadata = celltype.predictions)

# Update cluster identities and generate UMAP plots
Idents(object = mp0_N_UE.ATAC) <- 'common_ATAC_snn_res.1.2'
p1 <- UMAPPlot(mp0_N_UE.ATAC, label = TRUE, repel = TRUE) + NoLegend()
Idents(object = mp0_N_UE.ATAC) <- 'predicted.id'
Idents(object = mp0_N_UE.ATAC) <- factor(x = Idents(mp0_N_UE.ATAC), levels = stringr::str_sort(levels(mp0_N_UE.ATAC), numeric = TRUE))

p2 <- UMAPPlot(mp0_N_UE.ATAC, label = TRUE, repel = TRUE) + NoLegend()
p3 <- UMAPPlot(transcriptome.integrated.mouse, label = TRUE, repel = TRUE) + NoLegend()

# Save combined UMAP plot
ggsave(p1 + p2 + p3, filename = "003_mp0_N_UE.predicted.id.pdf", width = 18, height = 6)

# Identify variable features
genes.use <- VariableFeatures(transcriptome.integrated.mouse)

# Transfer imputation data
refdata.integrated <- GetAssayData(transcriptome.integrated.mouse, assay = "integrated", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata.integrated, weight.reduction = mp0_N_UE.ATAC[["lsi"]])

# Add integrated imputation to ATAC data
mp0_N_UE.ATAC[["integrated"]] <- imputation

# Reset ATAC-seq's previous clustering information
mp0_N_UE.ATAC@meta.data$seurat_clusters <- NULL
transcriptome.integrated.mouse$tech <- "RNA"
mp0_N_UE.ATAC$tech <- "ATAC"

# Merge RNA-seq and ATAC-seq data
coembed.subset <- merge(x = transcriptome.integrated.mouse, y = mp0_N_UE.ATAC)

# Re-cluster and run PCA and UMAP
DefaultAssay(coembed.subset) <- "integrated"
coembed.subset <- ScaleData(coembed.subset, features = genes.use, do.scale = TRUE)
coembed.subset <- RunPCA(coembed.subset, features = genes.use, verbose = FALSE)
coembed.subset <- RunUMAP(coembed.subset, dims = 1:20)

# Update cluster identities
coembed.subset$seurat_clusters <- ifelse(!is.na(coembed.subset$seurat_clusters), coembed.subset$seurat_clusters, coembed.subset$predicted.id)

# Plotting UMAP results
p1 <- DimPlot(coembed.subset, group.by = "tech")
p4 <- DimPlot(coembed.subset, group.by = "orig.ident")
ggsave(p1 + p4, filename = "008_coembed.subset.mp0_N_UE.RNA.refdata.integrated.Rds.tech.orig.ident.pdf", width = 12, height = 6)

p2 <- DimPlot(coembed.subset, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p3 <- DimPlot(coembed.subset, split.by = "tech", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
ggsave(p2 + p3 + plot_layout(widths = c(1, 2)), filename = "009_coembed.subset.mp0_N_UE.RNA.refdata.integrated.Rds.clusters_by_tech.pdf", width = 18, height = 6)

p5 <- DimPlot(coembed.subset, group.by = "tech")
p6 <- DimPlot(coembed.subset, group.by = "replicate")
ggsave(p5 + p6, filename = "008_coembed.subset.mp0_N_UE.RNA.refdata.integrated.Rds.tech.replicate.pdf", width = 12, height = 6)




coembed.subset.human <- readRDS("021_coembed.subset.hp0_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds")
coembed.subset.mouse <- readRDS("007_coembed.subset.mp0_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds")

coembed.subset.human$species <- "human"
coembed.subset.mouse$species <- "mouse"

coembed.subset.both <- merge(x=coembed.subset.human, y=coembed.subset.mouse)
coembed.subset.both %>% saveRDS("coembed.subset.both.Rds", compress=F)

genes.use <- readRDS("genes.use.transcriptome.integrated.human.Rds")

DefaultAssay(coembed.subset.both) <- "integrated"
coembed.subset.both <- ScaleData(coembed.subset.both, features = genes.use, do.scale = T)
coembed.subset.both <- RunPCA(coembed.subset.both, features = genes.use, verbose = FALSE)
coembed.subset.both <- RunUMAP(coembed.subset.both, dims = 1:20)

coembed.subset.both$seurat_clusters <- ifelse(!is.na(coembed.subset.both$seurat_clusters), coembed.subset.both$seurat_clusters, coembed.subset.both$predicted.id)
saveRDS(coembed.subset.both, "0070_coembed.subset.both_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds", compress=F)

coembed.subset.both <- readRDS("0070_coembed.subset.both_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds")
my_level <- c(0,12,17,10,2,16,13,19,15,1,11,7,8,18,9,14,20,6,22,4,24,5,3,21,23)
Idents(coembed.subset.both) <- factor(Idents(coembed.subset.both), levels= my_level)

# Load previously merged and processed human and mouse data
coembed.subset.human <- readRDS("021_coembed.subset.hp0_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds")
coembed.subset.mouse <- readRDS("007_coembed.subset.mp0_N_UE.RNA.refdata.integrated.Rds.UMAP.Rds")

# Add species metadata
coembed.subset.human$species <- "human"
coembed.subset.mouse$species <- "mouse"

# Merge human and mouse datasets
coembed.subset.both <- merge(x = coembed.subset.human, y = coembed.subset.mouse)

# Load genes to be used for integration
genes.use <- readRDS("genes.use.transcriptome.integrated.human.Rds")

# Set the default assay to 'integrated'
DefaultAssay(coembed.subset.both) <- "integrated"

# Scale data, run PCA and UMAP
coembed.subset.both <- ScaleData(coembed.subset.both, features = genes.use, do.scale = TRUE)
coembed.subset.both <- RunPCA(coembed.subset.both, features = genes.use, verbose = FALSE)
coembed.subset.both <- RunUMAP(coembed.subset.both, dims = 1:20)

# Update cluster identities based on RNA-seq clustering or predicted IDs from label transfer
coembed.subset.both$seurat_clusters <- ifelse(!is.na(coembed.subset.both$seurat_clusters), coembed.subset.both$seurat_clusters, coembed.subset.both$predicted.id)

# Define cluster levels and update identities
my_level <- c(0, 12, 17, 10, 2, 16, 13, 19, 15, 1, 11, 7, 8, 18, 9, 14, 20, 6, 22, 4, 24, 5, 3, 21, 23)
Idents(coembed.subset.both) <- factor(Idents(coembed.subset.both), levels = my_level)

# Plotting UMAP results
p1 <- DimPlot(coembed.subset.both, group.by = "tech")
p4 <- DimPlot(coembed.subset.both, group.by = "orig.ident")
ggsave(p1 + p4, filename = "0080_coembed.subset.both_N_UE.RNA.refdata.integrated.Rds.tech.orig.ident.pdf", width = 12, height = 6)

p2 <- DimPlot(coembed.subset.both, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p3 <- DimPlot(coembed.subset.both, split.by = "tech", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
ggsave(p2 + p3 + plot_layout(widths = c(1, 2)), filename = "0090_coembed.subset.both_N_UE.RNA.refdata.integrated.Rds.clusters_by_tech.pdf", width = 18, height = 6)

p2 <- DimPlot(coembed.subset.both, label = TRUE, repel = TRUE) + NoLegend()
p4 <- DimPlot(coembed.subset.both, split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
ggsave(p2 + p4 + plot_layout(widths = c(1, 2)), filename = "0100_coembed.subset.both_N_UE.RNA.refdata.integrated.Rds.clusters_by_species.pdf", width = 18, height = 7)

# Display the combined object
print(coembed.subset.both)

