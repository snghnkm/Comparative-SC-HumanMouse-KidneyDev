# 002_Mouse_p0_ATAC.r

# Load required libraries
lapply(c("Signac", "Seurat", "dplyr", "GenomeInfoDb", "EnsDb.Mmusculus.v79", "ggplot2", "patchwork", "readr"), library, character.only = TRUE)
library(Signac)          # For single-cell chromatin accessibility analysis
library(chromVARmotifs)  # For chromatin variability analysis
data("mouse_pwms_v1")    # Load mouse position weight matrices
library(JASPAR2020)      # For accessing the JASPAR database of transcription factor binding profiles

# Creating a unified peak set from multiple scATAC-seq datasets

# Define file paths for peak sets
m1_peak <- "scATAC/Matac1d_m1/outs/peaks.bed"
m2_peak <- "scATAC/Matac2b_m2/outs/peaks.bed"
f1_peak <- "scATAC/Matac3c_f1/outs/peaks.bed"
f2_peak <- "scATAC/Matac4_f2/outs/peaks.bed"

# Read in peak sets and convert to genomic ranges
gr.m1 <- read.table(file = m1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.m2 <- read.table(file = m2_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.f1 <- read.table(file = f1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.f2 <- read.table(file = f2_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

# Create a unified set of peaks to quantify in each dataset by reducing overlapping peaks
combined.peaks <- reduce(x = c(gr.m1, gr.m2, gr.f1, gr.f2))

# Filter out peaks based on length (retain peaks between 20 and 10,000 base pairs)
peakwidths <- width(combined.peaks)
combined.peaks.filtered <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks.filtered

# Save the filtered combined peaks to a file
saveRDS(combined.peaks.filtered, "combined.peaks.filtered.signac.1.0.0.Rds", compress = FALSE)

# Creating an intersectioned peak set (common peaks present in all datasets)
intersection.peaks <- Reduce(subsetByOverlaps, list(gr.m1, gr.m2, gr.f1, gr.f2))
intersection.peaks

# Filter out peaks based on length (retain peaks between 20 and 10,000 base pairs)
peakwidths <- width(intersection.peaks)
intersection.peaks.filtered <- intersection.peaks[peakwidths < 10000 & peakwidths > 20]
intersection.peaks.filtered

# Save the filtered intersection peaks to a file
saveRDS(intersection.peaks.filtered, "intersection.peaks.filtered.signac.1.0.0.Rds", compress = FALSE)

# Load the saved peak sets
combined.peaks.filtered <- readRDS("combined.peaks.filtered.signac.1.0.0.Rds")
intersection.peaks.filtered <- readRDS("intersection.peaks.filtered.signac.1.0.0.Rds")

# Define file paths for fragments and metadata
m1_fragment <- "scATAC/Matac1d_m1/outs/fragments.tsv.gz"
m2_fragment <- "scATAC/Matac2b_m2/outs/fragments.tsv.gz"
f1_fragment <- "scATAC/Matac3c_f1/outs/fragments.tsv.gz"
f2_fragment <- "scATAC/Matac4_f2/outs/fragments.tsv.gz"

m1_metadata <- "scATAC/Matac1d_m1/outs/singlecell.csv"
m2_metadata <- "scATAC/Matac2b_m2/outs/singlecell.csv"
f1_metadata <- "scATAC/Matac3c_f1/outs/singlecell.csv"
f2_metadata <- "scATAC/Matac4_f2/outs/singlecell.csv"

# Create a custom function to generate Seurat objects for scATAC-seq data
createSignacObject <- function(metadata_path, fragment_path, combined.peaks, intersection.peaks) {
    # Read in metadata
    md <- read.table(
        file = metadata_path,
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ] # Remove the first row (column names)
    print(dim(md))
    
    # Perform initial filtering of low count cells (retain cells with more than 500 fragments)
    md.filtered <- md[md$passed_filters > 500, ]
    print(dim(md.filtered))

    # Create fragment objects
    frags <- CreateFragmentObject(
        path = fragment_path,
        cells = rownames(md.filtered)
    )

    # Quantify peaks in each dataset using the combined and intersection peak sets
    counts <- FeatureMatrix(
        fragments = frags,
        features = combined.peaks,
        cells = rownames(md.filtered)
    )
    counts_common <- FeatureMatrix(
        fragments = frags,
        features = intersection.peaks,
        cells = rownames(md.filtered)
    )

    # Create Chromatin Assays and Seurat objects
    assay <- CreateChromatinAssay(counts, fragments = frags)
    seuratObj <- CreateSeuratObject(
        assay, assay = "ATAC", meta.data = md.filtered
    )
    assay_common <- CreateChromatinAssay(counts_common, fragments = frags)
    seuratObj[["common_ATAC"]] <- assay_common

    return(seuratObj)
}


# Create Signac objects for each dataset and label them
atac.M1 <- createSignacObject(m1_metadata, m1_fragment, combined.peaks.filtered, intersection.peaks.filtered)
atac.M1$dataset <- 'M1'  # Label dataset as 'M1'

atac.M2 <- createSignacObject(m2_metadata, m2_fragment, combined.peaks.filtered, intersection.peaks.filtered)
atac.M2$dataset <- 'M2'  # Label dataset as 'M2'

atac.F1 <- createSignacObject(f1_metadata, f1_fragment, combined.peaks.filtered, intersection.peaks.filtered)
atac.F1$dataset <- 'F1'  # Label dataset as 'F1'

atac.F2 <- createSignacObject(f2_metadata, f2_fragment, combined.peaks.filtered, intersection.peaks.filtered)
atac.F2$dataset <- 'F2'  # Label dataset as 'F2'


# Annotation, QC
# Annotation

annotateUCSCstyle <- function(seuratObj) {
    # Extract gene annotations from EnsDb
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

    # Change to UCSC style since the data was mapped to hg19
    seqlevelsStyle(annotations) <- 'UCSC'
    genome(annotations) <- "mm10"

    # Add the gene information to the object
    Annotation(seuratObj) <- annotations
    return(seuratObj)
}

# QC (Quality Control)
defaultQC <- function(seuratObj) {
    # Compute nucleosome signal score per cell
    seuratObj <- NucleosomeSignal(object = seuratObj)

    # Compute TSS enrichment score per cell
    seuratObj <- TSSEnrichment(object = seuratObj, fast = TRUE)

    # Add blacklist ratio and fraction of reads in peaks
    seuratObj$pct_reads_in_peaks <- seuratObj$peak_region_fragments / seuratObj$passed_filters * 100
    seuratObj$blacklist_ratio <- seuratObj$blacklist_region_fragments / seuratObj$peak_region_fragments

    return(seuratObj)
}

applyingQC <- function(seuratObj) {
    # Subset the Seurat object based on QC metrics
    seuratObj <- subset(
        x = seuratObj,
        subset = peak_region_fragments > 2000 &  # Minimum peak region fragments
            peak_region_fragments < 40000 &  # Maximum peak region fragments
            pct_reads_in_peaks > 25 &  # Minimum percentage of reads in peaks
            blacklist_ratio < 0.05 &  # Maximum blacklist ratio
            nucleosome_signal < 2 &  # Maximum nucleosome signal
            TSS.enrichment > 2  # Minimum TSS enrichment
    )
    return(seuratObj)
}

# Merge Objects
# Load the pre-saved Signac objects
atac.M1 <- readRDS("atac.M1.signac.1.0.0.new.Rds")
atac.M2 <- readRDS("atac.M2.signac.1.0.0.new.Rds")
atac.F1 <- readRDS("atac.F1.signac.1.0.0.new.Rds")
atac.F2 <- readRDS("atac.F2.signac.1.0.0.new.Rds")

# Merge all datasets, adding a cell ID to ensure unique cell names
combined <- merge(
  x = atac.M1,
  y = list(atac.M2, atac.F1, atac.F2),
  add.cell.ids = c("M1", "M2", "F1", "F2")
)
combined <- RunTFIDF(combined)  # Run Term Frequency-Inverse Document Frequency (TF-IDF) on combined object
combined <- FindTopFeatures(combined, min.cutoff = "q0")  # Find top features with no minimum cutoff
combined <- RunSVD(combined)  # Run Singular Value Decomposition (SVD)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')  # Run UMAP for dimensionality reduction

# Save the combined object
saveRDS(combined, "combined.ATAC.defaultQC.signac.1.0.0.QC.TFIDF.SVD.UMAP.Rds", compress = FALSE)

# Gene Activity

# Since our object has F1_, F2_, M1_, M2_ suffix
# Use this combined fragment with modified barcodes
fragment.path <- "atac.combined.fragments.tsv.gz"

# Create Fragment object
fragments <- CreateFragmentObject(
  path = fragment.path,
  cells = colnames(combined),
  validate.fragments = FALSE
)

# Extract gene coordinates for protein-coding genes
gene.coords <- genes(EnsDb.Mmusculus.v79, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'  # Convert sequence levels to UCSC style
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')  # Keep standard chromosomes
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)  # Extend coordinates to include promoters

# Build a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragments,
  features = genebodyandpromoter.coords,
  cells = colnames(combined),
  verbose = TRUE
)

# Convert row names from chromosomal coordinates to gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
gene.activities <- gene.activities[rownames(gene.activities) != "",]

# Add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

# Save the final combined object with gene activities
saveRDS(combined, "combined.ATAC.QC.Anno.signac.1.0.0.TFIDF.q5.SVD.n50.UMAP.GA.Rds", compress = FALSE)

# Set the default assay to 'RNA' and perform data scaling
DefaultAssay(combined) <- 'RNA'
combined <- ScaleData(combined, verbose = TRUE)

# Find variable features using the 'vst' method
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000, verbose = TRUE)

# Run Principal Component Analysis (PCA)
combined <- RunPCA(combined, npcs = 30, verbose = TRUE)

# Find neighbors and clusters based on PCA
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20, k.param = 15)
combined <- FindClusters(combined, resolution = 0.8)

# Run UMAP for dimensionality reduction
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)

# Load the combined object from RDS file
combined <- readRDS("combined.ATAC.QC.Anno.signac.1.0.0.TFIDF.q5.SVD.n50.UMAP.GA.Rds")

# Set the default assay to 'RNA'
DefaultAssay(combined) <- 'RNA'

# Generate and save a feature plot for Six2 and Cited1
ggsave(FeaturePlot(combined, features = c("Six2", "Cited1"), label = TRUE), filename = "combined.six2.cited1.pdf")

# Define a function to apply quality control filters
applyingQC <- function(seuratObj) {
    seuratObj <- subset(
        x = seuratObj,
        subset = peak_region_fragments > 2000 &
            peak_region_fragments < 50000 &  # Adjusted maximum peak region fragments
            pct_reads_in_peaks > 30 &  # Adjusted minimum percentage of reads in peaks
            blacklist_ratio < 0.05 &
            nucleosome_signal < 2 &
            TSS.enrichment > 2
    )
    return(seuratObj)
}

# Apply QC filters to each dataset
atac.M1.QC.filtered <- applyingQC(atac.M1)
atac.M2.QC.filtered <- applyingQC(atac.M2)
atac.F1.QC.filtered <- applyingQC(atac.F1)
atac.F2.QC.filtered <- applyingQC(atac.F2)

# Merge QC filtered datasets and ensure unique cell names
unintegrated.anno.QC.filtered <- merge(
  x = atac.M1.QC.filtered,
  y = list(atac.M2.QC.filtered, atac.F1.QC.filtered, atac.F2.QC.filtered),
  add.cell.ids = c("M1", "M2", "F1", "F2")
)

# Set the default assay to 'common_ATAC' and perform further analysis
DefaultAssay(unintegrated.anno.QC.filtered) <- "common_ATAC"
unintegrated.anno.QC.filtered <- RunTFIDF(unintegrated.anno.QC.filtered)
unintegrated.anno.QC.filtered <- FindTopFeatures(unintegrated.anno.QC.filtered, min.cutoff = 50)
unintegrated.anno.QC.filtered <- RunSVD(unintegrated.anno.QC.filtered)
unintegrated.anno.QC.filtered <- RunUMAP(unintegrated.anno.QC.filtered, reduction = 'lsi', dims = 2:30)

# Generate and save a UMAP plot grouped by dataset
p1 <- DimPlot(unintegrated.anno.QC.filtered, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("unintegrated.anno.QC.filtered")
ggsave(p1, filename = "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.pdf")

# Calculate gene activities
gene.activities <- GeneActivity(unintegrated.anno.QC.filtered)

# Add the gene activity matrix to the Seurat object as a new assay and normalize it
unintegrated.anno.QC.filtered[['RNA']] <- CreateAssayObject(counts = gene.activities)
unintegrated.anno.QC.filtered <- NormalizeData(
  object = unintegrated.anno.QC.filtered,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(unintegrated.anno.QC.filtered$nCount_RNA)
)

# Save the updated object
saveRDS(unintegrated.anno.QC.filtered, "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Rds", compress = FALSE)

# Perform clustering on the 'common_ATAC' assay
DefaultAssay(unintegrated.anno.QC.filtered) <- "common_ATAC"
unintegrated.anno.QC.filtered <- unintegrated.anno.QC.filtered %>%
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
    FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Reload the updated object and perform clustering again
unintegrated.anno.QC.filtered <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Rds")

DefaultAssay(unintegrated.anno.QC.filtered) <- "common_ATAC" unintegrated.anno.QC.filtered <- unintegrated.anno.QC.filtered %>% FindNeighbors(reduction = 'lsi', dims = 2:30) %>% FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Save the final clustered object
saveRDS(unintegrated.anno.QC.filtered, "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.Rds", compress = FALSE)

# Run chromVAR analysis on unintegrated data
unintegrated.anno.QC.filtered <- RunChromVAR(
  object = unintegrated.anno.QC.filtered,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = mouse_pwms_v1
)

# Load JASPAR motif data
pfm.JASPAR2020 <- readRDS("JASPAR2020.9606.pfm.Rds")

# Load the processed Seurat object
unintegrated.anno.QC.filtered <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.added.Rds")

# Copy the common_ATAC assay to a new assay common_ATAC2
unintegrated.anno.QC.filtered[['common_ATAC2']] <- unintegrated.anno.QC.filtered[['common_ATAC']]

# Create motif matrix using JASPAR2020 motif data
motif.matrix <- CreateMotifMatrix(
  features = Signac::StringToGRanges(regions = rownames(unintegrated.anno.QC.filtered[['common_ATAC2']]@meta.features), sep = c(":", "-")),
  pwm = pfm.JASPAR2020,
  genome = 'mm10',
  use.counts = FALSE,
  sep = c("-", "-")
)

# Find motif positions in the genome
motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm.JASPAR2020, 
  subject = Signac::StringToGRanges(regions = rownames(unintegrated.anno.QC.filtered[['common_ATAC2']]@meta.features), sep = c(":", "-")), 
  out = 'positions', 
  genome = 'mm10'
)

# Create Motif object and add it to the common_ATAC2 assay
motif <- CreateMotifObject(data = motif.matrix, positions = motif.positions, pwm = pfm.JASPAR2020)
unintegrated.anno.QC.filtered[['common_ATAC2']] <- SetAssayData(object = unintegrated.anno.QC.filtered[['common_ATAC2']], slot = 'motifs', new.data = motif)

# Calculate region statistics
library(BSgenome.Mmusculus.UCSC.mm10)
unintegrated.anno.QC.filtered <- RegionStats(
  object = unintegrated.anno.QC.filtered,
  assay = 'common_ATAC2',
  genome = BSgenome.Mmusculus.UCSC.mm10,
  verbose = TRUE
)

# Set the default assay to common_ATAC2
DefaultAssay(unintegrated.anno.QC.filtered) <- "common_ATAC2"

# Run chromVAR analysis again on the updated object
unintegrated.anno.QC.filtered.ChromVAR2 <- RunChromVAR(
  object = unintegrated.anno.QC.filtered[['common_ATAC2']],
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# Save the updated object with chromVAR results
saveRDS(unintegrated.anno.QC.filtered.ChromVAR2, "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.added.JASPARmotif.RegionStats.chromVAR.added.Rds", compress = FALSE)

# Load the cisBP chromVAR assay and the Seurat object
signac <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.Rds")
unintegrated.anno.QC.filtered <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.Rds")

# Add cisBP chromVAR assay to the Seurat object
unintegrated.anno.QC.filtered[['common_ATAC_cisBP_chromVAR']] <- signac

# Load the Seurat object with JASPAR motifs
signac <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.added.JASPARmotif.RegionStats.chromVAR.added.Rds")

# Add JASPAR chromVAR assay to the Seurat object
unintegrated.anno.QC.filtered <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.added.Rds")
unintegrated.anno.QC.filtered[['common_ATAC_JASPAR_chromVAR']] <- signac

# Load the final Seurat object with both cisBP and JASPAR motifs
unintegrated.anno.QC.filtered <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.both.chromVAR.added.final.Rds")
signac <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.chromVAR.added.JASPARmotif.all.RegionStats.chromVAR.added.Rds")
unintegrated.anno.QC.filtered[['common_ATAC_allJASPAR_chromVAR']] <- signac # all JASPAR motifs

saveRDS(unintegrated.anno.QC.filtered, "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.triple.chromVAR.added.final.of.final.Rds", compress=F)

# Mouse N-UE subclustering

# Load the processed ATAC-seq data
atac <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.triple.chromVAR.added.final.of.final.Rds")

# Subset the ATAC-seq data for specific clusters
atac.subset <- subset(atac, idents = c(0, 11, 18, 24, 19, 13, 32, 12, 16, 3, 28, 26, 5, 9, 15, 14, 20, 30, 31))

# Set the default assay to 'common_ATAC'
DefaultAssay(atac.subset) <- "common_ATAC"
atac.subset <- RunTFIDF(atac.subset)
atac.subset <- FindTopFeatures(atac.subset, min.cutoff = 50)
atac.subset <- RunSVD(atac.subset)
atac.subset <- FindNeighbors(atac.subset, reduction = 'lsi', dims = 2:30)
atac.subset <- FindClusters(atac.subset, verbose = TRUE, algorithm = 3, resolution = 1.2)

# Run UMAP for dimensionality reduction
atac.subset <- RunUMAP(atac.subset, reduction = 'lsi', dims = 2:30)


saveRDS(atac.subset, "atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.triple.chromVAR.added.final.of.final.N_UE.newUMAP.Rds", compress=F)
atac.newUMAP <- readRDS("atac.signac.1.0.0.new.unintegrated.anno.QC.filtered.TFIDF.SVD.UMAP.GeneActivity.Norm.Clustered.cisBPmotif.RegionStats.triple.chromVAR.added.final.of.final.N_UE.newUMAP.Rds")

# Remove specific clusters from the new UMAP data
atac.newUMAP.rm01091113 <- subset(atac.newUMAP, idents = c(1, 9, 11, 13), invert = TRUE)

# Set the default assay to 'common_ATAC' and run further analysis
DefaultAssay(atac.newUMAP.rm01091113) <- "common_ATAC"
atac.newUMAP.rm01091113 <- atac.newUMAP.rm01091113 %>% 
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff = "q0") %>% 
    RunSVD() %>% 
    RunUMAP(reduction = 'lsi', dims = 2:30) %>% 
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>% 
    FindClusters(algorithm = 3, resolution = 1.2)

# Further subset the data to remove additional clusters
atac.newUMAP.rm01091113.rm0825 <- subset(atac.newUMAP.rm01091113, idents = c(8, 25), invert = TRUE)

# Set the default assay and run the analysis again
DefaultAssay(atac.newUMAP.rm01091113.rm0825) <- "common_ATAC"
atac.newUMAP.rm01091113.rm0825 <- atac.newUMAP.rm01091113.rm0825 %>% 
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff = "q0") %>% 
    RunSVD() %>% 
    RunUMAP(reduction = 'lsi', dims = 2:30) %>% 
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>% 
    FindClusters(algorithm = 3, resolution = 1.2)

# Further subset the data to remove another cluster
atac.newUMAP.rm01091113.rm0825.rm24 <- subset(atac.newUMAP.rm01091113.rm0825, idents = c(24), invert = TRUE)

# Set the default assay and run the analysis again
DefaultAssay(atac.newUMAP.rm01091113.rm0825.rm24) <- "common_ATAC"
atac.newUMAP.rm01091113.rm0825.rm24 <- atac.newUMAP.rm01091113.rm0825.rm24 %>% 
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff = "q0") %>% 
    RunSVD() %>% 
    RunUMAP(reduction = 'lsi', dims = 2:30) %>% 
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>% 
    FindClusters(algorithm = 3, resolution = 1.2)

# Save the UMAP plot to a PDF file
ggsave(DimPlot(atac.newUMAP.rm01091113.rm0825.rm24, label = TRUE, repel = TRUE), filename = "atac.newUMAP.rm01091113.rm0825.rm24.UMAP.pdf")

# Save violin plots of quality control metrics to a PDF file
ggsave(
    VlnPlot(atac.newUMAP.rm01091113.rm0825.rm24, c("nCount_ATAC", "nFeature_ATAC", "nCount_common_ATAC", "nFeature_common_ATAC", "duplicate", "mitochondrial", "TSS_fragments", "peak_region_fragments", "on_target_fragments", "enhancer_region_fragments", "blacklist_region_fragments", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_percentile", "nCount_RNA", "nFeature_RNA"), pt.size = 0), 
    filename = "atac.newUMAP.rm01091113.rm0825.rm24.QC.pdf", 
    width = 20, height = 16, limitsize = FALSE
)