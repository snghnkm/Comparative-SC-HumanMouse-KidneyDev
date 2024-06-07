# 004_Human_fetal_ATAC.r

# Load necessary libraries
lapply(c("Signac", "Seurat", "dplyr", "GenomeInfoDb", "EnsDb.Hsapiens.v86", "ggplot2", "patchwork", "readr"), library, character.only = TRUE)

# Creating a common "UNIFIED" peak set

# File paths for peak files
Hp0_BGI_M1_peak = "snATACseq/10_6_M_snATAC/outs/peaks.bed"
Hp0_BGI_M2_peak = "snATACseq/11_1_M_snATAC/outs/peaks.bed"
Hp0_BGI_M3_peak = "snATACseq/11_6_M_snATAC/outs/peaks.bed"
Hp0_BGI_F1_peak = "snATACseq/11_3_F_snATAC/outs/peaks.bed"

Hp0_NOV_F1_peak = "snATACseq/11_2_F_snATAC/outs/peaks.bed"
Hp0_NOV_M1_peak = "snATACseq/17_1_M_snATAC/outs/peaks.bed"
Hp0_NOV_M2_peak = "snATACseq/17_5_M_snATAC/outs/peaks.bed"
Hp0_NOV_F2_peak = "snATACseq/17_1_F_snATAC/outs/peaks.bed"
Hp0_NOV_F3_peak = "snATACseq/17_6_F_snATAC/outs/peaks.bed"

# Read in peak sets and convert to genomic ranges
gr.BGI_M1_peak <- read.table(file = Hp0_BGI_M1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.BGI_M2_peak <- read.table(file = Hp0_BGI_M2_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.BGI_M3_peak <- read.table(file = Hp0_BGI_M3_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.BGI_F1_peak <- read.table(file = Hp0_BGI_F1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

gr.NOV_F1_peak <- read.table(file = Hp0_NOV_F1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.NOV_M1_peak <- read.table(file = Hp0_NOV_M1_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.NOV_M2_peak <- read.table(file = Hp0_NOV_M2_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.NOV_F2_peak <- read.table(file = Hp0_NOV_F2_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
gr.NOV_F3_peak <- read.table(file = Hp0_NOV_F3_peak, col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()

# Create a unified set of peaks to quantify in each dataset
BGI.combined.peaks <- reduce(x = c(gr.BGI_M1_peak, gr.BGI_M2_peak, gr.BGI_M3_peak, gr.BGI_F1_peak))
NOV.combined.peaks <- reduce(x = c(gr.NOV_F1_peak, gr.NOV_M1_peak, gr.NOV_M2_peak, gr.NOV_F2_peak, gr.NOV_F3_peak))

# Filter out bad peaks based on length
BGI.peakwidths <- width(BGI.combined.peaks)
NOV.peakwidths <- width(NOV.combined.peaks)

# Arbitrary peak width filters
BGI.combined.peaks.filtered <- BGI.combined.peaks[BGI.peakwidths < 10000 & BGI.peakwidths > 20]
NOV.combined.peaks.filtered <- NOV.combined.peaks[NOV.peakwidths < 10000 & NOV.peakwidths > 20]

# Creating a common "INTERSECTIONED" peak set

BGI.intersection.peaks <- Reduce(subsetByOverlaps, list(gr.BGI_M1_peak, gr.BGI_M2_peak, gr.BGI_M3_peak, gr.BGI_F1_peak))
NOV.intersection.peaks <- Reduce(subsetByOverlaps, list(gr.NOV_F1_peak, gr.NOV_M1_peak, gr.NOV_M2_peak, gr.NOV_F2_peak, gr.NOV_F3_peak))

# Filter out bad peaks based on length
BGI.peakwidths <- width(BGI.intersection.peaks)
NOV.peakwidths <- width(NOV.intersection.peaks)

# Arbitrary peak width filters
BGI.intersection.peaks.filtered <- BGI.intersection.peaks[BGI.peakwidths < 10000 & BGI.peakwidths > 20]
NOV.intersection.peaks.filtered <- NOV.intersection.peaks[NOV.peakwidths < 10000 & NOV.peakwidths > 20]

# Create the objects

# Define file paths for fragment and metadata files
Hp0_BGI_M1_fragment = "snATACseq/10_6_M_snATAC/outs/fragments.tsv.gz"
Hp0_BGI_M2_fragment = "snATACseq/11_1_M_snATAC/outs/fragments.tsv.gz"
Hp0_BGI_M3_fragment = "snATACseq/11_6_M_snATAC/outs/fragments.tsv.gz"
Hp0_BGI_F1_fragment = "snATACseq/11_3_F_snATAC/outs/fragments.tsv.gz"

Hp0_NOV_F1_fragment = "snATACseq/11_2_F_snATAC/outs/fragments.tsv.gz"
Hp0_NOV_M1_fragment = "snATACseq/17_1_M_snATAC/outs/fragments.tsv.gz"
Hp0_NOV_M2_fragment = "snATACseq/17_5_M_snATAC/outs/fragments.tsv.gz"
Hp0_NOV_F2_fragment = "snATACseq/17_1_F_snATAC/outs/fragments.tsv.gz"
Hp0_NOV_F3_fragment = "snATACseq/17_6_F_snATAC/outs/fragments.tsv.gz"

Hp0_BGI_M1_metadata = "snATACseq/10_6_M_snATAC/outs/singlecell.csv"
Hp0_BGI_M2_metadata = "snATACseq/11_1_M_snATAC/outs/singlecell.csv"
Hp0_BGI_M3_metadata = "snATACseq/11_6_M_snATAC/outs/singlecell.csv"
Hp0_BGI_F1_metadata = "snATACseq/11_3_F_snATAC/outs/singlecell.csv"

Hp0_NOV_F1_metadata = "snATACseq/11_2_F_snATAC/outs/singlecell.csv"
Hp0_NOV_M1_metadata = "snATACseq/17_1_M_snATAC/outs/singlecell.csv"
Hp0_NOV_M2_metadata = "snATACseq/17_5_M_snATAC/outs/singlecell.csv"
Hp0_NOV_F2_metadata = "snATACseq/17_1_F_snATAC/outs/singlecell.csv"
Hp0_NOV_F3_metadata = "snATACseq/17_6_F_snATAC/outs/singlecell.csv"

# Custom function to create Signac objects
createSignacObject <- function(metadata_path, fragment_path, combined.peaks, intersection.peaks) {
    md <- read.table(
        file = metadata_path,
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ] # remove the first row
    print(dim(md))
    
    # Perform an initial filtering of low count cells (500 is arbitrary)
    md.filtered <- md[md$passed_filters > 500, ]
    print(dim(md.filtered))

    # Create fragment objects
    frags <- CreateFragmentObject(
        path = fragment_path,
        cells = rownames(md.filtered)
    )

    # Quantify peaks in each dataset
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

    assay <- CreateChromatinAssay(counts, fragments = frags)
    seuratObj <- CreateSeuratObject(
        assay, assay = "ATAC", meta.data = md.filtered
    )
    assay_common <- CreateChromatinAssay(counts_common, fragments = frags)
    seuratObj[["common_ATAC"]] <- assay_common

    return(seuratObj)
}

# Create Signac objects for each dataset
Hp0_BGI_M1 <- createSignacObject(Hp0_BGI_M1_metadata, Hp0_BGI_M1_fragment, BGI.combined.peaks.filtered, BGI.intersection.peaks.filtered)
Hp0_BGI_M2 <- createSignacObject(Hp0_BGI_M2_metadata, Hp0_BGI_M2_fragment, BGI.combined.peaks.filtered, BGI.intersection.peaks.filtered)
Hp0_BGI_M3 <- createSignacObject(Hp0_BGI_M3_metadata, Hp0_BGI_M3_fragment, BGI.combined.peaks.filtered, BGI.intersection.peaks.filtered)
Hp0_BGI_F1 <- createSignacObject(Hp0_BGI_F1_metadata, Hp0_BGI_F1_fragment, BGI.combined.peaks.filtered, BGI.intersection.peaks.filtered)

Hp0_BGI_M1$dataset <- 'BGI_M1'
Hp0_BGI_M2$dataset <- 'BGI_M2'
Hp0_BGI_M3$dataset <- 'BGI_M3'
Hp0_BGI_F1$dataset <- 'BGI_F1'

Hp0_NOV_F1 <- createSignacObject(Hp0_NOV_F1_metadata, Hp0_NOV_F1_fragment, NOV.combined.peaks.filtered, NOV.intersection.peaks.filtered)
Hp0_NOV_F1$dataset <- 'NOV_F1'

Hp0_NOV_M1 <- createSignacObject(Hp0_NOV_M1_metadata, Hp0_NOV_M1_fragment, NOV.combined.peaks.filtered, NOV.intersection.peaks.filtered)
Hp0_NOV_M1$dataset <- 'NOV_M1'

Hp0_NOV_M2 <- createSignacObject(Hp0_NOV_M2_metadata, Hp0_NOV_M2_fragment, NOV.combined.peaks.filtered, NOV.intersection.peaks.filtered)
Hp0_NOV_M2$dataset <- 'NOV_M2'

Hp0_NOV_F2 <- createSignacObject(Hp0_NOV_F2_metadata, Hp0_NOV_F2_fragment, NOV.combined.peaks.filtered, NOV.intersection.peaks.filtered)
Hp0_NOV_F2$dataset <- 'NOV_F2'

Hp0_NOV_F3 <- createSignacObject(Hp0_NOV_F3_metadata, Hp0_NOV_F3_fragment, NOV.combined.peaks.filtered, NOV.intersection.peaks.filtered)
Hp0_NOV_F3$dataset <- 'NOV_F3'

# Annotation and QC functions

annotateUCSCstyle <- function(seuratObj) {
    # Extract gene annotations from EnsDb
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

    # Change to UCSC style since the data was mapped
    seqlevelsStyle(annotations) <- 'UCSC'
    genome(annotations) <- "hg38"

    # Add the gene information to the object
    Annotation(seuratObj) <- annotations
    return(seuratObj)
}

defaultQC <- function(seuratObj) {
    # Compute TSS enrichment score per cell
    seuratObj <- TSSEnrichment(seuratObj, fast = FALSE, verbose = TRUE, n = 10000)

    # Compute nucleosome signal score per cell
    seuratObj <- NucleosomeSignal(object = seuratObj)

    # Add blacklist ratio and fraction of reads in peaks
    seuratObj$pct_reads_in_peaks <- seuratObj$peak_region_fragments / seuratObj$passed_filters * 100
    seuratObj$blacklist_ratio <- seuratObj$blacklist_region_fragments / seuratObj$peak_region_fragments

    return(seuratObj)
}

applyingQC <- function(seuratObj) {
    seuratObj <- subset(
        x = seuratObj,
        subset = peak_region_fragments > 2000 & peak_region_fragments < 50000 & pct_reads_in_peaks > 30 & blacklist_ratio < 0.05 & nucleosome_signal < 2 & TSS.enrichment > 2
    )
    return(seuratObj)
}

# Apply annotation and QC to each dataset
DefaultAssay(Hp0_BGI_M1) <- "common_ATAC"
Hp0_BGI_M1 <- annotateUCSCstyle(Hp0_BGI_M1)
Hp0_BGI_M1 <- defaultQC(Hp0_BGI_M1)
Hp0_BGI_M1.QC.filtered <- applyingQC(Hp0_BGI_M1)

DefaultAssay(Hp0_BGI_M2) <- "common_ATAC"
Hp0_BGI_M2.QC.filtered <- annotateUCSCstyle(Hp0_BGI_M2) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_BGI_M3) <- "common_ATAC"
Hp0_BGI_M3.QC.filtered <- annotateUCSCstyle(Hp0_BGI_M3) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_BGI_F1) <- "common_ATAC"
Hp0_BGI_F1.QC.filtered <- annotateUCSCstyle(Hp0_BGI_F1) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_NOV_F1) <- "common_ATAC"
Hp0_NOV_F1.QC.filtered <- annotateUCSCstyle(Hp0_NOV_F1) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_NOV_M1) <- "common_ATAC"
Hp0_NOV_M1.QC.filtered <- annotateUCSCstyle(Hp0_NOV_M1) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_NOV_M2) <- "common_ATAC"
Hp0_NOV_M2.QC.filtered <- annotateUCSCstyle(Hp0_NOV_M2) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_NOV_F2) <- "common_ATAC"
Hp0_NOV_F2.QC.filtered <- annotateUCSCstyle(Hp0_NOV_F2) %>% defaultQC() %>% applyingQC()
DefaultAssay(Hp0_NOV_F3) <- "common_ATAC"
Hp0_NOV_F3.QC.filtered <- annotateUCSCstyle(Hp0_NOV_F3) %>% defaultQC() %>% applyingQC()


# BGI Dataset Processing

# Merge BGI datasets
Hp0_BGI_filtered_four <- merge(x = Hp0_BGI_M1,
      y = list(Hp0_BGI_M2, Hp0_BGI_M3, Hp0_BGI_F1),
      add.cell.ids = c("BGI_M1", "BGI_M2", "BGI_M3", "BGI_F1"))

# Set the default assay to 'common_ATAC'
DefaultAssay(Hp0_BGI_filtered_four) <- "common_ATAC"

# Run TF-IDF normalization
Hp0_BGI_filtered_four <- RunTFIDF(Hp0_BGI_filtered_four)

# Find top features with a minimum cutoff of 50
Hp0_BGI_filtered_four <- FindTopFeatures(Hp0_BGI_filtered_four, min.cutoff = 50)

# Perform Singular Value Decomposition (SVD)
Hp0_BGI_filtered_four <- RunSVD(Hp0_BGI_filtered_four)

# Run UMAP for dimensionality reduction
Hp0_BGI_filtered_four <- RunUMAP(Hp0_BGI_filtered_four, reduction = 'lsi', dims = 2:30)

# Plot UMAP
p1 <- DimPlot(Hp0_BGI_filtered_four, group.by = 'dataset', pt.size = 0.1)
ggsave(p1, filename = "022_Hp0_BGI_filtered_four.UMAP.pdf")

# Find neighbors and clusters
Hp0_BGI_filtered_four <- Hp0_BGI_filtered_four %>%
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
    FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Plot UMAP with clusters
ggsave(UMAPPlot(Hp0_BGI_filtered_four, repel = TRUE, label = TRUE), filename = "026_Hp0_BGI_filtered_four.UMAP.res1.6.pdf")

# Generate violin plots for QC metrics
ggsave(
    VlnPlot(Hp0_BGI_filtered_four, c(
        "nCount_ATAC", "nFeature_ATAC", "nCount_common_ATAC", "nFeature_common_ATAC", "duplicate", "mitochondrial", "TSS_fragments", "peak_region_fragments", "on_target_fragments", "enhancer_region_fragments", "blacklist_region_fragments", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_percentile", "nCount_RNA", "nFeature_RNA"),
    pt.size = 0),
    filename = "027_Hp0_BGI_filtered_four.QC.pdf",
    width = 20,
    height = 16,
    limitsize = FALSE
)

# NOV Dataset Processing

# Merge NOV datasets
Hp0_NOV_filtered_five <- merge(x = Hp0_NOV_F1,
      y = list(Hp0_NOV_M1, Hp0_NOV_M2, Hp0_NOV_F2, Hp0_NOV_F3),
      add.cell.ids = c("NOV_F1", "NOV_M1", "NOV_M2", "NOV_F2", "NOV_F3"))

# Set the default assay to 'common_ATAC'
DefaultAssay(Hp0_NOV_filtered_five) <- "common_ATAC"

# Run TF-IDF normalization
Hp0_NOV_filtered_five <- RunTFIDF(Hp0_NOV_filtered_five)

# Find top features with a minimum cutoff of 50
Hp0_NOV_filtered_five <- FindTopFeatures(Hp0_NOV_filtered_five, min.cutoff = 50)

# Perform Singular Value Decomposition (SVD)
Hp0_NOV_filtered_five <- RunSVD(Hp0_NOV_filtered_five)

# Run UMAP for dimensionality reduction
Hp0_NOV_filtered_five <- RunUMAP(Hp0_NOV_filtered_five, reduction = 'lsi', dims = 2:30)

# Plot UMAP
p1 <- DimPlot(Hp0_NOV_filtered_five, group.by = 'dataset', pt.size = 0.1)
ggsave(p1, filename = "024_Hp0_NOV_filtered_five.UMAP.pdf")

# Find neighbors and clusters
Hp0_NOV_filtered_five <- Hp0_NOV_filtered_five %>%
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>%
    FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.2)

# Plot UMAP with clusters
ggsave(UMAPPlot(Hp0_NOV_filtered_five, repel = TRUE, label = TRUE), filename = "029_Hp0_NOV_filtered_five.UMAP.res1.2.pdf")

# Generate violin plots for QC metrics
ggsave(
    VlnPlot(Hp0_NOV_filtered_five, c(
        "nCount_ATAC", "nFeature_ATAC", "nCount_common_ATAC", "nFeature_common_ATAC", "duplicate", "mitochondrial", "TSS_fragments", "peak_region_fragments", "on_target_fragments", "enhancer_region_fragments", "blacklist_region_fragments", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_percentile", "nCount_RNA", "nFeature_RNA"),
    pt.size = 0),
    filename = "031_Hp0_NOV_filtered_five.QC.pdf",
    width = 20,
    height = 16,
    limitsize = FALSE
)


# Merging fragments.tsv files for BGI and NOV datasets

# Decompress files and add the same cell prefix as was added to the Seurat object for BGI samples
gzip -dc snATACseq/10_6_M_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"BGI_M1_"$4,$5}' - > BGI_M1_fragments.tsv &
gzip -dc snATACseq/11_1_M_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"BGI_M2_"$4,$5}' - > BGI_M2_fragments.tsv &
gzip -dc snATACseq/11_6_M_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"BGI_M3_"$4,$5}' - > BGI_M3_fragments.tsv &
gzip -dc snATACseq/11_3_F_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"BGI_F1_"$4,$5}' - > BGI_F1_fragments.tsv &

# Merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n BGI_M1_fragments.tsv BGI_M2_fragments.tsv BGI_M3_fragments.tsv BGI_F1_fragments.tsv > BGI.combined.fragments.tsv &

# Block gzip compress the merged file
bgzip -@ 20 BGI.combined.fragments.tsv # -@ 20 uses 20 threads

# Index the bgzipped file
tabix -p bed BGI.combined.fragments.tsv.gz


# Decompress files and add the same cell prefix as was added to the Seurat object for NOV samples
gzip -dc snATACseq/11_2_F_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"NOV_F1_"$4,$5}' - > NOV_F1_fragments.tsv &
gzip -dc snATACseq/17_1_M_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"NOV_M1_"$4,$5}' - > NOV_M1_fragments.tsv &
gzip -dc snATACseq/17_5_M_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"NOV_M2_"$4,$5}' - > NOV_M2_fragments.tsv &
gzip -dc snATACseq/17_1_F_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"NOV_F2_"$4,$5}' - > NOV_F2_fragments.tsv &
gzip -dc snATACseq/17_6_F_snATAC/outs/fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"NOV_F3_"$4,$5}' - > NOV_F3_fragments.tsv &

# Merge files (avoids having to re-sort)
sort -m -k 1,1V -k2,2n NOV_F1_fragments.tsv NOV_M1_fragments.tsv NOV_M2_fragments.tsv NOV_F2_fragments.tsv NOV_F3_fragments.tsv > NOV.combined.fragments.tsv &

# Block gzip compress the merged file
bgzip -@ 20 NOV.combined.fragments.tsv # -@ 20 uses 20 threads

# Index the bgzipped file
tabix -p bed NOV.combined.fragments.tsv.gz


# Merge all BGI and NOV fragment files into a single file
sort -m -k 1,1V -k2,2n BGI_M1_fragments.tsv BGI_M2_fragments.tsv BGI_M3_fragments.tsv BGI_F1_fragments.tsv NOV_F1_fragments.tsv NOV_M1_fragments.tsv NOV_M2_fragments.tsv NOV_F2_fragments.tsv NOV_F3_fragments.tsv > BGI.NOV.combined.fragments.tsv &

# Block gzip compress the merged file
bgzip -@ 64 BGI.NOV.combined.fragments.tsv # -@ 64 uses 64 threads

# Index the bgzipped file
tabix -p bed BGI.NOV.combined.fragments.tsv.gz

# Define a function to add gene activity to the Seurat object
AddGeneActivity <- function(combined, fragment_path, ensembl) {
    fragments <- CreateFragmentObject(
        path = fragment_path,
        cells = colnames(combined),
        validate.fragments = FALSE
    )

    gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
    seqlevelsStyle(gene.coords) <- 'UCSC'
    genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
    genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

    # Build a gene by cell matrix
    gene.activities <- FeatureMatrix(
        fragments = fragments,
        features = genebodyandpromoter.coords,
        cells = colnames(combined),
        verbose = TRUE
    )

    # Convert rownames from chromosomal coordinates into gene names
    gene.key <- genebodyandpromoter.coords$gene_name
    names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
    rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
    gene.activities <- gene.activities[rownames(gene.activities) != "", ]

    # Add the gene activity matrix to the Seurat object as a new assay and normalize it
    combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
    combined <- NormalizeData(
        object = combined,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(combined$nCount_RNA)
    )

    return(combined)
}

# Load Seurat objects and add gene activity
Hp0_BGI_filtered_four <- AddGeneActivity(Hp0_BGI_filtered_four, fragment_path = "BGI.combined.fragments.tsv.gz", ensembl = EnsDb.Hsapiens.v86)
Hp0_NOV_filtered_five <- AddGeneActivity(Hp0_NOV_filtered_five, fragment_path = "NOV.combined.fragments.tsv.gz", ensembl = EnsDb.Hsapiens.v86)

# Set the default assay to 'common_ATAC'
DefaultAssay(Hp0_BGI_filtered_four) <- "common_ATAC"
DefaultAssay(Hp0_NOV_filtered_five) <- "common_ATAC"

# Find overlapping peaks between BGI and NOV datasets
Hp0_BGI.intersecting.regions <- findOverlaps(query = Hp0_BGI_filtered_four, subject = Hp0_NOV_filtered_five)
Hp0_BGI.intersections <- unique(queryHits(Hp0_BGI.intersecting.regions))

Hp0_NOV.intersecting.regions <- findOverlaps(query = Hp0_NOV_filtered_five, subject = Hp0_BGI_filtered_four)
Hp0_NOV.intersections <- unique(queryHits(Hp0_NOV.intersecting.regions))

Hp0_BGI_peaks.use <- sort(granges(Hp0_BGI_filtered_four)[Hp0_BGI.intersections])
Hp0_NOV_peaks.use <- sort(granges(Hp0_NOV_filtered_five)[Hp0_NOV.intersections])

# Count fragments per cell overlapping the set of peaks
NOV_peaks_BGI <- FeatureMatrix(
    fragments = Fragments(Hp0_BGI_filtered_four),
    features = Hp0_NOV_peaks.use,
    cells = colnames(Hp0_BGI_filtered_four)
)

# Create a new assay and add it to the BGI dataset
Hp0_BGI_filtered_four[['novIntPeaks']] <- CreateChromatinAssay(
    counts = NOV_peaks_BGI,
    min.cells = 1,
    ranges = Hp0_NOV_peaks.use,
    genome = 'hg38'
)

DefaultAssay(Hp0_BGI_filtered_four) <- 'novIntPeaks'
Hp0_BGI_filtered_four <- RunTFIDF(Hp0_BGI_filtered_four)

# Prepare NOV dataset for integration
peaknames <- GRangesToString(grange = Hp0_NOV_peaks.use)
Seqinfo.hg38 <- readRDS("Seqinfo.hg38.Rds")
Hp0_NOV_filtered_five[['novIntPeaks']] <- CreateChromatinAssay(
    counts = GetAssayData(Hp0_NOV_filtered_five, assay = "common_ATAC", slot = "counts")[peaknames, ],
    ranges = Hp0_NOV_peaks.use,
    genome = Seqinfo.hg38
)

# Run TF-IDF for the new assay
DefaultAssay(Hp0_NOV_filtered_five) <- "novIntPeaks"
Hp0_NOV_filtered_five <- RunTFIDF(Hp0_NOV_filtered_five)

# Find top features for BGI and NOV datasets
Hp0_BGI_filtered_four <- FindTopFeatures(Hp0_BGI_filtered_four, min.cutoff = 50)
Hp0_NOV_filtered_five <- FindTopFeatures(Hp0_NOV_filtered_five, min.cutoff = 50)

# Merge BGI and NOV datasets
unintegrated <- merge(Hp0_BGI_filtered_four, Hp0_NOV_filtered_five)
DefaultAssay(unintegrated) <- "novIntPeaks"
unintegrated <- RunTFIDF(unintegrated)
unintegrated <- FindTopFeatures(unintegrated, min.cutoff = 50)
unintegrated <- RunSVD(unintegrated)
unintegrated <- RunUMAP(unintegrated, reduction = 'lsi', dims = 2:30)

# Plot UMAP for unintegrated data
p1 <- DimPlot(unintegrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
ggsave(p1, filename = "045_Hp0_BGI_NOV_filtered_unintegrated.novIntPeaks.pdf")

# Run Harmony integration
library(harmony)

hm.integrated <- RunHarmony(
    object = unintegrated,
    group.by.vars = 'dataset',
    reduction = 'lsi',
    assay.use = 'novIntPeaks',
    project.dim = FALSE
)

# Re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
ggsave(
    DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration"),
    filename = "050_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.UMAP.pdf"
)

# Find neighbors and clusters
hm.integrated <- hm.integrated %>%
    FindNeighbors(reduction = 'harmony', dims = 2:30) %>%
    FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

hm.integrated <- hm.integrated %>%
    FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.0)

# Plot UMAP with clusters
ggsave(UMAPPlot(hm.integrated, label = TRUE, repel = TRUE), filename = "051_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.res1.0.UMAP.pdf")
ggsave(UMAPPlot(hm.integrated, label = TRUE, repel = TRUE), filename = "051_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.res1.6.UMAP.pdf")

# Generate violin plots for QC metrics
ggsave(
    VlnPlot(hm.integrated, c(
        "nCount_ATAC", "nFeature_ATAC", "nCount_common_ATAC", "nFeature_common_ATAC", "duplicate", "mitochondrial", "TSS_fragments", "peak_region_fragments", "on_target_fragments", "enhancer_region_fragments", "blacklist_region_fragments", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_percentile", "nCount_RNA", "nFeature_RNA"),
    pt.size = 0),
    filename = "052_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.QC.pdf",
    width = 20,
    height = 16,
    limitsize = FALSE
)

# Plot UMAP split by clusters
p5 <- DimPlot(hm.integrated, split.by = 'seurat_clusters', pt.size = 0.1, label = TRUE)
ggsave(p5, filename = "053_1_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.UMAP.by.cluster.pdf", width = 7 * 28, height = 7, limitsize = FALSE)

# Save the final integrated object
saveRDS(hm.integrated, "053_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.res1.6.UMAP.Rds", compress = FALSE)


# Load necessary libraries
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(motifmatchr)
library(future)

# Load the integrated Seurat object and JASPAR motifs
hm.integrated <- readRDS("053_Hp0_BGI_NOV_filtered_hm_integrated.novIntPeaks.res1.6.UMAP.Rds")
pfm.JASPAR2020 <- readRDS("new_signac_Mp0/JASPAR2020.9606.10090.pfm.Rds")
Seqinfo.hg38 <- readRDS("Seqinfo.hg38.Rds")

# Create the motif matrix
motif.matrix <- CreateMotifMatrix(
  features = Signac::StringToGRanges(regions = rownames(hm.integrated[['novIntPeaks']]@meta.features), sep = c(":", "-")),
  pwm = pfm.JASPAR2020,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  use.counts = FALSE,
  sep = c("-", "-")
)

# Match motifs to the genomic regions
motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm.JASPAR2020,
  subject = Signac::StringToGRanges(regions = rownames(hm.integrated[['novIntPeaks']]@meta.features), sep = c(":", "-")),
  out = 'positions',
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Create the motif object
motif <- CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pfm.JASPAR2020
)

# Add the motif object to the Seurat object
hm.integrated[['novIntPeaks']] <- SetAssayData(
  object = hm.integrated[['novIntPeaks']],
  slot = 'motifs',
  new.data = motif
)

# Calculate region statistics
hm.integrated <- RegionStats(
  object = hm.integrated,
  assay = 'novIntPeaks',
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = TRUE
)

# Run ChromVAR
DefaultAssay(hm.integrated) <- "novIntPeaks"
plan("multiprocess", workers = 16)
hm.integrated.chromVAR <- RunChromVAR(
  object = hm.integrated[['novIntPeaks']],
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Add ChromVAR results to the Seurat object
hm.integrated[['novIntPeaks_JASPAR_chromVAR']] <- hm.integrated.chromVAR

# Remove poorly clustered cells (cluster 21)
hm.integrated.rm21 <- subset(hm.integrated, idents = c(21), invert = TRUE)
DefaultAssay(hm.integrated.rm21) <- 'novIntPeaks'
hm.integrated.rm21 <- RunUMAP(hm.integrated.rm21, dims = 2:30, reduction = 'harmony')
hm.integrated.rm21 <- hm.integrated.rm21 %>%
  FindNeighbors(reduction = 'harmony', dims = 2:30) %>%
  FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Further remove additional clusters (25 and 26)
hm.integrated.rm21.rm2526 <- subset(hm.integrated.rm21, idents = c(25, 26), invert = TRUE)
DefaultAssay(hm.integrated.rm21.rm2526) <- 'novIntPeaks'
hm.integrated.rm21.rm2526 <- RunUMAP(hm.integrated.rm21.rm2526, dims = 2:30, reduction = 'harmony')
hm.integrated.rm21.rm2526 <- hm.integrated.rm21.rm2526 %>%
  FindNeighbors(reduction = 'harmony', dims = 2:30) %>%
  FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Normalize RNA data
DefaultAssay(hm.integrated.rm21.rm2526) <- "RNA"
hm.integrated.rm21.rm2526 <- NormalizeData(hm.integrated.rm21.rm2526)

# Set cell identity levels
Idents(hm.integrated.rm21.rm2526) <- factor(Idents(hm.integrated.rm21.rm2526), levels= c(0, 2, 16, 10, 15, 6, 20, 9, 19, 7, 17, 11, 1, 4, 3, 5, 21, 25))

# Subset N_UE cells
hm.integrated.rm21.rm2526.N_UE <- subset(hm.integrated.rm21.rm2526, idents = c(0, 2, 16, 10, 15, 6, 20, 9, 19, 7, 17, 11, 1, 4, 3, 5, 21, 25))

# Process without Harmony integration
DefaultAssay(hm.integrated.rm21.rm2526.N_UE) <- "novIntPeaks"
hm.integrated.rm21.rm2526.N_UE <- RunTFIDF(hm.integrated.rm21.rm2526.N_UE)
hm.integrated.rm21.rm2526.N_UE <- FindTopFeatures(hm.integrated.rm21.rm2526.N_UE, min.cutoff = 50)
hm.integrated.rm21.rm2526.N_UE <- RunSVD(hm.integrated.rm21.rm2526.N_UE)
hm.integrated.rm21.rm2526.N_UE <- RunUMAP(hm.integrated.rm21.rm2526.N_UE, reduction = 'lsi', dims = 2:30)

# Save UMAP plot
p1 <- DimPlot(hm.integrated.rm21.rm2526.N_UE, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("hm.integrated.rm21.rm2526.N_UE")
ggsave(p1, filename = "074_Hp0_BGI_NOV_filtered_hm.integrated.rm21.rm2526.N_UE.novIntPeaks.pdf")

# Subset and identify cells for each dataset
hm.integrated.rm21.rm2526.N_UE.BGI_F1 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "BGI_F1")
hm.integrated.rm21.rm2526.N_UE.BGI_M1 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "BGI_M1")
hm.integrated.rm21.rm2526.N_UE.BGI_M2 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "BGI_M2")
hm.integrated.rm21.rm2526.N_UE.BGI_M3 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "BGI_M3")
hm.integrated.rm21.rm2526.N_UE.BGI.cells <- c(
  WhichCells(hm.integrated.rm21.rm2526.N_UE.BGI_F1),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.BGI_M1),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.BGI_M2),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.BGI_M3)
)

hm.integrated.rm21.rm2526.N_UE.NOV_F1 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "NOV_F1")
hm.integrated.rm21.rm2526.N_UE.NOV_F2 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "NOV_F2")
hm.integrated.rm21.rm2526.N_UE.NOV_F3 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "NOV_F3")
hm.integrated.rm21.rm2526.N_UE.NOV_M1 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "NOV_M1")
hm.integrated.rm21.rm2526.N_UE.NOV_M2 <- subset(hm.integrated.rm21.rm2526.N_UE, subset = dataset == "NOV_M2")
hm.integrated.rm21.rm2526.N_UE.NOV.cells <- c(
  WhichCells(hm.integrated.rm21.rm2526.N_UE.NOV_F1),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.NOV_F2),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.NOV_F3),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.NOV_M1),
  WhichCells(hm.integrated.rm21.rm2526.N_UE.NOV_M2)
)

# Load BGI and NOV filtered datasets
DefaultAssay(Hp0_BGI_filtered_four) <- "common_ATAC"
DefaultAssay(Hp0_NOV_filtered_five) <- "common_ATAC"

# Subset N_UE cells from BGI and NOV datasets
Hp0_BGI_filtered_four.N_UE <- subset(Hp0_BGI_filtered_four, cells = hm.integrated.rm21.rm2526.N_UE.BGI.cells)
Hp0_NOV_filtered_five.N_UE <- subset(Hp0_NOV_filtered_five, cells = hm.integrated.rm21.rm2526.N_UE.NOV.cells)

# Remove intermediate objects to free up memory
rm(hm.integrated.rm21.rm2526.N_UE.BGI_F1, hm.integrated.rm21.rm2526.N_UE.BGI_M1, hm.integrated.rm21.rm2526.N_UE.BGI_M2, hm.integrated.rm21.rm2526.N_UE.BGI_M3, hm.integrated.rm21.rm2526.N_UE.NOV_F1, hm.integrated.rm21.rm2526.N_UE.NOV_F2, hm.integrated.rm21.rm2526.N_UE.NOV_F3, hm.integrated.rm21.rm2526.N_UE.NOV_M1, hm.integrated.rm21.rm2526.N_UE.NOV_M2)
gc()

# Identify intersecting peaks
Hp0_BGI.N_UE.intersecting.regions <- findOverlaps(query = Hp0_BGI_filtered_four.N_UE, subject = Hp0_NOV_filtered_five.N_UE)
Hp0_BGI.N_UE.intersections <- unique(queryHits(Hp0_BGI.N_UE.intersecting.regions))
Hp0_NOV.N_UE.intersecting.regions <- findOverlaps(query = Hp0_NOV_filtered_five.N_UE, subject = Hp0_BGI_filtered_four.N_UE)
Hp0_NOV.N_UE.intersections <- unique(queryHits(Hp0_NOV.N_UE.intersecting.regions))

# Subset intersecting peaks
Hp0_BGI_peaks.N_UE.use <- sort(granges(Hp0_BGI_filtered_four.N_UE)[Hp0_BGI.N_UE.intersections])
Hp0_NOV_peaks.N_UE.use <- sort(granges(Hp0_NOV_filtered_five.N_UE)[Hp0_NOV.N_UE.intersections])

# Merge BGI and NOV subsets
hm.integrated.rm21.rm2526.N_UE.BGI <- subset(hm.integrated.rm21.rm2526.N_UE, cells = hm.integrated.rm21.rm2526.N_UE.BGI.cells)
hm.integrated.rm21.rm2526.N_UE.NOV <- subset(hm.integrated.rm21.rm2526.N_UE, cells = hm.integrated.rm21.rm2526.N_UE.NOV.cells)
unintegrated <- merge(hm.integrated.rm21.rm2526.N_UE.BGI, hm.integrated.rm21.rm2526.N_UE.NOV)

# Save the unintegrated object
DefaultAssay(unintegrated) <- "novIntPeaks"

# Process the unintegrated object
unintegrated <- RunTFIDF(unintegrated)
unintegrated <- FindTopFeatures(unintegrated, min.cutoff = 50)
unintegrated <- RunSVD(unintegrated)

# Run Harmony integration
hm.integrated <- RunHarmony(
  object = unintegrated,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'novIntPeaks',
  project.dim = FALSE
)

# Re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')

# Save UMAP plot
p5 <- DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
ggsave(p5, filename = "082_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.novIntPeaks.UMAP.pdf")

# Find neighbors and clusters
hm.integrated <- FindNeighbors(hm.integrated, reduction = 'harmony', dims = 2:30)
hm.integrated <- hm.integrated %>%
  FindClusters(verbose = TRUE, algorithm = 3, resolution = 1.6)

# Save UMAP plot with clusters
ggsave(UMAPPlot(hm.integrated, label = TRUE, repel = TRUE), filename = "083_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.novIntPeaks.res1.6.UMAP.pdf")

#Remove cluster 8 (Vasculature)
hm.integrated.rm8 <- subset(hm.integrated, idents=c(8), invert=T)