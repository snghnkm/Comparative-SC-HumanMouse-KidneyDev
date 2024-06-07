# 010_Cicero_analysis.r

# Load necessary libraries
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(cicero)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set working directory
setwd("newest_all_inclusive_analysis")

# Load Seurat object
coembed.subset.both.human.atac <- readRDS("coembed.subset.both.human.atac.MACS2_reduced_500.motif.Rds")

# Set the default assay
DefaultAssay(coembed.subset.both.human.atac) <- "MACS2_reduced_500"

# Define a function to run Cicero connections
cicero_conns <- function(atac){
    # Convert to CellDataSet format and make the Cicero object
    atac.cds <- as.cell_data_set(x = atac)
    atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)

    # Get the chromosome sizes from the Seurat object
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)

    # Convert chromosome sizes to a dataframe
    genome.df <- data.frame("chr" = names(genome), "length" = genome)

    # Run Cicero to find connections
    conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 100)

    return(conns)
}

# Convert to CellDataSet format and make the Cicero object
atac.cds <- SeuratWrappers::as.cell_data_set(x = coembed.subset.both.human.atac)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP)

# Get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)

# Convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# Run Cicero to find connections
conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 100)

# Generate Cicero cis-co-accessibility networks (CCANs)
ccans <- generate_ccans(conns)

# Save the results
saveRDS(conns, "coembed.subset.both.human.atac.MACS2_reduced_500.motif.cicero.conns.Rds")
saveRDS(ccans, "coembed.subset.both.human.atac.MACS2_reduced_500.motif.cicero.ccans.Rds")
