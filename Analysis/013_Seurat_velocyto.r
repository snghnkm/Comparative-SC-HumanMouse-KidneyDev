# 013_Seurat_velocyto.r

# Step 1. Seurat bam file to loom file format conversion
# Velocyto_combined_loom.py

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import loompy
import matplotlib
matplotlib.use('pdf')

# Reading loom files from various directories
loom_files = [
    "RNAVelocity/old_Mp0_F1/possorted_genome_bam_X5OC1.loom",
    "RNAVelocity/old_Mp0_F2/possorted_genome_bam_K09W9.loom",
    "RNAVelocity/old_Mp0_F3/possorted_genome_bam_6TO5C.loom",
    "RNAVelocity/old_Mp0_M4/possorted_genome_bam_RZFWW.loom",
    "RNAVelocity/old_Mp0_M5/possorted_genome_bam_GBJ2R.loom",
    "RNAVelocity/old_Mp0_M6/possorted_genome_bam_DB4HH.loom",
    "RNAVelocity/old_Mp0_M7/possorted_genome_bam_EN01H.loom",
    "RNAVelocity/new_Mp0_F4/possorted_genome_bam_JVDD4.loom",
    "RNAVelocity/new_Mp0_F5/possorted_genome_bam_0BK5O.loom",
    "RNAVelocity/new_Mp0_F6/possorted_genome_bam_1CNYP.loom",
    "RNAVelocity/new_Mp0_M8/possorted_genome_bam_6OOFP.loom",
    "RNAVelocity/new_Mp0_M9/possorted_genome_bam_XT6Y9.loom",
    "RNAVelocity/new_Mp0_M10/possorted_genome_bam_I5T6Y.loom"
]

# Combine loom files into a single file
loompy.combine(loom_files, "RNAVelocity/Mp0_combined.loom", key="Accession")

# Step 2. In R, extract cell barcode, UMAP, and cluster information

# (R code snippet, to be executed in an R environment)
# Extract cell barcode, UMAP, and cluster information from Seurat object
seurat_object.cells <- Cells(seurat_object)
seurat_object.cells <- seurat_object.cells %>% str_replace("F1_", "possorted_genome_bam_X5OC1:")
seurat_object.cells <- seurat_object.cells %>% str_replace("F2_", "possorted_genome_bam_K09W9:")
seurat_object.cells <- seurat_object.cells %>% str_replace("F3_", "possorted_genome_bam_6TO5C:")
seurat_object.cells <- seurat_object.cells %>% str_replace("F4_", "possorted_genome_bam_JVDD4:")
seurat_object.cells <- seurat_object.cells %>% str_replace("F5_", "possorted_genome_bam_0BK5O:")
seurat_object.cells <- seurat_object.cells %>% str_replace("F6_", "possorted_genome_bam_1CNYP:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M4_", "possorted_genome_bam_RZFWW:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M5_", "possorted_genome_bam_GBJ2R:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M6_", "possorted_genome_bam_DB4HH:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M7_", "possorted_genome_bam_EN01H:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M8_", "possorted_genome_bam_6OOFP:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M9_", "possorted_genome_bam_XT6Y9:")
seurat_object.cells <- seurat_object.cells %>% str_replace("M10_", "possorted_genome_bam_I5T6Y:")
seurat_object.cells <- seurat_object.cells %>% paste("x", sep="")
write.csv(seurat_object.cells, file = "seurat_object.cellID_obs.2.csv", row.names = FALSE)
write.csv(Embeddings(seurat_object, reduction = "umap"), file = "seurat_object.cell_embeddings.2.csv")
write.csv(seurat_object@meta.data$seurat_clusters, file = "seurat_object.clusters.2.csv")

# Step 3. In Python, load loom file and cluster information, then run scvelo

# Load combined loom file
loom = anndata.read_loom("RNAVelocity/Mp0_combined.loom")

# Read cell information and embeddings from CSV files
sample_obs = pd.read_csv("final_trajectory_analysis/seurat_object.cellID_obs.2.csv")
umap_cord = pd.read_csv("final_trajectory_analysis/seurat_object.cell_embeddings.2.csv")
cell_clusters = pd.read_csv("final_trajectory_analysis/seurat_object.clusters.2.csv")

# Filter loom object to include only relevant cells
loom_filtered = loom[np.isin(loom.obs.index.to_series(), sample_obs.x)]

# Rename columns and merge data
umap_cord = umap_cord.rename(columns={'Unnamed: 0': 'Cell ID'})
umap_cord['Cell ID'] = sample_obs["x"]
umap_cord['cluster'] = cell_clusters["x"]

loom_filtered_index = pd.DataFrame(loom_filtered.obs.index).rename(columns={'CellID': 'Cell ID'})
umap_ordered = loom_filtered_index.merge(umap_cord, on="Cell ID").set_index("Cell ID")

# Ensure indices match
assert loom_filtered.obs_names.equals(umap_ordered.index)

# Add UMAP and cluster information to the loom object
loom_filtered.obsm['X_umap'] = np.array(umap_ordered.iloc[:, 0:2])
loom_filtered.obsm['seurat_clusters'] = np.array(umap_ordered.iloc[:, 2])

# Preprocess and analyze RNA velocity
scv.pp.filter_and_normalize(loom_filtered)
scv.pp.moments(loom_filtered)
scv.tl.velocity(loom_filtered, mode="stochastic")
scv.tl.velocity_graph(loom_filtered)

# Plot velocity embedding grid and save as PDF
scv.pl.velocity_embedding_grid(loom_filtered, basis='umap', dpi=800, size=6, color_map='Set3', color=loom_filtered.obsm['seurat_clusters'], arrow_length=4, arrow_color="grey")
plt.pyplot.savefig("final_trajectory_analysis/seurat_object.velocity_arrows.pdf")

scv.pl.velocity_embedding_grid(loom_filtered, basis='umap', dpi=1200, size=20, color_map='Set3', color=loom_filtered.obsm['seurat_clusters'], arrow_length=5, arrow_color="grey")
plt.pyplot.savefig("final_trajectory_analysis/seurat_object.velocity_arrows.2.pdf")
