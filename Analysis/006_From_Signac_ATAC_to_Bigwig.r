# 006_From_Signac_ATAC_to_Bigwig.R
# An example script for an integrated snATACseq data with multiple batches

# Step 1 in R environment: Extract barcodes 

# Load necessary libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)

# Read in the dataset
coembed.subset.both.human.atac <- readRDS("all_inclusive_analysis/coembed.subset.both.human.atac.Rds")
table(coembed.subset.both.human.atac@meta.data$dataset)

# Initialize a dataframe to store barcodes
barcodes <- data.frame()

# Loop through each cluster and extract barcodes
for (each_cluster in unique(coembed.subset.both.human.atac@active.ident)) {
  print(each_cluster)
  cells_each_cluster <- WhichCells(object = coembed.subset.both.human.atac, idents = each_cluster)
  
  for (each_pattern in c("BGI_F1", "BGI_M1", "BGI_M2", "BGI_M3", "NOV_F1", "NOV_F2", "NOV_F3", "NOV_M1", "NOV_M2")) {
    print(each_pattern)
    barcodes_each_cluster_each_pattern <- grep(pattern = each_pattern, x = cells_each_cluster, perl = TRUE, value = TRUE)
    clusters <- rep(each_cluster, length(barcodes_each_cluster_each_pattern))
    patterns <- rep(each_pattern, length(barcodes_each_cluster_each_pattern))
    barcodes <- rbind(barcodes, data.frame(barcodes = barcodes_each_cluster_each_pattern, cluster = clusters, pattern = patterns))
  }
}

# Trim barcodes
barcodes.trimmed <- barcodes %>% tidyr::separate(barcodes, c("company", "sample", "barcode"), "_")

# Write barcodes to CSV files
for (each_pattern in c("BGI_F1", "BGI_M1", "BGI_M2", "BGI_M3", "NOV_F1", "NOV_F2", "NOV_F3", "NOV_M1", "NOV_M2")) {
  write.table(barcodes.trimmed %>% filter(pattern == each_pattern) %>% dplyr::select(barcode, cluster),
              file = paste0(each_pattern, ".coembed.subset.both.human.atac.barcodes.csv"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# Step 2 in bash script: Filter BAM files by using Cell barcode for each cluster/sample

# Bash commands to filter BAM files using sinto
# Activate sinto environment
conda activate sinto

# Run sinto filterbarcodes for each BAM file and barcode set
nohup sinto filterbarcodes -b snATACseq/11_3_F_snATAC/outs/possorted_bam.bam -c ../BGI_F1.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/10_6_M_snATAC/outs/possorted_bam.bam -c ../BGI_M1.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/11_1_M_snATAC/outs/possorted_bam.bam -c ../BGI_M2.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/11_6_M_snATAC/outs/possorted_bam.bam -c ../BGI_M3.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/11_2_F_snATAC/outs/possorted_bam.bam -c ../NOV_F1.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/17_1_F_snATAC/outs/possorted_bam.bam -c ../NOV_F2.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/17_6_F_snATAC/outs/possorted_bam.bam -c ../NOV_F3.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/17_1_M_snATAC/outs/possorted_bam.bam -c ../NOV_M1.coembed.subset.both.human.atac.barcodes -p 64 &
nohup sinto filterbarcodes -b snATACseq/17_5_M_snATAC/outs/possorted_bam.bam -c ../NOV_M2.coembed.subset.both.human.atac.barcodes -p 64 &

# Step 3 python generated bash scripts: Process BAM files and generate bigwig files

# Python: Merge and index BAM files, then generate BigWig files

# Merge BAM files for each cluster/sample
for each in range(0, 24):
  print(f"samtools merge {each}.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.merged.bam \
  ./BGI_F1.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./BGI_M1.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./BGI_M2.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./BGI_M3.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./NOV_F1.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./NOV_F2.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./NOV_F3.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./NOV_M1.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam \
  ./NOV_M2.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8/{each}.bam")

# Index the merged BAM files
for each in range(0, 24):
  print(f"samtools index {each}.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.merged.bam")

# Generate BigWig files using deeptools bamCoverage
for each in range(0, 24):
  print(f"bamCoverage -b {each}.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.merged.bam -o {each}.Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.bw -p max --binSize 100 --ignoreForNormalization chrX --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --extendReads")
