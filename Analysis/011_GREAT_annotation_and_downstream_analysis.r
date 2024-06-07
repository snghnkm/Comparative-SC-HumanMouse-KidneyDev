# 011_GREAT_annotation_and_downstream_analysis.r

# Load required libraries
library(tidyverse)
library(Signac)
library(Seurat)
library(GenomicRanges)

# Load GREAT annotation files
common_bed_GREAT <- read_tsv("newest_all_inclusive_analysis/coembed.subset.both.human.atac.DAR.AllMarkers.common.bed.GREAT.gene_association.txt", col_names = FALSE)
human_only_bed_GREAT <- read_tsv("newest_all_inclusive_analysis/coembed.subset.both.human.atac.DAR.AllMarkers.human_only.bed.GREAT.gene_association.txt", col_names = FALSE)

# Load RNA DEGs
rna_DEG_all_markers <- readRDS("coembed.subset.both.human.rna.DEG.AllMarkers.Rds")

# Load gene-to-region mapping
gene2region <- read_tsv("newest_all_inclusive_analysis/coembed.subset.both.human.atac.DAR.AllMarkers.merged.p_val.0.05.named.bed.gene2region.txt", col_names = FALSE)

# Filter DEGs
deg_genes <- rna_DEG_all_markers %>% filter(p_val_adj < 0.05) %>% pull(gene) %>% unique()
deg_with_GREATanno <- intersect(gene2region %>% pull(X1), deg_genes)

# Add count of regions
gene2region_count <- gene2region %>% mutate(count = str_count(X2, pattern = fixed(",")) + 1)

# Summarize counts
deg_gene_counts <- gene2region_count %>% filter(X1 %in% deg_genes) %>% summarise(total_count = sum(count))

# Extract peaks for DEG with GREAT annotation
deg_peaks <- gene2region_count %>% 
  filter(X1 %in% deg_genes) %>% 
  pull(X2) %>% 
  str_extract_all('chr(\\d*|[X-Y]):\\d*-\\d*') %>% 
  unlist() %>% 
  unique() %>% 
  StringToGRanges(sep = c(":", "-"))

saveRDS(deg_peaks, "human_DEG_with_GREATanno_peaks.Rds")

# Extract peaks for conserved and human-specific DEGs
conserved_peaks <- gene2region_count %>% 
  filter(X1 %in% (conservedmarkers.pct_diff.filtered %>% pull(rowname) %>% unique())) %>% 
  pull(X2) %>% 
  str_extract_all('chr(\\d*|[X-Y]):\\d*-\\d*') %>% 
  unlist() %>% 
  unique() %>% 
  StringToGRanges(sep = c(":", "-"))

saveRDS(conserved_peaks, "human_conserved_DEG_with_GREATanno_peaks.Rds")

human_only_peaks <- gene2region_count %>% 
  filter(X1 %in% human_specific_DEGs) %>% 
  pull(X2) %>% 
  str_extract_all('chr(\\d*|[X-Y]):\\d*-\\d*') %>% 
  unlist() %>% 
  unique() %>% 
  StringToGRanges(sep = c(":", "-"))

saveRDS(human_only_peaks, "human_only_DEG_with_GREATanno_peaks.Rds")

# Offset peaks by 1
offset_peaks <- function(peaks) {
  peaks %>% 
    data.frame() %>% 
    mutate(start = start + 1) %>% 
    makeGRangesFromDataFrame() %>% 
    GRangesToString()
}

saveRDS(offset_peaks(deg_peaks), "human_DEG_with_GREATanno_peaks.offset1.Rds")
saveRDS(offset_peaks(conserved_peaks), "human_conserved_DEG_with_GREATanno_peaks.offset1.Rds")
saveRDS(offset_peaks(human_only_peaks), "human_only_DEG_with_GREATanno_peaks.offset1.Rds")

# Load ATAC and RNA datasets
atac_data <- readRDS("coembed.subset.both.human.atac.MACS2_reduced_500.Rds")
rna_data <- readRDS("coembed.subset.both.human.rna.Rds")

# Run TF-IDF on ATAC data and normalize RNA data
atac_data <- RunTFIDF(atac_data)
rna_data <- NormalizeData(rna_data)

# Remove clusters 12 and 17
atac_data <- subset(atac_data, idents = c(12, 17), invert = TRUE)
rna_data <- subset(rna_data, idents = c(12, 17), invert = TRUE)

saveRDS(atac_data, "coembed.subset.both.human.atac.rm1217.Rds")
saveRDS(rna_data, "coembed.subset.both.human.rna.rm1217.Rds")

# Compute average expression
rna_avg <- AverageExpression(rna_data, assays = "RNA", return.seurat = FALSE, slot = "data")
atac_avg <- AverageExpression(atac_data, assays = "MACS2_reduced_500", return.seurat = FALSE, slot = "data")

saveRDS(rna_avg, "coembed.subset.both.human.rna.rm1217.avg.RNA.Rds")
saveRDS(atac_avg, "coembed.subset.both.human.atac.rm1217.avg.MACS2_reduced_500.Rds")

# Save VEGFA enhancer data
vegfa_data <- rbind(
  rna_avg$RNA["VEGFA", ],
  atac_avg$MACS2_reduced_500["chr6-43744329-43745600", ]
) %>% t() %>% data.frame()

write_excel_csv(vegfa_data, "VEGFA.enhancer1.csv")

# Create and save plots for enhancers and promoter
plot_and_save <- function(region, filename) {
  p <- atac_avg$MACS2_reduced_500[region, ] %>% 
    data.frame() %>% 
    select(c(`X0`, `X10`, `X2`, `X16`, `X13`, `X19`, `X1`, `X11`, `X7`, `X8`, `X18`, `X9`, `X14`, `X20`, `X6`, `X22`, `X4`, `X5`, `X3`, `X21`, `X23`)) %>% 
    t() %>% 
    data.frame() %>% 
    rownames_to_column() %>% 
    ggbarplot(x = "rowname", y = region)
  
  p2 <- p + coord_flip()
  ggsave(p2, filename = filename, width = 2)
}

plot_and_save("chr6-43744329-43745600", "coembed.subset.both.human.atac.rm1217.avg.chr6-43744329-43745600.pdf")
plot_and_save("chr6-43746183-43748431", "coembed.subset.both.human.atac.rm1217.avg.chr6-43746183-43748431.pdf")
plot_and_save("chr6-43769373-43772158", "coembed.subset.both.human.atac.rm1217.avg.chr6-43769373-43772158.pdf")

# Check dimensions of average expression data
dim(rna_avg$RNA)  # Expected: 4178
dim(atac_avg$MACS2_reduced_500)  # Check number of peaks

# RNA data imputation and peak analysis

# Load required libraries
library(tidyverse)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set working directory
setwd("newest_all_inclusive_analysis")

# Load necessary data
transcriptome_integrated_human <- readRDS("newest_seurat_Mp0_Hp0_integrated/010_Hp0_Mp0_integrated_trial5.human.RNA_norm.Rds")
hp0_N_UE_ATAC <- readRDS("new_signac_Hp0/092_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.Rds")
transfer_anchors <- readRDS("001_transfer.anchors.all_inclusive_analysis.human.Rds")

# Perform data imputation
refdata_RNA <- GetAssayData(transcriptome_integrated_human, assay = "RNA", slot = "data")
imputation <- TransferData(anchorset = transfer_anchors, refdata = refdata_RNA, weight.reduction = hp0_N_UE_ATAC[["harmony"]])

# Load ATAC data and add imputed RNA data
coembed_subset_both_human_atac <- readRDS("coembed.subset.both.human.atac.MACS2_reduced_500.Rds")
coembed_subset_both_human_atac[["imputation_RNA"]] <- imputation
saveRDS(coembed_subset_both_human_atac, "coembed.subset.both.human.atac.MACS2_reduced_500.imputation_RNA.2.Rds", compress = FALSE)

# Normalize ATAC data using TF-IDF
DefaultAssay(coembed_subset_both_human_atac) <- "MACS2_reduced_500"
coembed_subset_both_human_atac <- RunTFIDF(coembed_subset_both_human_atac)
saveRDS(coembed_subset_both_human_atac, "coembed.subset.both.human.atac.MACS2_reduced_500.imputation_RNA.TFIDF.Rds", compress = FALSE)

# Remove clusters 12 and 17
coembed_subset_both_human_atac_rm1217 <- subset(coembed_subset_both_human_atac, idents = c(12, 17), invert = TRUE)

# Extract and save normalized peak data and RNA expression data
norm_peak_data <- GetAssayData(coembed_subset_both_human_atac_rm1217, assay = "MACS2_reduced_500", slot = "data")
norm_rna_exp <- GetAssayData(coembed_subset_both_human_atac_rm1217, assay = "imputation_RNA", slot = "data")
peak_ranges <- granges(coembed_subset_both_human_atac_rm1217)

saveRDS(norm_peak_data, "norm_peak_data.2.Rds", compress = FALSE)
saveRDS(norm_rna_exp, "norm_rna_exp.2.Rds", compress = FALSE)

# Extract and save peak counts and ranges
DefaultAssay(coembed_subset_both_human_atac_rm1217) <- "MACS2_reduced_500"
peak_counts <- GetAssayData(coembed_subset_both_human_atac_rm1217, slot = "counts")
peak_counts <- peak_counts[rowSums(x = peak_counts) > 0, ]
saveRDS(peak_counts, "peak_counts.2.Rds", compress = FALSE)

peak_ranges <- StringToGRanges(regions = rownames(peak_counts), sep = c(":", "-"))
peak_ranges <- peak_ranges[rowSums(x = peak_counts) > 0]
saveRDS(peak_ranges, "peak_ranges.2.Rds", compress = FALSE)
peak_ranges_str <- GRangesToString(peak_ranges, sep = c("-", "-"))
saveRDS(peak_ranges_str, "peak_ranges_str.2.Rds", compress = FALSE)

# Generate background peaks with GC bias
bgpeaks_all <- SummarizedExperiment(
  assays = list(counts = peak_counts), 
  rowRanges = peak_ranges
) %>% 
  addGCBias(genome = BSgenome.Hsapiens.UCSC.hg38) %>%
  getBackgroundPeaks(niterations = 100) %>% 
  as.data.frame()
rownames(bgpeaks_all) <- peak_ranges_str

saveRDS(bgpeaks_all, "bgpeaks_all.2.Rds", compress = FALSE)

# Load GREAT annotation and RNA DEG data
gene2region <- read_tsv("newest_all_inclusive_analysis/coembed.subset.both.human.atac.DAR.AllMarkers.merged.p_val.0.05.named.bed.gene2region.txt", col_names = FALSE)
rna_DEG_all_markers <- readRDS("coembed.subset.both.human.rna.DEG.AllMarkers.Rds")

# Identify human DEGs with GREAT annotation
human_DEG_with_GREATanno <- intersect(gene2region %>% pull(`X1`), rna_DEG_all_markers %>% filter(p_val_adj < 0.05) %>% pull(gene) %>% unique())
saveRDS(human_DEG_with_GREATanno, "human_DEG_with_GREATanno.Rds")

# Load offset peaks for human DEGs with GREAT annotation
human_DEG_with_GREATanno_peaks_offset1 <- readRDS("human_DEG_with_GREATanno_peaks.offset1.Rds")

# Select background peaks
bgpeaks_all_selected <- bgpeaks_all[human_DEG_with_GREATanno_peaks_offset1, ]
saveRDS(bgpeaks_all_selected, "bgpeaks_all.selected.2.Rds", compress = FALSE)

# Save peak counts and ranges for background peak analysis
peak_ranges_str <- GRangesToString(peak_ranges, sep = c("-", "-"))
saveRDS(peak_ranges_str, "peak_ranges_str.Rds", compress = FALSE)

bgpeaks_all <- SummarizedExperiment(
  assays = list(counts = peak_counts), 
  rowRanges = peak_ranges
) %>% 
  addGCBias(genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) %>%
  getBackgroundPeaks(niterations = 100) %>% 
  as.data.frame()
rownames(bgpeaks_all) <- peak_ranges_str

saveRDS(bgpeaks_all, "bgpeaks_all.Rds", compress = FALSE)
bgpeaks_all <- readRDS("bgpeaks_all.Rds")

# Select and save background peaks for human DEGs with GREAT annotation
bgpeaks_all_selected <- bgpeaks_all[human_DEG_with_GREATanno_peaks_offset1, ]
saveRDS(bgpeaks_all_selected, "bgpeaks_all.selected.Rds", compress = FALSE)

# Core DORC analysis

library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)
library(dplyr)
library(Signac)
library(GenomicRanges)
library(stringr)
library(readr)

# Set up parallel processing
cl <- makeCluster(64)
registerDoSNOW(cl)

# Set working directory
setwd("newest_all_inclusive_analysis/DORC")

# Load necessary data
human_DEG_with_GREATanno <- readRDS("human_DEG_with_GREATanno.Rds")
already_done <- list.files(path="./", pattern=".Rds", all.files=TRUE, full.names=F) %>% 
  str_split(pattern=fixed(".")) %>% 
  sapply(head, 1)
human_DEG_with_GREATanno <- setdiff(human_DEG_with_GREATanno, already_done)

# Print remaining genes to be processed
message(human_DEG_with_GREATanno %>% length())

# Perform DORC analysis in parallel
DORC_parallel_res <- foreach(i = seq(human_DEG_with_GREATanno %>% length()), .combine = dplyr::bind_rows, .packages = "dplyr") %dopar% {
  gene_name <- human_DEG_with_GREATanno[i]
  
  message(gene_name)
  gene_exp <- coembed.subset.both.human.rna.rm1217.avg$RNA[gene_name, ] %>% as.numeric()
  chromatin_acc <- coembed.subset.both.human.atac.rm1217.avg$MACS2_reduced_500
  
  GREATannotated_peaks.offset1 <- coembed.subset.both.human.atac.DAR.AllMarkers.merged.p_val.0.05.named.bed.gene2region %>%
    filter(`X1` == gene_name) %>%
    pull(`X2`) %>%
    stringr::str_extract_all('chr(\\d*|[X-Y]):\\d*-\\d*') %>%
    unlist() %>%
    unique() %>%
    Signac::StringToGRanges(sep = c(":", "-")) %>%
    data.frame() %>%
    mutate(start = start + 1) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    Signac::GRangesToString()
  
  norm_obs_list <- tibble()
  
  for (each_peak in GREATannotated_peaks.offset1) {
    obs <- cor(x = gene_exp, y = chromatin_acc[each_peak, ] %>% as.numeric(), method = "pearson")
    bgpeaks_indices <- bgpeaks_all.selected[each_peak, ] %>% as.numeric()
    
    pop_list <- c()
    for (each_bgpeak_index in bgpeaks_indices[1:100]) {
      each_bgpeak <- peak_ranges_str[each_bgpeak_index]
      norm_bgpeak_counts <- as.numeric(chromatin_acc[each_bgpeak, ])
      pop <- cor(x = gene_exp, y = norm_bgpeak_counts, method = "pearson")
      pop_list <- c(pop_list, pop)
    }
    
    pop.mean <- mean(pop_list)
    pop.sd <- sd(pop_list)
    norm_obs <- (obs - pop.mean) / pop.sd
    norm_obs_p <- round(pnorm(abs(norm_obs), lower.tail = FALSE), 100)
    
    norm_obs_list <- bind_rows(norm_obs_list, tibble(
      gene = gene_name,
      peak = each_peak,
      obs_pearson = obs,
      pop_mean_pearson = pop.mean,
      pop_sd_pearson = pop.sd,
      norm_obs_pearson = norm_obs,
      norm_obs_p_pearson = norm_obs_p
    ))
  }
  
  norm_obs_list_sorted <- norm_obs_list %>% arrange(norm_obs_p_pearson)
  saveRDS(norm_obs_list_sorted, paste0("newest_all_inclusive_analysis/DORC/", gene_name, ".norm_obs_list_sorted.Rds"), compress = FALSE)
  return(norm_obs_list_sorted)
}

# Merge results
setwd("newest_all_inclusive_analysis/DORC/")
DORC_res_merged <- tibble()
for (each_file in list.files(path = ".", pattern = ".Rds", all.files = TRUE, full.names = FALSE)) {
  message(each_file)
  DORC_res_merged <- rbind(DORC_res_merged, readRDS(each_file))
}

# Save and analyze merged results
saveRDS(DORC_res_merged, "newest_all_inclusive_analysis/DORC/DORC_res_merged.Rds", compress = FALSE)
write_excel_csv(DORC_res_merged, "newest_all_inclusive_analysis/DORC/DORC_res_merged.csv")

DORC_res_merged <- readRDS("newest_all_inclusive_analysis/DORC/DORC_res_merged.Rds")

# Summary statistics
DORC_res_merged %>% pull(peak) %>% unique() %>% length() # 32782
DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% select(peak, gene) %>% unique() %>% dim() # 11156, 11156/32782
DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% pull(gene) %>% unique() %>% length() # 3426
DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% filter(norm_obs_pearson < 0) %>% select(peak, gene) %>% unique() %>% dim() # 307/11156
DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% filter(norm_obs_pearson > 0) %>% select(peak, gene) %>% unique() %>% dim() # 10891/11156

# Save peaks with both positive and negative correlation
DORC_res_merged %>%
  filter(peak %in% intersect(
    DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% filter(norm_obs_pearson < 0) %>% pull(peak) %>% unique(),
    DORC_res_merged %>% filter(norm_obs_p_pearson < 0.05) %>% filter(norm_obs_pearson > 0) %>% pull(peak) %>% unique()
  )) %>%
  arrange(peak) %>%
  write_excel_csv("DORC_res_merged.selected.both.pos.neg.peaks.csv")

# Caret Random Forest and Volcano Plot Analysis

library(Seurat)
library(Signac)
library(caret)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(tidyr)

# Set working directory
setwd("newest_all_inclusive_analysis")

# Load datasets
coembed.subset.both.human.atac <- readRDS("coembed.subset.both.human.atac.MACS2_reduced_500.Rds")
coembed.subset.both.human.rna <- readRDS("coembed.subset.both.human.rna.Rds")

# Run TF-IDF and Normalize Data
coembed.subset.both.human.atac <- RunTFIDF(coembed.subset.both.human.atac)
coembed.subset.both.human.rna <- NormalizeData(coembed.subset.both.human.rna)

# Remove clusters 12 and 17 in RNA
coembed.subset.both.human.rna.rm1217 <- subset(coembed.subset.both.human.rna, idents = c(12,17), invert=T)

# Compute Average Expression
coembed.subset.both.human.rna.rm1217.avg <- AverageExpression(coembed.subset.both.human.rna.rm1217, assays="RNA", return.seurat = FALSE, slot = "data")
coembed.subset.both.human.atac.avg <- AverageExpression(coembed.subset.both.human.atac, assays="MACS2_reduced_500", return.seurat = FALSE, slot = "data")

# Load DORC annotated data
DORC_res_merged.annotated <- readRDS("DORC_res_merged.annotated.coembed.subset.both.human.rna.DEG.AllMarkers.DEGs.all.Rds")

# Perform Caret Random Forest
# Assuming DORC_res_merged.annotated contains features and labels for training
# Prepare data for Caret
set.seed(123)
train_data <- DORC_res_merged.annotated %>% select(-gene_peak, -label)
train_labels <- DORC_res_merged.annotated$label

train_control <- trainControl(method="cv", number=5)
rf_model <- train(train_data, train_labels, method="rf", trControl=train_control)

# Save the Random Forest model
saveRDS(rf_model, "caret_random_forest_model.Rds")



