# 008_Comparative_open_chromatin_accessibility_analysis.R

# Load required libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(rtracklayer)

# MACS2 Peak Calling - Human
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
qValue_cutoff <- 2
gr_narrowPeaks <- GRangesList()

# Import human narrowPeak files
for (i in setdiff(seq(0, 24), c(12, 17))) {
  file_narrowPeak <- paste0("newest_all_inclusive_analysis/MACS2/human/", i, ".coembed.subset.both.human.atac.bam_peaks.narrowPeak")
  gr_narrowPeak <- rtracklayer::import(file_narrowPeak, format = "BED", extraCols = extraCols_narrowPeak)
  gr_narrowPeak <- gr_narrowPeak[gr_narrowPeak$qValue > qValue_cutoff]
  gr_narrowPeaks <- c(gr_narrowPeaks, gr_narrowPeak)
}

# Reduce peaks
gr_narrowPeaks.reduced <- reduce(x = unlist(gr_narrowPeaks))
export.bed(gr_narrowPeaks.reduced, con='gr_narrowPeaks.human.reduced.bed')

# Reduce peaks with a minimum gap width of 500
gr_narrowPeaks.reduced.500 <- reduce(x = unlist(gr_narrowPeaks), min.gapwidth=500L)
export.bed(gr_narrowPeaks.reduced.500, con='gr_narrowPeaks.human.reduced.500.bed')

# Set the human fragment path
human_fragment_path <- "new_signac_Hp0/BGI.NOV.combined.fragments.tsv.gz"

# Load the Seurat object and rename cells
hp0_N_UE.ATAC <- readRDS("new_signac_Hp0/092_Hp0_BGI_NOV_filtered_N_UE_hm_integrated.rm8.novIntPeaks.Rds")
hp0_N_UE.ATAC.cells <- colnames(hp0_N_UE.ATAC)

coembed.subset.both.human.atac <- readRDS("coembed.subset.both.human.atac.Rds")
coembed.subset.both.human.atac <- RenameCells(coembed.subset.both.human.atac, new.names=hp0_N_UE.ATAC.cells)
coembed.subset.both.human.atac.cells <- colnames(coembed.subset.both.human.atac)

# Create Fragment Object and Feature Matrix
human_frags <- CreateFragmentObject(
  path = human_fragment_path,
  cells = coembed.subset.both.human.atac.cells,
  validate.fragments = TRUE
)

human_counts <- FeatureMatrix(
  fragments = human_frags,
  features = granges(gr_narrowPeaks.reduced.500),
  cells = coembed.subset.both.human.atac.cells
)

# Create Chromatin Assay and save the Seurat object
human_assay <- CreateChromatinAssay(human_counts, fragments = human_fragment_path)
coembed.subset.both.human.atac[['MACS2_reduced_500']] <- human_assay
DefaultAssay(coembed.subset.both.human.atac) <- "MACS2_reduced_500"

# Perform Differential Accessibility Analysis
coembed.subset.both.human.atac.DAR.AllMarkers <- FindAllMarkers(coembed.subset.both.human.atac, assay="MACS2_reduced_500", test.use = 'LR', latent.vars = 'peak_region_fragments', logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)

# Save unique peaks and overlaps
coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks <- coembed.subset.both.human.atac.DAR.AllMarkers %>% pull(gene) %>% unique()
coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks.granges <- coembed.subset.both.human.atac.DAR.AllMarkers %>% pull(gene) %>% unique() %>% StringToGRanges(sep=c('-','-'))

# Import ENCODE annotation and save
extraCols_encode <- c(fullname = "character", shortname = "character", histone = "numeric", abbreviation = "character", accession="character", description="character")
gr_encode <- rtracklayer::import("encodeCcreCombined/hg38/encodeCcreCombined.bb.bed", format = "BED", extraCols = extraCols_encode)

# Find specific peaks
human_only_DAR_with_encode_annotation <- tibble()

coembed.subset.both.human.atac.DAR.AllMarkers.human_only.p.0.05.peaks <- coembed.subset.both.human.atac.DAR.AllMarkers.human_only %>% filter(`p_val_adj` < 0.05) %>% pull(gene) %>% unique()

# Parallel processing setup
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)

closeAllConnections();gc()
cl <- parallel::makeCluster(20)
doSNOW::registerDoSNOW(cl)

human_DAR_with_encode_annotation <- foreach (i = seq(length(coembed.subset.both.human.atac.DAR.AllMarkers.p.0.05.peaks)), .combine=dplyr::bind_rows, .packages = c("dplyr","rtracklayer")) %dopar% {
  peakofinterest = coembed.subset.both.human.atac.DAR.AllMarkers.p.0.05.peaks[i]
  coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks.granges.encode.overlaps <- findOverlaps(coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks.granges, gr_encode) %>% as.data.frame()
  peakofinterest.index <- coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks %>% match(peakofinterest, .)
  return(gr_encode[coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks.granges.encode.overlaps %>% filter(queryHits == peakofinterest.index) %>% pull(subjectHits),] %>% as_tibble() %>% mutate(peak = peakofinterest))
}

# Import mouse narrowPeak files
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
qValue_cutoff <- 2
gr_narrowPeaks <- GRangesList()

for (i in setdiff(seq(0, 24), c(15, 24))) {
  file_narrowPeak <- paste0("newest_all_inclusive_analysis/MACS2/mouse/", i, ".coembed.subset.both.mouse.atac.bam_peaks.narrowPeak")
  gr_narrowPeak <- rtracklayer::import(file_narrowPeak, format = "BED", extraCols = extraCols_narrowPeak)
  gr_narrowPeak <- gr_narrowPeak[gr_narrowPeak$qValue > qValue_cutoff]
  gr_narrowPeaks <- c(gr_narrowPeaks, gr_narrowPeak)
}

# Reduce peaks
gr_narrowPeaks.reduced <- reduce(x = unlist(gr_narrowPeaks))
export.bed(gr_narrowPeaks.reduced, con='gr_narrowPeaks.mouse.reduced.bed')

# Reduce peaks with a minimum gap width of 500
gr_narrowPeaks.reduced.500 <- reduce(x = unlist(gr_narrowPeaks), min.gapwidth=500L)
export.bed(gr_narrowPeaks.reduced.500, con='gr_narrowPeaks.mouse.reduced.500.bed')

# Load and rename cells in the Seurat object
coembed.subset.both.mouse.atac <- readRDS("coembed.subset.both.mouse.atac.Rds")
mouse_fragment_path <- "atac.combined.fragments.tsv.gz"

# Create Fragment Object and Feature Matrix
mouse_frags <- CreateFragmentObject(
  path = mouse_fragment_path,
  cells = colnames(coembed.subset.both.mouse.atac),
  validate.fragments = TRUE
)

mouse_counts <- FeatureMatrix(
  fragments = mouse_frags,
  features = granges(gr_narrowPeaks.reduced.500),
  cells = colnames(coembed.subset.both.mouse.atac)
)

# Create Chromatin Assay and save the Seurat object
mouse_assay <- CreateChromatinAssay(mouse_counts, fragments = mouse_fragment_path)
coembed.subset.both.mouse.atac[['MACS2_reduced_500']] <- mouse_assay
DefaultAssay(coembed.subset.both.mouse.atac) <- "MACS2_reduced_500"

# Perform Differential Accessibility Analysis
coembed.subset.both.mouse.atac.DAR.AllMarkers <- FindAllMarkers(coembed.subset.both.mouse.atac, assay="MACS2_reduced_500", test.use = 'LR', latent.vars = 'peak_region_fragments', logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)

# Save and load results
coembed.subset.both.mouse.atac.DAR.AllMarkers <- readRDS("coembed.subset.both.mouse.atac.DAR.MACS2_reduced_500.AllMarkers.Rds")
coembed.subset.both.mouse.atac.DAR.AllMarkers %>% write_excel_csv("coembed.subset.both.mouse.atac.csv")

# Human ATAC side common and divergent analysis
coembed.subset.both.mouse.atac.DAR.AllMarkers.granges <- coembed.subset.both.mouse.atac.DAR.AllMarkers %>% filter(`p_val_adj` < 0.05) %>% pull(gene) %>% unique() %>% StringToGRanges()
mm10_hg38 <- rtracklayer::import.chain("liftover/mm10ToHg38.over.chain")
coembed.subset.both.mouse.atac.DAR.AllMarkers.granges_hg38 <- rtracklayer::liftOver(x = coembed.subset.both.mouse.atac.DAR.AllMarkers.granges, chain = mm10_hg38)
coembed.subset.both.mouse.atac.DAR.AllMarkers.granges_hg38.merge_2000 <- endoapply(coembed.subset.both.mouse.atac.DAR.AllMarkers.granges_hg38, FUN=reduce_2000L)

# Find overlaps between human and mouse peaks
coembed.subset.both.human.atac.DAR.AllMarkers <- readRDS('coembed.subset.both.human.atac.DAR.AllMarkers.merged.Rds')
coembed.subset.both.human.atac.DAR.AllMarkers_hg38 <- coembed.subset.both.human.atac.DAR.AllMarkers %>% filter(`p_val_adj` < 0.05)%>% pull(gene) %>% unique() %>% StringToGRanges()

findoverlaps.human_peaks_hg38.mouse_peaks_hg38 <- findOverlaps(coembed.subset.both.human.atac.DAR.AllMarkers_hg38, coembed.subset.both.mouse.atac.DAR.AllMarkers.granges_hg38.merge_2000.granges.extend_2kb, select="all")
human_common_peaks_hg38 <- coembed.subset.both.human.atac.DAR.AllMarkers_hg38[queryHits(findoverlaps.human_peaks_hg38.mouse_peaks_hg38)] %>% as.data.frame()

# Save common peaks
coembed.subset.both.human.atac.DAR.AllMarkers.common <- coembed.subset.both.human.atac.DAR.AllMarkers %>% filter(`gene` %in% (coembed.subset.both.human.atac.DAR.AllMarkers_hg38.intersection.coords %>% pull(hg38_coords)))
coembed.subset.both.human.atac.DAR.AllMarkers.common %>% filter(avg_logFC > 0.5 & `p_val_adj` < 0.05) %>% pull(gene) %>% unique() %>% StringToGRanges() %>% export.bed(con='coembed.subset.both.human.atac.DAR.AllMarkers.common.logFC.0.5.bed')

# Save human-specific peaks
coembed.subset.both.human.atac.DAR.AllMarkers.human_only <- coembed.subset.both.human.atac.DAR.AllMarkers %>% filter(!(`gene` %in% (coembed.subset.both.human.atac.DAR.AllMarkers_hg38.intersection.coords %>% pull(hg38_coords))))
coembed.subset.both.human.atac.DAR.AllMarkers.human_only %>% filter(avg_logFC > 0.5 & `p_val_adj` < 0.05) %>% pull(gene) %>% unique() %>% StringToGRanges() %>% export.bed(con='coembed.subset.both.human.atac.DAR.AllMarkers.human_only.logFC.0.5.bed')
