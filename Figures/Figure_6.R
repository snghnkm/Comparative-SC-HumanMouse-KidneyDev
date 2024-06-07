# Figure_6.R

# Load required libraries
library(ggplot2)
library(ggpubr)
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)
library(dplyr)
library(rtracklayer)

# Figure 6A-B

# Load and filter data for plotting
DEG.common_divergent.distribution <- read_csv("DEG.common_divergent.distribution.csv")
DEG.common_divergent.distribution$group <- factor(DEG.common_divergent.distribution$group, levels = rev(c("Cons. DA", "Human-specific DA", "Both", "No linked DA")))

# Plot and save distribution of conserved DEGs
ggp_cons <- ggplot(DEG.common_divergent.distribution %>% filter(x == "Cons."), aes(x = x, y = y, fill = group, label = y)) + 
  geom_bar(stat = "identity") + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "left")
ggsave(ggp_cons, filename="DEG.common_divergent.distribution.cons.pdf", width=4, height=7)

# Plot and save distribution of human-specific DEGs
ggp_human <- ggplot(DEG.common_divergent.distribution %>% filter(x == "Human-specific"), aes(x = x, y = y, fill = group, label = y)) + 
  geom_bar(stat = "identity") + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "left")
ggsave(ggp_human, filename="DEG.common_divergent.distribution.Human-specific.pdf", width=4, height=7)

# Initialize parallel cluster
cl <- parallel::makeCluster(20)
doSNOW::registerDoSNOW(cl)

# Process human DAR with ENCODE annotation in parallel
human_DAR_with_encode_annotation <- foreach(i = seq(length(coembed.subset.both.human.atac.DAR.AllMarkers.p.0.05.peaks)), .combine = dplyr::bind_rows, .packages = c("dplyr", "rtracklayer")) %dopar% {
  peakofinterest <- coembed.subset.both.human.atac.DAR.AllMarkers.p.0.05.peaks[i]
  peakofinterest.index <- coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks %>% match(peakofinterest, .)
  overlaps <- findOverlaps(coembed.subset.both.human.atac.DAR.AllMarkers.unique.peaks.granges, gr_encode) %>% as.data.frame()
  gr_encode[overlaps %>% filter(queryHits == peakofinterest.index) %>% pull(subjectHits), ] %>% 
    as_tibble() %>% 
    mutate(peak = peakofinterest)
}

# Calculate and display proportions of human-specific DARs with ENCODE annotation
prop_table <- human_DAR_with_encode_annotation %>% 
  filter(peak %in% coembed.subset.both.human.atac.DAR.AllMarkers.human_only.p.0.05.peaks) %>% 
  pull(shortname) %>% 
  table() %>% 
  prop.table()
print(prop_table)

# Plot pie chart of human-specific DARs with ENCODE annotation
data <- data.frame(
  group = c("CTCF-only", "DNase-H3K4me3", "PLS", "dELS", "pELS"),
  value = c(909, 889, 5698, 45516, 16398)
)
p <- ggplot(data, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0)
ggsave(p, filename="human_DAR_with_encode_annotation.human_only.pdf")

# Figure 6C-F: Combined feature and violin plots
featureplot_vlnplot_combined <- function(my_feature) {
  p1 <- FeaturePlot(coembed.subset.both.human.rna, features = c(my_feature), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.5, blend = FALSE, order = TRUE)
  p2 <- VlnPlot(coembed.subset.both.human.rna, features = my_feature, pt.size = 0) + NoLegend()
  p3 <- FeaturePlot(coembed.subset.both.mouse.rna, features = c(my_feature), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.5, blend = FALSE, order = TRUE)
  p4 <- VlnPlot(coembed.subset.both.mouse.rna, features = my_feature, pt.size = 0) + NoLegend()
  
  combined_plot <- (p1 / p2 / p3 / p4) + plot_layout(heights = c(18, 3, 18, 3))
  
  ggsave(combined_plot, filename = paste0("featureplot.vlnplot.combined.", my_feature, ".2.png"), width = 12, height = 30, limitsize = FALSE)
  ggsave(combined_plot, filename = paste0("featureplot.vlnplot.combined.", my_feature, ".2.pdf"), width = 12, height = 30, limitsize = FALSE)
}

# Generate combined plots for selected features
featureplot_vlnplot_combined("PTPRQ")
featureplot_vlnplot_combined("DLG2")
featureplot_vlnplot_combined("LMX1B")
featureplot_vlnplot_combined("ZBTB7C")
