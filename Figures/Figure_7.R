# Figure_7.R

# Fig 7B
# Merge DORC results and annotate
DORC_res_merged <- readRDS("newest_all_inclusive_analysis/DORC/DORC_res_merged.Rds")

# Prepare data for volcano plot
DORC_data.significant.volcano <- DORC_res_merged %>%
  unite("gene_peak", c(gene, peak), remove = FALSE) %>%
  mutate(`neg_log10_norm_obs_p` = -log10(norm_obs_p_pearson)) 

# Select top and bottom significant gene-peaks
DORC_data.significant.top50 <- top_n(DORC_data.significant.volcano, wt=`norm_obs_pearson`, n=100) %>% pull(`gene_peak`)
DORC_data.significant.bottom50 <- top_n(DORC_data.significant.volcano, wt=-`norm_obs_pearson`, n=100) %>% pull(`gene_peak`)

DORC_data.significant.volcano <- DORC_data.significant.volcano %>% 
  mutate(label = case_when(
    `gene_peak` %in% DORC_data.significant.top50 ~ gene,
    `gene_peak` %in% DORC_data.significant.bottom50 ~ gene,
    TRUE ~ ""
  ))

# Generate volcano plot
volcano_plot <- ggplot(DORC_data.significant.volcano, aes(obs_pearson, neg_log10_norm_obs_p, label = label)) +
  geom_point(alpha=0.4) +
  ggrepel::geom_text_repel(size=3) +
  theme_pubr() +
  labs(title = "Volcano Plot of DORCs", x = "Observed Pearson Correlation", y = "-log10 Normalized Observed p-value")

# Save volcano plot
ggsave(volcano_plot, filename="DORC_data.significant.volcano.pdf", width=12, height=12)

# Fig G-H.
featureplot_vlnplot_combined <- function(my_feature){
    p1 <- FeaturePlot(coembed.subset.both.human.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, pt.size=0.5, blend=F, order=T)
    p2 <- VlnPlot(coembed.subset.both.human.rna, features=my_feature, pt.size=0) + NoLegend()
    p3 <- FeaturePlot(coembed.subset.both.mouse.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, , pt.size=0.5, blend=F, order=T)
    p4 <- VlnPlot(coembed.subset.both.mouse.rna, features=my_feature, pt.size=0) + NoLegend()
    ggsave(
        (p1 / p2 / p3 / p4 ) + plot_layout(heights=c(18,3,18,3)),
        filename=paste0("./GWAS_final_genesets/featureplot.vlnplot.combined.",my_feature,".3.pdf"), width=12, height=30, limitsize = FALSE
    )

}

GWAS_kidney_related_final_genesets <- c("CDA", "NR", "FADS2", "MIR1908", "WIPF3", "DPY19L2P3", "PLA2R1", "LOC100132354", "LINC01512", "VEGFA", "WDR81", "MYRF", "HAS2", "MFSD4A", "SLC16A9", "HOXD8", "STC1", "GBAP1", "TRIM46", "MPPED2", "PSD4", "AK123617", "PAX8", "FBN3", "CERS4", "HNF1B", "SIM1", "PHYHD1", "FMO4", "LOC101928126", "ATP50", "LOC105372790", "WWC1", "FGFR2", "CLDN14", "PDILT", "UMOD", "ACSM5", "GP2", "IRX1", "TFEB", "RP11-429K17.1", "KIRREL2", "SLC7A9", "DAB2", "LRP2", "SLC47A1", "GATM", "LOC154092", "ALPL", "DPEP1", "BCAS3", "TBX2", "C17orf82", "A4GALT", "PCK1", "KIAA1614", "SLC34A1", "ZEB2", "RGS14", "PFN3", "F12", "TARID", "FADS1", "LINC02537", "TMEM258", "RNU2-19P", "HOXD-AS2", "SPIRE2", "MUC1", "MPPED2-AS1", "PAX8-AS1", "SRP14P4", "GM2AP2", "LINC00649", "MRPS6", "LINC01153", "KNG1", "SMAD7", "PPFIA2", "NPHS1", "C9", "CT69", "LINC01010", "ACSM2A", "PDLIM4", "P4HA2", "LINC02068", "ABO", "LNC-LBCS", "LINC02413")

featureplot_vlnplot_combined("GATM")
for (each_human_specific_gene in GWAS_kidney_related_final_genesets){
  tryCatch({
    featureplot_vlnplot_combined(each_human_specific_gene)
  }, error=function(e){})
}


