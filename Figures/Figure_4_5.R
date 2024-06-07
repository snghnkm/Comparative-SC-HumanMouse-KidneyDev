# Figure 4.R

# Figure 4B-C

dotplot_features = c("EYA1", "BMPER", "WNT4", "LHX1", "MAFB", "NPHS2", "AKAP12", "LRP2", "EPHA7", "LGR5", "TFAP2B", "SLC12A1", "TMEM52B", "RET", "AQP2", "UPK3A", "FOXI1")
length(dotplot_features)
p1 <- DotPlot(coembed.subset.both.human.rna, features=rev(dotplot_features), cols = c("lightgrey", "red"), col.min = 0) + coord_flip()
p2 <- DotPlot(coembed.subset.both.mouse.rna, features=rev(dotplot_features), cols = c("lightgrey", "blue"), col.min = 0) + coord_flip()
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.RNA_norm.dotplot.markers.divided.eps", width=6*4*1.3, height=6*1.8*1.3)
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.RNA_norm.dotplot.markers.divided.pdf", width=6*4*1.3, height=6*1.8*1.3)

p1 <- DotPlot(coembed.subset.both.human.atac, features=rev(dotplot_features), cols = c("lightgrey", "red"), col.min = 0) + coord_flip()
p2 <- DotPlot(coembed.subset.both.mouse.atac, features=rev(dotplot_features), cols = c("lightgrey", "blue"), col.min = 0) + coord_flip()
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.atac_norm.dotplot.markers.divided.eps", width=6*4*1.3, height=6*1.8*1.3)
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.atac_norm.dotplot.markers.divided.pdf", width=6*4*1.3, height=6*1.8*1.3)


dotplot_features2 = c(
"PCDH15", "SLC9C1", "NRXN1", "ROR2", "MAP2", "SYT1", "FANCA", "HES4", "KLHL29", "TARID", "PTPRQ", "PLA2R1", "NTNG1", "NPHS2", "FRMD1", "DPYS", "TRPM3", "ACPP", "CLCN5", "INPP4B", "SLC8A1", "NMU", "KCTD8", "CACNA2D3", "PCDH7", "PIP5K1B",  "PIK3C2G", "CEMIP2", "CD96", "DHRS2", "BCAS1", "COLCA1", "ZNF66")
length(dotplot_features2)
p1 <- DotPlot(coembed.subset.both.human.rna, features=rev(dotplot_features2), cols = c("lightgrey", "red"), col.min = 0) + coord_flip()
p2 <- DotPlot(coembed.subset.both.mouse.rna, features=rev(dotplot_features2), cols = c("lightgrey", "blue"), col.min = 0) + coord_flip()
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.RNA_norm.dotplot2.markers.divided.eps", width=6*4*1.3, height=6*1.8*1.3)
ggsave(p1|p2, filename="newest_all_inclusive_analysis.N_UE.RNA_norm.dotplot2.markers.divided.pdf", width=6*4*1.3, height=6*1.8*1.3)

# Fig 4D.

human_GO_Biological_Process <- read_csv("newest_all_inclusive_analysis/DEG_GO_fianl.csv")
human_GO_Biological_Process$Term <- factor(human_GO_Biological_Process$Term, levels = rev(unique(human_GO_Biological_Process$Term)))

ggsave(ggplot(human_GO_Biological_Process, aes(x = Term, y = minuslogpvalue, fill = Group)) + 
         geom_bar(stat = "identity", color = "black", position = position_dodge()) + 
         coord_flip() + 
         ggpubr::theme_pubr(), 
filename = "human_GO_Biological_Process.pdf", width = 15, height = 7)

# Figure 4E-I

featureplot_vlnplot_combined <- function(my_feature){
    p1 <- FeaturePlot(coembed.subset.both.human.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, pt.size=0.5, blend=F, order=T)
    p2 <- VlnPlot(coembed.subset.both.human.rna, features=my_feature, pt.size=0) + NoLegend()
    p3 <- FeaturePlot(coembed.subset.both.mouse.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, , pt.size=0.5, blend=F, order=T)
    p4 <- VlnPlot(coembed.subset.both.mouse.rna, features=my_feature, pt.size=0) + NoLegend()
    ggsave(
        (p1 / p2 / p3 / p4 ) + plot_layout(heights=c(18,3,18,3)),
        filename=paste0("featureplot.vlnplot.combined.",my_feature,".2.png"), width=12, height=30, limitsize = FALSE
    )
    ggsave(
        (p1 / p2 / p3 / p4 ) + plot_layout(heights=c(18,3,18,3)),
        filename=paste0("featureplot.vlnplot.combined.",my_feature,".2.pdf"), width=12, height=30, limitsize = FALSE
    )
}

featureplot_vlnplot_combined("NRXN1")
featureplot_vlnplot_combined("NTNG1")
featureplot_vlnplot_combined("PLA2R1")
featureplot_vlnplot_combined("PLXNA2")
featureplot_vlnplot_combined("LMX1B")

# CAKUT genes
featureplot_vlnplot_combined("AHI1")
featureplot_vlnplot_combined("CD96")
featureplot_vlnplot_combined("GLI2")
featureplot_vlnplot_combined("GLI3")
featureplot_vlnplot_combined("GREB1L")
featureplot_vlnplot_combined("KANK4")
featureplot_vlnplot_combined("MID1")
featureplot_vlnplot_combined("PCSK5")
featureplot_vlnplot_combined("ROR2")
featureplot_vlnplot_combined("SEMA3A")
featureplot_vlnplot_combined("TRPC6")

# human-specific DEGs and DARs
featureplot_vlnplot_combined("PTPRQ")
featureplot_vlnplot_combined("DLG2")
featureplot_vlnplot_combined("NKAIN3")

# conserved DEGs and DARs
featureplot_vlnplot_combined("BMPER")
featureplot_vlnplot_combined("LMX1B")
featureplot_vlnplot_combined("CLIC5")
featureplot_vlnplot_combined("SLC12A2")
featureplot_vlnplot_combined("ITGB6")
featureplot_vlnplot_combined("RXRA")
featureplot_vlnplot_combined("COL5A1")
featureplot_vlnplot_combined("CTIF")
featureplot_vlnplot_combined("ZBTB7C")
featureplot_vlnplot_combined("SLC39A11")

featureplot_vlnplot_combined <- function(my_feature){
    p1 <- FeaturePlot(coembed.subset.both.human.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, pt.size=0.5, blend=F, order=T)
    p2 <- VlnPlot(coembed.subset.both.human.rna, features=my_feature, pt.size=0) + NoLegend()
    p3 <- FeaturePlot(coembed.subset.both.mouse.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, , pt.size=0.5, blend=F, order=T)
    p4 <- VlnPlot(coembed.subset.both.mouse.rna, features=my_feature, pt.size=0) + NoLegend()
    ggsave(
        (p1 / p2 / p3 / p4 ) + plot_layout(heights=c(18,3,18,3)),
        filename=paste0("./human_specific_DEGs/featureplot.vlnplot.combined.",my_feature,".2.png"), width=12, height=30, limitsize = FALSE
    )
}

diff_local_enriched_genes <- setdiff(intersect(full_join.human.vs.mouse.DEGs.human_specific.geneTable %>% mutate(`pct_diff` = `pct.1.x` - `pct.2.x`) %>% filter(`avg_logFC.x` > 0.5 & `pct_diff` > 0.3) %>% pull(`gene.x`) %>% unique(), coembed.subset.both.mouse.rna.DEG.AllMarkers.cluster_gene %>% pull(`gene`) %>% unique()), conservedmarkers %>% pull(rowname) %>% unique() %>% sort())

for (each_human_specific_gene in diff_local_enriched_genes){
  featureplot_vlnplot_combined(each_human_specific_gene)
}

diff_local_enriched_genes <- read_csv("diff_local_enriched_genes.full_join.csv") %>% pull(`gene.x`) %>% unique()

DefaultAssay(coembed.subset.both.human.rna) <- "RNA"
DefaultAssay(coembed.subset.both.mouse.rna) <- "RNA"

featureplot_vlnplot_combined <- function(my_feature){
    p1 <- FeaturePlot(coembed.subset.both.human.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, pt.size=0.5, blend=F, order=T)
    p2 <- VlnPlot(coembed.subset.both.human.rna, features=my_feature, pt.size=0) + NoLegend()
    p3 <- FeaturePlot(coembed.subset.both.mouse.rna, features=c(my_feature), min.cutoff="q10", max.cutoff="q90", label=T, repel=T, raster=T, , pt.size=0.5, blend=F, order=T)
    p4 <- VlnPlot(coembed.subset.both.mouse.rna, features=my_feature, pt.size=0) + NoLegend()
    ggsave(
        (p1 / p2 / p3 / p4 ) + plot_layout(heights=c(18,3,18,3)),
        filename=paste0("./diff_local_enrichment/featureplot.vlnplot.combined.",my_feature,".3.pdf"), width=12, height=30, limitsize = FALSE
    )
}

for (each_human_specific_gene in diff_local_enriched_genes){
  tryCatch({
    featureplot_vlnplot_combined(each_human_specific_gene)
  }, error=function(e){})
}

