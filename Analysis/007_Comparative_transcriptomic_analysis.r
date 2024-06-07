# 007_Comparative_transcriptomic_analysis.r

# Load required libraries
lapply(c("Signac", "Seurat", "ggplot2", "patchwork", "readr", "tidyr", "dplyr", "tibble", "rtracklayer", "stringr"), library, character.only = TRUE)

# Mouse RNA Analysis
coembed.subset.both.mouse.rna.DEG.AllMarkers <- coembed.subset.both.mouse.rna %>% FindAllMarkers(assay="RNA", only.pos = TRUE)

# Human RNA Analysis
coembed.subset.both.human.rna.DEG.AllMarkers <- coembed.subset.both.human.rna %>% FindAllMarkers(assay="RNA", only.pos = TRUE)

# Filter human RNA markers based on adjusted p-value
num_unique_genes <- coembed.subset.both.human.rna.DEG.AllMarkers %>% filter(p_val_adj < 0.05) %>% pull(gene) %>% unique() %>% length()
print(num_unique_genes) # Expected output: 5426

# Subsetting RNA data for both human and mouse
coembed.subset.both.rna <- subset(coembed.subset.both, subset = tech == "RNA")

# Finding conserved markers between human and mouse RNA data
conservedmarkers <- data.frame()
clusters <- c(0, 12, 17, 10, 2, 16, 13, 19, 15, 1, 11, 7, 8, 18, 9, 14, 20, 6, 22, 4, 24, 5, 3, 21, 23)
for (each in clusters) {
  res <- FindConservedMarkers(coembed.subset.both.rna, grouping.var = "species", ident.1 = each, only.pos = TRUE)
  res$cluster <- each
  conservedmarkers <- rbind(conservedmarkers, res %>% rownames_to_column())
}

# Calculate percentage differences and filter conserved markers
conservedmarkers.pct_diff <- conservedmarkers %>% 
  mutate(`human_pct_diff` = `human_pct.1` - `human_pct.2`) %>% 
  mutate(`mouse_pct_diff` = `mouse_pct.1` - `mouse_pct.2`)

conservedmarkers.pct_diff.filtered <- conservedmarkers.pct_diff %>% 
  filter(mouse_avg_logFC > 0.5 & human_avg_logFC > 0.5 & human_pct_diff > 0.3 & mouse_pct_diff > 0.3)

conservedmarkers.pct_diff.filtered.2 <- conservedmarkers.pct_diff %>% 
  filter(mouse_avg_logFC > 0.25 & human_avg_logFC > 0.25 & human_pct_diff > 0.25 & mouse_pct_diff > 0.25)


# Specific gene expression for human or mouse
human_specific_genes <- intersect(
  coembed.subset.both.human.rna.DEG.AllMarkers %>% 
    mutate(`pct_diff` = `pct.1` - `pct.2`) %>% 
    filter(avg_logFC > 0.5 & pct_diff > 0.3 & p_val_adj < 0.05) %>% 
    pull(gene) %>% unique(),
  geneTable %>% pull(`Human.Symbol`) %>% unique()
)

# Joining human and mouse DEGs
conservedmarkers.cluster_gene <- conservedmarkers %>% unite("cluster_gene", c("cluster", "rowname"), remove = FALSE)
coembed.subset.both.mouse.rna.DEG.AllMarkers.cluster_gene <- coembed.subset.both.mouse.rna.DEG.AllMarkers %>% unite("cluster_gene", c("cluster", "gene"), remove = FALSE)
coembed.subset.both.human.rna.DEG.AllMarkers.cluster_gene <- coembed.subset.both.human.rna.DEG.AllMarkers %>% unite("cluster_gene", c("cluster", "gene"), remove = FALSE)

# Full join to compare human vs mouse DEGs
full_join.human.vs.mouse.DEGs <- full_join(
  full_join(coembed.subset.both.human.rna.DEG.AllMarkers.cluster_gene, conservedmarkers.cluster_gene, by = "cluster_gene"),
  coembed.subset.both.mouse.rna.DEG.AllMarkers.cluster_gene, by = "cluster_gene"
)

# Filter for human-specific genes
full_join.human.vs.mouse.DEGs.human_specific <- full_join.human.vs.mouse.DEGs %>% 
  filter(`cluster.x` %in% clusters) %>% 
  filter(is.na(`cluster`)) %>% 
  mutate(`pct_diff` = `pct.1.x` - `pct.2.x`) %>% 
  filter(`avg_logFC.x` > 0.5 & `pct_diff` > 0.3)

full_join.human.vs.mouse.DEGs.human_specific.2 <- full_join.human.vs.mouse.DEGs %>% 
  filter(`cluster.x` %in% clusters) %>% 
  filter(is.na(`cluster`)) %>% 
  mutate(`pct_diff` = `pct.1.x` - `pct.2.x`) %>% 
  filter(`p_val_adj.x` < 0.05)

# Further filtering and analysis
final_human_specific_genes <- intersect(
  setdiff(full_join.human.vs.mouse.DEGs.human_specific.geneTable %>% pull(`gene.x`) %>% unique(), 
          coembed.subset.both.mouse.rna.DEG.AllMarkers.cluster_gene %>% pull(`gene`) %>% unique()),
  geneTable %>% pull(`Human.Symbol`) %>% unique()
)

# Annotations and GO Biological Process analysis
GREAT_annotation.human_peaks.divergent.gene2regions <- read_csv("GREAT_annotation.human_peaks.divergent.gene2regions.csv")

GREAT_annotated_divergent_genes <- intersect(final_human_specific_genes, GREAT_annotation.human_peaks.divergent.gene2regions %>% pull(`# GREAT version 4.0.4`)) %>% tibble()
GREAT_annotated_divergent_genes$divergent_peak <- "yes"

human_specific_GO_Biological_Process <- read_tsv("newest_all_inclusive_analysis/Human_specific_GO_Biological_Process_2021_table (2).txt") %>% 
  mutate(`minuslogpvalue` = -log10(`P-value`)) %>% 
  select(c(`Term`, `P-value`, `minuslogpvalue`))

human_common_GO_Biological_Process <- read_tsv("newest_all_inclusive_analysis/Human_common_GO_Biological_Process_2021_table (2).txt") %>% 
  mutate(`minuslogpvalue` = -log10(`P-value`)) %>% 
  select(c(`Term`, `P-value`, `minuslogpvalue`))

full_join(human_specific_GO_Biological_Process %>% filter(`P-value` < 0.05), human_common_GO_Biological_Process %>% filter(`P-value` < 0.05), by = "Term") %>% 
  write_excel_csv("full_join.human_GO.csv")
