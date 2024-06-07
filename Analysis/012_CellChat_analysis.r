# 012_CellChat_analysis.r
# Load required libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(CellChat)

# Process Nephrogenic subset
Nephrogenic <- subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, ident=c(4,2,14,25,6,11))
Nephrogenic <- RenameIdents(object = Nephrogenic, `4` = "NPC", `2` = "NPC", `14` = "NPC", `25` = "UPC", `6` = "UPC", `11` = "IPC")
Nephrogenic[['active.ident']] <- Nephrogenic@active.ident
Nephrogenic <- createCellChat(object = Nephrogenic, group.by = "active.ident")

# Process Glomerulus subsets
Glomerulus_Pod <- subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, ident=c(17,0))
Glomerulus_gEC <- subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, ident=c(13))
Glomerulus_gEC_EHD3 <- subset(Glomerulus_gEC, EHD3 > 0, invert = TRUE)
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, cells.highlight = WhichCells(Glomerulus_gEC_EHD3), label = TRUE, repel = TRUE), filename = "Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.Glomerulus_gEC_EHD3.pdf")
Glomerulus_Mes <- subset(Hp0_10_SCT.combined.rm022.rm3134.rm303132, ident=c(28))
Glomerulus_Mes_noREN <- subset(Glomerulus_Mes, REN > 0, invert = TRUE)
ggsave(UMAPPlot(Hp0_10_SCT.combined.rm022.rm3134.rm303132, cells.highlight = WhichCells(Glomerulus_Mes_noREN), label = TRUE, repel = TRUE), filename = "Hp0_10_SCT.combined.rm022.rm3134.rm303132.UMAP.Glomerulus_Mes_noREN.pdf")

# Merge and rename Glomerulus subsets
Glomerulus <- merge(x = Glomerulus_Pod, y = c(Glomerulus_gEC_EHD3, Glomerulus_Mes_noREN))
Glomerulus <- RenameIdents(object = Glomerulus, `17` = "Pod", `0` = "Pod", `28` = "Mes", `13` = "gEC")
table(Glomerulus@active.ident)
Glomerulus[['active.ident']] <- Glomerulus@active.ident
Glomerulus <- createCellChat(object = Glomerulus, group.by = "active.ident")

# Load CellChat database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# Process Nephrogenic interactions
Nephrogenic@DB <- CellChatDB.use
Nephrogenic <- subsetData(Nephrogenic)
Nephrogenic <- identifyOverExpressedGenes(Nephrogenic)
Nephrogenic <- identifyOverExpressedInteractions(Nephrogenic)
Nephrogenic <- computeCommunProb(Nephrogenic, type = "triMean")
Nephrogenic <- filterCommunication(Nephrogenic, min.cells = 10)
Nephrogenic <- computeCommunProbPathway(Nephrogenic)
Nephrogenic <- aggregateNet(Nephrogenic)
Nephrogenic <- netAnalysis_computeCentrality(Nephrogenic, slot.name = "netP")

# Process Glomerulus interactions
Glomerulus@DB <- CellChatDB.use
Glomerulus <- subsetData(Glomerulus)
Glomerulus <- identifyOverExpressedGenes(Glomerulus)
Glomerulus <- identifyOverExpressedInteractions(Glomerulus)
Glomerulus <- computeCommunProb(Glomerulus, type = "triMean")
Glomerulus <- filterCommunication(Glomerulus, min.cells = 10)
Glomerulus <- computeCommunProbPathway(Glomerulus)
Glomerulus <- aggregateNet(Glomerulus)
Glomerulus <- netAnalysis_computeCentrality(Glomerulus, slot.name = "netP")

# Visualize Nephrogenic signaling
ht1 <- netAnalysis_signalingRole_heatmap(Nephrogenic, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(Nephrogenic, pattern = "incoming")
ht1 + ht2
netAnalysis_signalingRole_network(Nephrogenic, signaling = "NRXN", width = 8, height = 2.5, font.size = 10)
netVisual_bubble(Nephrogenic, sources.use = c(1, 2, 3), targets.use = c(1, 2, 3), remove.isolate = FALSE, thresh = 0.01)
plotGeneExpression(Nephrogenic, signaling = "NRXN")

# Visualize Glomerulus signaling
ht1 <- netAnalysis_signalingRole_heatmap(Glomerulus, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(Glomerulus, pattern = "incoming")
ht1 + ht2
netAnalysis_signalingRole_network(Glomerulus, signaling = "NRG", width = 8, height = 2.5, font.size = 10)
netVisual_bubble(Nephrogenic, sources.use = c(1, 2, 3), targets.use = c(1, 2, 3), remove.isolate = FALSE, thresh = 0.01)

# Search for specific ligand-receptor pairs in WNT signaling
pairLR <- searchPair(signaling = "WNT", pairLR.use = CellChatDB.human$interaction, key = "pathway_name", matching.exact = TRUE, pair.only = FALSE)