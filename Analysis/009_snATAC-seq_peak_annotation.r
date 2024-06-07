# 009_snATAC-seq_peak_annotation.r

# Load required libraries
library(tidyverse)
library(Signac)
library(Seurat)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)

# Set working directory and initialize
closeAllConnections()
gc()

# Filenames for annotation
filenames <- c("newest_all_inclusive_analysis/DORC/DORC_res_merged.csv")

# Configuration
cluster_n <- detectCores()
tssRegion_window <- c(-2000, 0)
flankDistance_max <- 5 * 10^6 # 5M flanking distance
n_return <- 10
coord_column_name <- "peak"

# Annotate peaks for each file
for (each_filename in filenames) {
    cl <- parallel::makeCluster(cluster_n)
    doSNOW::registerDoSNOW(cl)
    
    Peaks_csv_filename <- each_filename
    output_annotated_rds <- paste0(Peaks_csv_filename, ".annotated.rds")
    output_annotated_csv <- paste0(Peaks_csv_filename, ".annotated.csv")
    
    message("Reading file: ", each_filename)
    Peaks <- read_csv(Peaks_csv_filename)
    Peaks.coord.granges <- StringToGRanges(Peaks %>% pull(as.name(coord_column_name)))
    
    # Annotate peaks
    message("Annotating peaks...")
    granges <- Peaks.coord.granges
    iterations <- length(granges)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    gc()
    
    Peaks.coord.granges.annotated <- foreach(i = 1:iterations, .combine = dplyr::bind_rows, .options.snow = opts) %dopar% {
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        library(org.Hs.eg.db)
        
        annotated <- tibble::as_tibble(
            ChIPseeker::annotatePeak(
                granges[i],
                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                annoDb = "org.Hs.eg.db",
                addFlankGeneInfo = TRUE,
                flankDistance = flankDistance_max, 
                tssRegion = tssRegion_window,
                level = 'transcript',
                overlap = "TSS",
                verbose = FALSE
            )
        )
        return(annotated)
    }
    close(pb)
    
    # Get nearest up/downstream genes
    message("Getting nearest up/downstream genes...")
    annotatedPeaks <- Peaks.coord.granges.annotated
    
    Peaks.coord.granges.annotated.updown <- foreach(i = seq(nrow(annotatedPeaks)), .combine = dplyr::bind_rows, .options.snow = opts) %dopar% {
        library(dplyr)
        library(tibble)
        library(org.Hs.eg.db)
        
        target <- annotatedPeaks[i, ]
        x <- tibble::tibble(
            flank_txIds = (target$flank_txIds %>% stringr::str_split(";"))[[1]],
            flank_geneIds = (target$flank_geneIds %>% stringr::str_split(";"))[[1]] %>% as.numeric() %>% as.character(),
            flank_gene_distances = (target$flank_gene_distances %>% stringr::str_split(";"))[[1]] %>% as.numeric()
        )
        
        y <- suppressMessages({
            AnnotationDbi::select(
                org.Hs.eg.db, 
                keys = x$flank_geneIds,
                columns = c("ENTREZID", "SYMBOL", "GENENAME"),
                keytype = "ENTREZID"
            )
        }) %>% unique()
        
        joined <- dplyr::left_join(x, y, by = c("flank_geneIds" = "ENTREZID"))
        joined <- joined %>% dplyr::select(-flank_txIds) %>% unique()
        
        joined.flank_gene_distance.top1.arranged <- joined %>% group_by(flank_geneIds) %>%
            top_n(n = 1, wt = -abs(flank_gene_distances)) %>% arrange(abs(flank_gene_distances)) %>% ungroup()
        
        joined.flank_gene_distance.top1.arranged.upstream <- joined.flank_gene_distance.top1.arranged %>% 
            filter(flank_gene_distances < 0) %>% top_n(n = n_return, wt = flank_gene_distances) %>% arrange(-flank_gene_distances)
        
        joined.flank_gene_distance.top1.arranged.downstream <- joined.flank_gene_distance.top1.arranged %>% 
            filter(flank_gene_distances >= 0) %>% top_n(n = n_return, wt = -flank_gene_distances) %>% arrange(flank_gene_distances)
        
        tmp1 <- annotatedPeaks[i, ] %>% dplyr::select(-c(transcriptId, flank_txIds, flank_geneIds, flank_gene_distances))
        tmp2 <- tibble::tibble(
            nearest_upstream_gene = joined.flank_gene_distance.top1.arranged.upstream %>% pull(SYMBOL) %>% toString(),
            nearest_upstream_gene_distance = joined.flank_gene_distance.top1.arranged.upstream %>% pull(flank_gene_distances) %>% toString(),
            nearest_upstream_geneID = joined.flank_gene_distance.top1.arranged.upstream %>% pull(flank_geneIds) %>% toString(),
            nearest_downstream_gene = joined.flank_gene_distance.top1.arranged.downstream %>% pull(SYMBOL) %>% toString(),
            nearest_downstream_gene_distance = joined.flank_gene_distance.top1.arranged.downstream %>% pull(flank_gene_distances) %>% toString(),
            nearest_downstream_geneID = joined.flank_gene_distance.top1.arranged.downstream %>% pull(flank_geneIds) %>% toString()
        )
        
        tmp <- dplyr::bind_cols(tmp1, tmp2)
        return(tmp)
    }
    close(pb)
    
    # Merge and save results
    Peaks.final <- dplyr::bind_cols(Peaks, Peaks.coord.granges.annotated.updown %>% dplyr::select(-c(seqnames, start, end, strand)))
    
    message("Saving annotated results...")
    Peaks.final %>% saveRDS(output_annotated_rds)
    Peaks.final %>% write_excel_csv(output_annotated_csv)
}

message("Annotation completed.")
