# src/functions/dimred.R
# ==============================================================================
# DIMENSIONALITY REDUCTION MODULE
# Description: Wrapper for UMAP generation and plotting (Pre/Post correction QC).
# Dependencies: scater, ggplot2, SingleCellExperiment
# ==============================================================================

library(scater)
library(ggplot2)
library(SingleCellExperiment)

#' @title Run UMAP
#' @description Calculates UMAP embedding on the SCE object with built-in downsampling for stability.
#' @param sce SingleCellExperiment object.
#' @param markers Character vector of markers to use for UMAP.
#' @param config List containing UMAP params (n_neighbors, min_dist, seed).
#' @return Subsampled SCE object with "UMAP" in reducedDims.
run_umap_generation <- function(sce, markers, config) {
  
  message("[DimRed] Preparing UMAP calculation...")
  set.seed(config$seed)
  
  if(!all(markers %in% rownames(sce))) stop("[DimRed Fatal] Markers for UMAP missing in SCE.")
  
  # Downsampling to prevent RAM exhaustion and visual overcrowding
  max_cells <- if (!is.null(config$max_cells)) as.integer(config$max_cells) else 100000
  total_cells <- ncol(sce)
  
  if (total_cells > max_cells) {
    message(sprintf("         Subsampling %d out of %d cells for UMAP visualization...", max_cells, total_cells))
    
    batches <- as.character(sce$batch)
    unique_batches <- sort(unique(batches))
    cells_per_batch <- floor(max_cells / length(unique_batches))
    
    message(sprintf("         Applying stratified sampling: ~%d cells per batch...", cells_per_batch))
    
    # Stratified selection to ensure identical cell indices across Pre and Post runs
    sampled_indices_list <- lapply(unique_batches, function(b) {
      idx <- which(batches == b)
      # Ensure reproducibility by selecting up to the available cells in small batches
      sample(idx, min(length(idx), cells_per_batch), replace = FALSE)
    })
    
    sampled_indices <- unlist(sampled_indices_list)
    sce_umap <- sce[, sampled_indices]
  } else {
    sce_umap <- sce
  }
  
  sce_umap <- scater::runUMAP(
    sce_umap,
    subset_row = markers,
    exprs_values = "exprs",
    n_neighbors = config$n_neighbors,
    min_dist = config$min_dist,
    name = "UMAP"
  )
  
  return(sce_umap)
}

#' @title Plot UMAP
#' @description Generates a static UMAP plot colored by a variable.
#' @param sce SCE object with UMAP computed.
#' @param color_by String, column name in colData to color by.
#' @param out_path String, path to save the PDF.
plot_umap <- function(sce, color_by, out_path) {
  
  message(sprintf("[DimRed] Plotting UMAP colored by %s...", color_by))
  
  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = color_by) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste("UMAP colored by", color_by))
  
  ggplot2::ggsave(out_path, plot = p, width = 8, height = 6)
}

#' @title Plot UMAP Split by Variable
#' @description Generates a faceted UMAP plot to evaluate batch correction.
#' @param sce SCE object with UMAP computed.
#' @param color_by String, column name in colData to color by.
#' @param split_by String, column name in colData to facet by.
#' @param out_path String, path to save the PDF.
plot_umap_split <- function(sce, color_by, split_by, out_path) {
  
  message(sprintf("[DimRed] Plotting UMAP colored by %s and split by %s...", color_by, split_by))
  
  # 'other_fields' is strictly required to pass the facet variable to scater's internal dataframe
  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = color_by, other_fields = split_by) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(as.formula(paste("~", split_by))) +
    ggplot2::ggtitle(paste("UMAP colored by", color_by, "| Split by", split_by))
  
  # Extended width to prevent squished panels in PDF
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 6)
}