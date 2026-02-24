# src/functions/dimred.R
# ==============================================================================
# DIMENSIONALITY REDUCTION MODULE
# Description: Wrapper for UMAP generation and plotting.
# Dependencies: scater, ggplot2, SingleCellExperiment
# ==============================================================================

library(scater)
library(ggplot2)
library(SingleCellExperiment)

#' @title Run UMAP
#' @description Calculates UMAP embedding on the SCE object, applying stratified downsampling if necessary.
#' @param sce SingleCellExperiment object.
#' @param markers Character vector of markers to use for UMAP.
#' @param config List containing UMAP params (n_neighbors, min_dist, seed, max_cells).
#' @return SCE object (subsampled) with "UMAP" in reducedDims.
run_umap_generation <- function(sce, markers, config) {
  
  message("[DimRed] Preparing UMAP execution...")
  set.seed(config$seed)
  
  # 1. Validation
  if(!all(markers %in% rownames(sce))) stop("[DimRed Fatal] Markers for UMAP missing in SCE.")
  
  # 2. Downsampling Setup
  # If max_cells is not configured, safely default to 100,000 cells to prevent memory crash
  max_cells <- if (!is.null(config$max_cells)) as.integer(config$max_cells) else 100000
  total_cells <- ncol(sce)
  
  if (total_cells > max_cells) {
    message(sprintf("[DimRed] Subsampling %d out of %d cells for UMAP calculation...", max_cells, total_cells))
    message("         (Note: This downsampled SCE is for visualization only. Statistics rely on the full clustered_sce.rds)")
    
    # Random sampling of columns (cells)
    sampled_indices <- sample(seq_len(total_cells), max_cells, replace = FALSE)
    sce_umap <- sce[, sampled_indices]
  } else {
    message(sprintf("[DimRed] Total cells (%d) is within limit (%d). Proceeding with all cells.", total_cells, max_cells))
    sce_umap <- sce
  }
  
  # 3. Execution
  message("[DimRed] Computing topologic manifold (this may take a while)...")
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
#' @param color_by String, column name in colData to color by (e.g., "batch", "metacluster_id").
#' @param out_path String, path to save the PDF/PNG.
plot_umap <- function(sce, color_by, out_path) {
  
  message(sprintf("[DimRed] Plotting UMAP colored by %s...", color_by))
  
  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = color_by) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste("UMAP colored by", color_by))
  
  ggplot2::ggsave(out_path, plot = p, width = 8, height = 6)
}

#' @title Plot UMAP Split by Variable
#' @description Generates a static UMAP plot colored by a variable and faceted by another.
#' @param sce SCE object with UMAP computed.
#' @param color_by String, column name in colData to color by (e.g., "metacluster_id").
#' @param split_by String, column name in colData to facet by (e.g., "batch").
#' @param out_path String, path to save the PDF.
plot_umap_split <- function(sce, color_by, split_by, out_path) {
  
  message(sprintf("[DimRed] Plotting UMAP colored by %s and split by %s...", color_by, split_by))
  
  # 'other_fields' is strictly required to pass the split variable to the internal plotting dataframe
  p <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = color_by, other_fields = split_by) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(as.formula(paste("~", split_by))) +
    ggplot2::ggtitle(paste("UMAP colored by", color_by, "split by", split_by))
  
  # Width is extended to 12 to accommodate multiple panels side-by-side without squishing
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 6)
}