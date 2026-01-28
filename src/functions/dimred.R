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
#' @description Calculates UMAP embedding on the SCE object.
#' @param sce SingleCellExperiment object.
#' @param markers Character vector of markers to use for UMAP.
#' @param config List containing UMAP params (n_neighbors, min_dist, seed).
#' @return SCE object with "UMAP" in reducedDims.
run_umap_generation <- function(sce, markers, config) {
  
  message("[DimRed] Running UMAP...")
  set.seed(config$seed)
  
  # Check if markers exist
  if(!all(markers %in% rownames(sce))) stop("Markers for UMAP missing in SCE.")
  
  sce <- scater::runUMAP(
    sce,
    subset_row = markers,
    exprs_values = "exprs",
    n_neighbors = config$n_neighbors,
    min_dist = config$min_dist,
    name = "UMAP"
  )
  
  return(sce)
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