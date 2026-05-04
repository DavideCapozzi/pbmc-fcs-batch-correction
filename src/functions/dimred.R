# src/functions/dimred.R
# ==============================================================================
# DIMENSIONALITY REDUCTION MODULE
# Description: Wrapper for UMAP generation and plotting (Pre/Post correction QC).
# Dependencies: scater, ggplot2, SingleCellExperiment
# ==============================================================================

library(scater)
library(ggplot2)
library(SingleCellExperiment)
library(patchwork)

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

#' @title Plot UMAP Colored by Matched Population
#' @description UMAP for a single batch colored by matched_population_id.
#'   Cells with NA (unmatched clusters) are shown in grey.
#' @param sce SCE with UMAP in reducedDims and matched_population_id in colData.
#' @param batch_id Character string for the plot title.
#' @param out_path Path to save PDF.
plot_umap_matched_populations <- function(sce, batch_id, out_path) {

  message(sprintf("[DimRed] Plotting matched-population UMAP for batch: %s", batch_id))

  p <- scater::plotReducedDim(
    sce, dimred = "UMAP",
    colour_by  = "matched_population_id",
    other_fields = "matched_population_id"
  ) +
    ggplot2::scale_colour_discrete(na.value = "grey80") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste("Matched Populations —", batch_id))

  ggplot2::ggsave(out_path, plot = p, width = 8, height = 6)
}

#' @title Plot Side-by-Side UMAP for Both Batches
#' @description Assembles per-batch UMAPs with a consistent population color
#'   palette so matched populations share the same color across panels.
#' @param sce_list Named list of per-batch SCE objects (each with UMAP computed).
#' @param out_path Path to save PDF.
plot_umap_batch_comparison <- function(sce_list, out_path) {

  message("[DimRed] Generating side-by-side batch comparison UMAP...")

  # Build a shared discrete color palette across all matched_population_id levels
  all_levels <- unique(unlist(lapply(sce_list, function(s) {
    as.character(s$matched_population_id)
  })))
  all_levels <- sort(all_levels[!is.na(all_levels)])

  n_levels  <- length(all_levels)
  palette   <- if (n_levels <= 20) scales::hue_pal()(n_levels) else
               colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                                  "#FF7F00","#A65628","#F781BF","#999999"))(n_levels)
  color_map <- setNames(palette, all_levels)

  plots <- lapply(names(sce_list), function(b) {
    sce_b <- sce_list[[b]]
    pop_factor <- factor(as.character(sce_b$matched_population_id),
                         levels = c(all_levels, NA))

    umap_coords <- as.data.frame(reducedDim(sce_b, "UMAP"))
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$pop <- pop_factor

    ggplot2::ggplot(umap_coords, ggplot2::aes(x = UMAP1, y = UMAP2, colour = pop)) +
      ggplot2::geom_point(size = 0.3, alpha = 0.5) +
      ggplot2::scale_colour_manual(values = color_map, na.value = "grey80",
                                   name = "Population", drop = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(b) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      ))
  })

  combined <- patchwork::wrap_plots(plots, ncol = length(sce_list)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = "Cross-Batch Matched Population Comparison")

  ggplot2::ggsave(out_path, plot = combined,
                  width = 8 * length(sce_list), height = 7)
}