# src/functions/correction.R
# ==============================================================================
# BATCH CORRECTION MODULE
# Description: Wrapper functions for cyCombine execution.
# Dependencies: cyCombine, SingleCellExperiment, magrittr, dplyr, tibble
# ==============================================================================

library(cyCombine)
library(SingleCellExperiment)
library(magrittr)
library(dplyr)
library(tibble)

#' @title Run cyCombine Batch Correction
#' @description Performs EMD-based batch correction using cyCombine.
#' @param sce SingleCellExperiment object containing arcsinh-transformed data.
#' @param markers Character vector of markers to use for correction.
#' @param batch_col String name of the column in colData containing batch IDs.
#' @param config List of correction parameters (xdim, ydim, rlen, norm_method).
#' @return A SingleCellExperiment object with corrected expression values.
run_batch_correction <- function(sce, markers, batch_col = "batch", config) {
  
  # 1. Input Validation
  if (!all(markers %in% rownames(sce))) {
    missing <- setdiff(markers, rownames(sce))
    stop(sprintf("[Correction Error] Markers missing in SCE: %s", paste(missing, collapse=", ")))
  }
  
  message("[Correction] Preparing data for cyCombine...")
  
  # Prepare Input Tibble
  # CRITICAL: Use check.names = FALSE to prevent R from renaming "CD4+" to "CD4."
  # This prevents mismatch errors downstream.
  input_df <- tryCatch({
    exprs_t <- t(assay(sce, "exprs"))
    meta    <- colData(sce) %>% as.data.frame()
    
    # Efficiently combine matrix and metadata without mangling column names
    dplyr::bind_cols(as.data.frame(exprs_t, check.names = FALSE), meta) %>%
      as_tibble()
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] Data preparation failed: %s", e$message))
  })
  
  # Pre-Flight Check: Verify markers exist in the tibble columns
  # This catches "CD4+" vs "CD4." issues before calling the library
  missing_cols <- setdiff(markers, colnames(input_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("[Correction Fatal] Column mismatch! Markers not found in input_df: %s", 
                 paste(missing_cols, collapse=", ")))
  }
  
  # 2. Execute Batch Correction
  # FIXED: 
  # - Passed 'input_df' as the FIRST POSITIONAL argument (no name) to avoid "unused argument" errors.
  # - Removed 'sample_col' (not required/supported by this wrapper).
  message(sprintf("[Correction] Running cyCombine::batch_correct (Grid: %dx%d, Rlen: %d)...", 
                  config$xdim, config$ydim, config$rlen))
  
  corrected_df <- tryCatch({
    cyCombine::batch_correct(
      input_df,                   # Positional argument (safest for API variations)
      markers = markers,
      covar = batch_col,
      norm_method = config$norm_method,
      xdim = config$xdim,
      ydim = config$ydim,
      rlen = config$rlen,
      seed = config$seed
    )
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] cyCombine execution failed: %s", e$message))
  })
  
  # 3. Reconstruct SingleCellExperiment
  message("[Correction] Rebuilding SingleCellExperiment object...")
  
  # Validation: Ensure row count integrity
  if (nrow(corrected_df) != ncol(sce)) {
    stop(sprintf("[Correction Fatal] Output rows (%d) do not match input cells (%d).", 
                 nrow(corrected_df), ncol(sce)))
  }
  
  # Extract corrected expression matrix
  # We select only the markers to avoid duplicating metadata
  corrected_exprs <- corrected_df %>%
    select(all_of(markers)) %>%
    as.matrix() %>%
    t()
  
  # Create new SCE retaining original metadata
  sce_corr <- SingleCellExperiment(
    assays = list(exprs = corrected_exprs),
    colData = colData(sce)
  )
  
  return(sce_corr)
}