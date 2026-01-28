# src/functions/correction.R
# ==============================================================================
# BATCH CORRECTION MODULE
# Description: Wrapper functions for cyCombine execution.
# Dependencies: cyCombine, SingleCellExperiment, magrittr, dplyr
# ==============================================================================

library(cyCombine)
library(SingleCellExperiment)
library(magrittr)
library(dplyr)

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
  
  # cyCombine expects a tibble with columns: markers + batch_col + sample_col
  # We extract data efficiently
  exprs_t <- t(assay(sce, "exprs"))
  meta    <- colData(sce) %>% as.data.frame()
  
  # Bind expression and metadata
  # Optimization: Use indices to avoid full data duplication if possible, 
  # but cyCombine needs a specific structure.
  input_df <- cbind(as.data.frame(exprs_t), meta) %>%
    as_tibble()
  
  # 2. Normalization & SOM Training (Learn the distribution shapes)
  message(sprintf("[Correction] Training SOM grid (%dx%d) with rlen=%d...", 
                  config$xdim, config$ydim, config$rlen))
  
  # Determine cofactor for denormalization if needed (usually 5 for cytof/flow)
  # Note: Since input is already asinh transformed, cyCombine uses it directly 
  # for determining anchors.
  
  labels <- tryCatch({
    cyCombine::prepare_data(
      data = input_df,
      markers = markers,
      covar = batch_col,
      sample_col = "sample_id", 
      norm_method = config$norm_method,
      xdim = config$xdim,
      ydim = config$ydim,
      rlen = config$rlen,
      seed = config$seed # Ensure reproducibility
    )
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] SOM Training failed: %s", e$message))
  })
  
  # 3. Apply Correction (EMD Matching)
  message("[Correction] Applying EMD correction...")
  
  corrected_df <- tryCatch({
    cyCombine::correct_data(
      data = input_df,
      label = labels,
      markers = markers,
      covar = batch_col,
      sample_col = "sample_id",
      prop = 0.1 # Default subsampling for EMD computation
    )
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] EMD Correction failed: %s", e$message))
  })
  
  # 4. Reconstruct SingleCellExperiment
  # The output of correct_data is a tibble. We need to put it back into SCE.
  message("[Correction] Rebuilding SingleCellExperiment object...")
  
  # Ensure order matches (cyCombine usually maintains order, but we verify dimensions)
  if (nrow(corrected_df) != ncol(sce)) {
    warning("[Correction Warning] Row count mismatch after correction. Checking alignment.")
  }
  
  # Extract corrected expression
  corrected_exprs <- corrected_df %>%
    select(all_of(markers)) %>%
    as.matrix() %>%
    t()
  
  # Create new SCE
  # We keep the original metadata
  sce_corr <- SingleCellExperiment(
    assays = list(exprs = corrected_exprs),
    colData = colData(sce)
  )
  
  return(sce_corr)
}