# src/legacy/functions/correction.R
# ==============================================================================
# [LEGACY] BATCH CORRECTION MODULE — MFI ALIGNMENT PATH
# Description: cyCombine EMD-based batch correction wrapper. Used only when
#              clustering.run_per_batch = false. Superseded by the stratified
#              ecological baseline approach (run_per_batch = true).
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

  if (!all(markers %in% rownames(sce))) {
    missing <- setdiff(markers, rownames(sce))
    stop(sprintf("[Correction Error] Markers missing in SCE: %s", paste(missing, collapse = ", ")))
  }

  message("[Correction] Preparing data for cyCombine...")

  input_df <- tryCatch({
    exprs_t <- t(assay(sce, "exprs"))
    meta    <- colData(sce) %>% as.data.frame()
    dplyr::bind_cols(as.data.frame(exprs_t, check.names = FALSE), meta) %>% as_tibble()
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] Data preparation failed: %s", e$message))
  })

  missing_cols <- setdiff(markers, colnames(input_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("[Correction Fatal] Column mismatch — markers not found in input_df: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  message(sprintf("[Correction] Running cyCombine::batch_correct (Grid: %dx%d, Rlen: %d)...",
                  config$xdim, config$ydim, config$rlen))

  corrected_df <- tryCatch({
    cyCombine::batch_correct(
      input_df,
      markers     = markers,
      covar       = batch_col,
      norm_method = config$norm_method,
      xdim        = config$xdim,
      ydim        = config$ydim,
      rlen        = config$rlen,
      seed        = config$seed
    )
  }, error = function(e) {
    stop(sprintf("[Correction Fatal] cyCombine execution failed: %s", e$message))
  })

  message("[Correction] Rebuilding SingleCellExperiment object...")

  if (nrow(corrected_df) != ncol(sce)) {
    stop(sprintf("[Correction Fatal] Output rows (%d) do not match input cells (%d).",
                 nrow(corrected_df), ncol(sce)))
  }

  corrected_exprs <- corrected_df %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    as.matrix() %>%
    t()

  SingleCellExperiment(
    assays  = list(exprs = corrected_exprs),
    colData = colData(sce)
  )
}
