# src/functions/qc.R
# ==============================================================================
# PRE-PROCESSING QC MODULE
# Programmatic cell filtering before FlowSOM clustering.
# Pipeline: PeacoQC (temporal drift + margin events) →
#           singlet gate (FSC-H/FSC-A slope) →
#           debris gate (FSC-A vs SSC-A window) →
#           viability gate (fluorescent dye exclusion, e.g. ViaKrome)
#
# All index arithmetic uses 1-based integer vectors remapped through the global
# frame at each step, so final indices are always relative to the raw FCS file.
# ==============================================================================

library(flowCore)

# ------------------------------------------------------------------------------
# 1. RAW FCS READER (scatter channels retained, no transformation)
# ------------------------------------------------------------------------------

#' Read an FCS file keeping only scatter channels for QC gating.
#' Uses transformation = FALSE and truncate_max_range = FALSE so that PeacoQC
#' receives unclipped values for margin-event detection.
#' @param file_path character(1) absolute path to FCS file
#' @param scatter_channels character vector of $PnN channel names to retain
#' @return list($flowframe, $scatter_mat [n_cells x n_channels], $n_raw)
read_fcs_with_scatter <- function(
    file_path,
    scatter_channels = c("FSC-A", "FSC-H", "SSC-A", "FSC-Width", "Time")
) {
  ff <- flowCore::read.FCS(
    file_path,
    transformation      = FALSE,
    truncate_max_range  = FALSE,
    emptyValue          = FALSE
  )

  all_cols      <- colnames(flowCore::exprs(ff))
  keep_cols     <- intersect(scatter_channels, all_cols)
  scatter_mat   <- flowCore::exprs(ff)[, keep_cols, drop = FALSE]

  list(
    flowframe   = ff,
    scatter_mat = scatter_mat,
    n_raw       = nrow(scatter_mat)
  )
}

# ------------------------------------------------------------------------------
# 2. PEACOQC — temporal drift and margin event removal
# ------------------------------------------------------------------------------

#' Apply PeacoQC to remove time-dependent anomalies and detector-saturated events.
#' @param flowframe flowFrame from read_fcs_with_scatter()$flowframe
#' @param channels_to_use character vector of fluorescent channel names, or NULL
#'   (auto-detected: all channels that are not FSC/SSC/Time/Width)
#' @param peacoqc_params named list with PeacoQC tuning parameters
#' @return list($valid_indices, $n_before, $n_after, $pct_removed, $peacoqc_output)
apply_peacoqc_filter <- function(
    flowframe,
    channels_to_use  = NULL,
    peacoqc_params   = list(
      MAD          = 6,
      IT_limit     = 0.6,
      remove_zeros = FALSE,
      min_cells    = 150
    )
) {
  if (!requireNamespace("PeacoQC", quietly = TRUE)) {
    stop(
      "[QC] PeacoQC is not installed. Install via:\n",
      "  BiocManager::install('PeacoQC')\n",
      "or add 'bioconductor-peacoqc' to env/environment.yml and update the conda env."
    )
  }

  all_cols <- colnames(flowCore::exprs(flowframe))

  if (is.null(channels_to_use)) {
    scatter_regex    <- "^FSC|^SSC|^Time$|Width"
    channels_to_use  <- all_cols[!grepl(scatter_regex, all_cols, ignore.case = TRUE)]
  }

  channels_to_use <- intersect(channels_to_use, all_cols)
  if (length(channels_to_use) == 0L) {
    stop("[QC] apply_peacoqc_filter: no valid fluorescent channels found in flowFrame.")
  }

  n_before <- nrow(flowCore::exprs(flowframe))

  pqc_out <- PeacoQC::PeacoQC(
    ff           = flowframe,
    channels     = channels_to_use,
    MAD          = as.numeric(peacoqc_params$MAD %||% 6),
    IT_limit     = as.numeric(peacoqc_params$IT_limit %||% 0.6),
    remove_zeros = isTRUE(peacoqc_params$remove_zeros),
    min_cells    = as.integer(peacoqc_params$min_cells %||% 150L),
    plot         = FALSE,
    save_fcs     = FALSE
  )

  good_cells    <- pqc_out$GoodCells
  valid_indices <- which(good_cells)
  n_after       <- length(valid_indices)

  list(
    valid_indices   = valid_indices,
    n_before        = n_before,
    n_after         = n_after,
    pct_removed     = 100 * (1 - n_after / n_before),
    peacoqc_output  = pqc_out
  )
}

# ------------------------------------------------------------------------------
# 3. SINGLET GATE — FSC-H / FSC-A slope method
# ------------------------------------------------------------------------------

#' Gate singlets using the FSC-H vs FSC-A ratio (slope method).
#' Doublets have a lower FSC-H/FSC-A ratio than singlets; the gate is a
#' symmetric window around the population median measured in robust MAD units.
#' @param scatter_mat numeric matrix with columns FSC-A and FSC-H
#' @param fsc_a_col character(1) column name for FSC area
#' @param fsc_h_col character(1) column name for FSC height
#' @param slope_sd_multiplier numeric(1) gate half-width in MAD units (default 3)
#' @return list($valid_indices, $n_before, $n_after, $pct_removed,
#'              $gate_center, $gate_mad, $gate_lower, $gate_upper)
apply_singlet_gate <- function(
    scatter_mat,
    fsc_a_col          = "FSC-A",
    fsc_h_col          = "FSC-H",
    slope_sd_multiplier = 3
) {
  stopifnot(
    fsc_a_col %in% colnames(scatter_mat),
    fsc_h_col %in% colnames(scatter_mat)
  )

  fsc_a  <- scatter_mat[, fsc_a_col]
  fsc_h  <- scatter_mat[, fsc_h_col]
  n_before <- nrow(scatter_mat)

  ratio <- fsc_h / fsc_a

  # Physical guard: ratios outside (0.1, 2.0) are instrument artefacts
  physical <- ratio > 0.1 & ratio < 2.0 & is.finite(ratio)
  center   <- median(ratio[physical], na.rm = TRUE)
  mad_val  <- mad(ratio[physical],    na.rm = TRUE)

  lower <- center - slope_sd_multiplier * mad_val
  upper <- center + slope_sd_multiplier * mad_val

  valid_indices <- which(ratio >= lower & ratio <= upper & is.finite(ratio))
  n_after       <- length(valid_indices)

  list(
    valid_indices = valid_indices,
    n_before      = n_before,
    n_after       = n_after,
    pct_removed   = 100 * (1 - n_after / n_before),
    gate_center   = center,
    gate_mad      = mad_val,
    gate_lower    = lower,
    gate_upper    = upper
  )
}

# ------------------------------------------------------------------------------
# 4. DEBRIS GATE — FSC-A vs SSC-A adaptive window
# ------------------------------------------------------------------------------

#' Remove small debris (low FSC-A) and large aggregates (high SSC-A).
#' Thresholds are estimated per-sample from quantiles, making the gate adaptive
#' to sample-to-sample variation in cell density and viability.
#' @param scatter_mat numeric matrix with columns FSC-A and SSC-A
#' @param fsc_a_col character(1)
#' @param ssc_a_col character(1)
#' @param fsc_min_quantile numeric(1) lower FSC-A quantile to discard (default 0.02)
#' @param ssc_max_quantile numeric(1) upper SSC-A quantile to discard (default 0.98)
#' @return list($valid_indices, $n_before, $n_after, $pct_removed,
#'              $fsc_min_thresh, $ssc_max_thresh)
apply_debris_gate <- function(
    scatter_mat,
    fsc_a_col        = "FSC-A",
    ssc_a_col        = "SSC-A",
    fsc_min_quantile  = 0.02,
    ssc_max_quantile  = 0.98
) {
  stopifnot(
    fsc_a_col %in% colnames(scatter_mat),
    ssc_a_col %in% colnames(scatter_mat)
  )

  fsc_a    <- scatter_mat[, fsc_a_col]
  ssc_a    <- scatter_mat[, ssc_a_col]
  n_before <- nrow(scatter_mat)

  fsc_min <- quantile(fsc_a, probs = fsc_min_quantile, na.rm = TRUE)
  ssc_max <- quantile(ssc_a, probs = ssc_max_quantile, na.rm = TRUE)

  valid_indices <- which(fsc_a >= fsc_min & ssc_a <= ssc_max & is.finite(fsc_a) & is.finite(ssc_a))
  n_after       <- length(valid_indices)

  list(
    valid_indices  = valid_indices,
    n_before       = n_before,
    n_after        = n_after,
    pct_removed    = 100 * (1 - n_after / n_before),
    fsc_min_thresh = as.numeric(fsc_min),
    ssc_max_thresh = as.numeric(ssc_max)
  )
}

# ------------------------------------------------------------------------------
# 5. VIABILITY GATE — dead cell exclusion via fluorescent viability dye
# ------------------------------------------------------------------------------

#' Exclude dead cells using a viability dye channel (e.g. ViaKrome, DAPI, 7-AAD).
#' Dead cells take up the dye and show high fluorescence; the gate keeps only cells
#' below median + n*MAD of the linear-scale signal.
#' The channel is located by case-insensitive pattern match against both $PnN
#' (channel name) and $PnS (description) in the flowFrame header, so the caller
#' can pass the semantic name ("VIAKROME") without knowing the raw detector name.
#' @param flowframe flowFrame — subset already filtered by upstream gates
#' @param channel_name character(1) semantic or partial channel name to search for
#' @param max_sd_multiplier numeric(1) threshold = median + n*MAD (default 3)
#' @return list($valid_indices, $n_before, $n_after, $pct_removed,
#'              $channel_found, $channel_raw, $threshold)
apply_viability_gate <- function(
    flowframe,
    channel_name,
    max_sd_multiplier = 3
) {
  all_pnn <- colnames(flowCore::exprs(flowframe))
  all_pns <- as.character(flowCore::pData(flowCore::parameters(flowframe))$desc)

  match_idx <- unique(c(
    grep(channel_name, all_pnn, ignore.case = TRUE),
    grep(channel_name, all_pns, ignore.case = TRUE)
  ))

  n_before <- nrow(flowCore::exprs(flowframe))

  if (length(match_idx) == 0L) {
    warning(sprintf(
      "[QC] Viability channel '%s' not found in FCS — dead cell gate skipped.", channel_name
    ), call. = FALSE)
    return(list(
      valid_indices = seq_len(n_before), n_before = n_before, n_after = n_before,
      pct_removed = 0, channel_found = FALSE, channel_raw = NA_character_, threshold = NA_real_
    ))
  }

  raw_channel <- all_pnn[match_idx[1L]]
  values      <- flowCore::exprs(flowframe)[, raw_channel]

  med       <- median(values, na.rm = TRUE)
  mad_val   <- mad(values,    na.rm = TRUE)
  threshold <- med + max_sd_multiplier * mad_val

  valid_indices <- which(values <= threshold & is.finite(values))
  n_after       <- length(valid_indices)

  list(
    valid_indices = valid_indices,
    n_before      = n_before,
    n_after       = n_after,
    pct_removed   = 100 * (1 - n_after / n_before),
    channel_found = TRUE,
    channel_raw   = raw_channel,
    threshold     = threshold
  )
}

# ------------------------------------------------------------------------------
# 6. ORCHESTRATOR — per-sample QC pipeline
# ------------------------------------------------------------------------------

#' Run the full QC pipeline for a single FCS file.
#' Applies PeacoQC → singlet gate → debris gate → viability gate in sequence.
#' Indices are remapped through the global cell frame at each step so that
#' the returned valid_indices are always relative to the original FCS row order.
#'
#' @param file_path character(1) absolute path to FCS file
#' @param qc_config list — the qc: section of global_params.yml
#' @param sample_id character(1) used only for log messages; defaults to basename
#' @return list($valid_indices, $qc_metrics, $file_path, $sample_id)
run_sample_qc <- function(file_path, qc_config, sample_id = NULL) {
  if (is.null(sample_id)) sample_id <- basename(file_path)

  raw <- read_fcs_with_scatter(file_path)
  global_valid <- seq_len(raw$n_raw)

  qc_metrics <- list(
    n_raw              = raw$n_raw,
    n_after_peacoqc    = raw$n_raw,
    n_after_singlet    = raw$n_raw,
    n_after_debris     = raw$n_raw,
    n_after_viability  = raw$n_raw,
    n_final            = raw$n_raw
  )

  # Step A: PeacoQC (with adaptive two-pass retry)
  # If PeacoQC removes > max_pct_per_step of events on first pass, retry with a
  # relaxed IT_limit. Both passes are logged for audit trail.
  peacoqc_cfg   <- qc_config$peacoqc %||% list()
  if (isTRUE(peacoqc_cfg$enabled %||% TRUE)) {
    max_pct_step  <- as.numeric(peacoqc_cfg$max_pct_per_step %||% 0.60)
    retry_IT_step <- as.numeric(peacoqc_cfg$retry_IT_limit_step %||% 0.1)

    pqc <- apply_peacoqc_filter(raw$flowframe, channels_to_use = NULL, peacoqc_params = peacoqc_cfg)
    global_valid                        <- global_valid[pqc$valid_indices]
    qc_metrics$pct_peacoqc_removed      <- pqc$pct_removed
    qc_metrics$peacoqc_IT_limit_used    <- as.numeric(peacoqc_cfg$IT_limit %||% 0.6)
    qc_metrics$peacoqc_retry_triggered  <- FALSE

    if (pqc$pct_removed > max_pct_step * 100) {
      global_valid_retry <- seq_len(raw$n_raw)
      new_IT    <- as.numeric(peacoqc_cfg$IT_limit %||% 0.6) + retry_IT_step
      retry_cfg <- modifyList(peacoqc_cfg, list(IT_limit = new_IT))
      pqc2      <- apply_peacoqc_filter(raw$flowframe, channels_to_use = NULL,
                                        peacoqc_params = retry_cfg)
      if (pqc2$pct_removed < pqc$pct_removed) {
        qc_metrics$pct_peacoqc_removed_first_pass <- pqc$pct_removed
        global_valid                       <- global_valid_retry[pqc2$valid_indices]
        qc_metrics$pct_peacoqc_removed     <- pqc2$pct_removed
        qc_metrics$peacoqc_IT_limit_used   <- new_IT
        qc_metrics$peacoqc_retry_triggered <- TRUE
        warning(sprintf(
          "[QC] Sample '%s': PeacoQC retry triggered (IT_limit %.2f→%.2f): %.1f%% → %.1f%% removed",
          sample_id, as.numeric(peacoqc_cfg$IT_limit %||% 0.6), new_IT,
          pqc$pct_removed, pqc2$pct_removed
        ), call. = FALSE)
      }
    }

    qc_metrics$n_after_peacoqc <- length(global_valid)
  } else {
    qc_metrics$n_after_peacoqc          <- length(global_valid)
    qc_metrics$pct_peacoqc_removed      <- 0
    qc_metrics$peacoqc_IT_limit_used    <- as.numeric(peacoqc_cfg$IT_limit %||% 0.6)
    qc_metrics$peacoqc_retry_triggered  <- FALSE
  }

  # Step B: Singlet gate
  singlet_cfg <- qc_config$singlet %||% list()
  if (isTRUE(singlet_cfg$enabled %||% TRUE)) {
    multiplier <- as.numeric(singlet_cfg$slope_sd_multiplier %||% 3)
    scatter_sub <- raw$scatter_mat[global_valid, , drop = FALSE]
    sg <- apply_singlet_gate(scatter_sub, slope_sd_multiplier = multiplier)
    global_valid               <- global_valid[sg$valid_indices]
    qc_metrics$n_after_singlet <- length(global_valid)
    qc_metrics$pct_singlet_removed <- sg$pct_removed
    qc_metrics$gate_params <- list(
      singlet_center = sg$gate_center,
      singlet_mad    = sg$gate_mad,
      singlet_lower  = sg$gate_lower,
      singlet_upper  = sg$gate_upper
    )
  } else {
    qc_metrics$n_after_singlet     <- length(global_valid)
    qc_metrics$pct_singlet_removed <- 0
  }

  # Step C: Debris gate
  debris_cfg <- qc_config$debris %||% list()
  if (isTRUE(debris_cfg$enabled %||% TRUE)) {
    fsc_q <- as.numeric(debris_cfg$fsc_min_quantile %||% 0.02)
    ssc_q <- as.numeric(debris_cfg$ssc_max_quantile %||% 0.98)
    scatter_sub2 <- raw$scatter_mat[global_valid, , drop = FALSE]
    dg <- apply_debris_gate(scatter_sub2, fsc_min_quantile = fsc_q, ssc_max_quantile = ssc_q)
    global_valid              <- global_valid[dg$valid_indices]
    qc_metrics$n_after_debris <- length(global_valid)
    qc_metrics$pct_debris_removed <- dg$pct_removed
    if (!is.null(qc_metrics$gate_params)) {
      qc_metrics$gate_params$fsc_min_thresh <- dg$fsc_min_thresh
      qc_metrics$gate_params$ssc_max_thresh <- dg$ssc_max_thresh
    } else {
      qc_metrics$gate_params <- list(
        fsc_min_thresh = dg$fsc_min_thresh,
        ssc_max_thresh = dg$ssc_max_thresh
      )
    }
  } else {
    qc_metrics$n_after_debris     <- length(global_valid)
    qc_metrics$pct_debris_removed <- 0
  }

  # Step D: Viability gate
  viability_cfg <- qc_config$viability %||% list()
  if (isTRUE(viability_cfg$enabled %||% FALSE)) {
    via_channel <- as.character(viability_cfg$channel %||% "VIAKROME")
    via_mult    <- as.numeric(viability_cfg$max_sd_multiplier %||% 3)
    ff_sub      <- raw$flowframe[global_valid, ]
    vg          <- apply_viability_gate(ff_sub, channel_name = via_channel,
                                        max_sd_multiplier = via_mult)
    global_valid                     <- global_valid[vg$valid_indices]
    qc_metrics$n_after_viability     <- length(global_valid)
    qc_metrics$pct_viability_removed <- vg$pct_removed
    qc_metrics$viability_channel_raw <- vg$channel_raw
    qc_metrics$viability_threshold   <- vg$threshold
  } else {
    qc_metrics$n_after_viability     <- length(global_valid)
    qc_metrics$pct_viability_removed <- 0
  }

  qc_metrics$n_final           <- length(global_valid)
  qc_metrics$pct_total_removed <- 100 * (1 - length(global_valid) / raw$n_raw)

  # Guard: catastrophic removal
  max_pct <- as.numeric(qc_config$max_pct_removed %||% 85)
  if (qc_metrics$pct_total_removed > max_pct) {
    stop(sprintf(
      "[QC] Sample '%s': %.1f%% of cells removed (threshold: %.0f%%). ",
      sample_id, qc_metrics$pct_total_removed, max_pct
    ), "Check data quality or raise qc.max_pct_removed in config.")
  }

  min_cells <- as.integer(qc_config$min_cells_final %||% 5000L)
  if (qc_metrics$n_final < min_cells) {
    warning(sprintf(
      "[QC] Sample '%s': only %d cells after QC (minimum recommended: %d).",
      sample_id, qc_metrics$n_final, min_cells
    ), call. = FALSE)
  }

  list(
    valid_indices = global_valid,
    qc_metrics    = qc_metrics,
    file_path     = file_path,
    sample_id     = sample_id
  )
}

# Null-coalescing operator (in case dplyr is not loaded in this module's scope)
`%||%` <- function(x, y) if (!is.null(x)) x else y
