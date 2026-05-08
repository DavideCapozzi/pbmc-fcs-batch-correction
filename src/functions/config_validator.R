# src/functions/config_validator.R
# ==============================================================================
# CONFIG VALIDATION
# Validates the global_params.yml config and normalises types that YAML may
# parse inconsistently (e.g., integer scalars read as lists in some R builds).
# Call validate_config(config) immediately after read_yaml() in main.R.
# ==============================================================================

#' Validate and normalise the pipeline configuration.
#' @param config list produced by yaml::read_yaml()
#' @return The validated (and type-corrected) config list, invisibly.
#' @details Stops with an informative error on the first violation found.
validate_config <- function(config) {

  .require_key <- function(obj, keys, context = "config") {
    for (k in keys) {
      if (is.null(obj[[k]])) {
        stop(sprintf("[Config] Missing required key: %s$%s", context, k), call. = FALSE)
      }
    }
  }

  .check_range <- function(val, lo, hi, name) {
    v <- suppressWarnings(as.numeric(val))
    if (is.na(v) || v < lo || v > hi) {
      stop(sprintf("[Config] %s must be in [%g, %g], got: %s", name, lo, hi, val),
           call. = FALSE)
    }
  }

  .check_positive <- function(val, name) {
    v <- suppressWarnings(as.numeric(val))
    if (is.na(v) || v <= 0) {
      stop(sprintf("[Config] %s must be > 0, got: %s", name, val), call. = FALSE)
    }
  }

  # ── Top-level required keys ─────────────────────────────────────────────────
  .require_key(config, c("directories", "markers", "qc", "clustering", "integration",
                          "auto_phenotyping", "scoring", "dimred"))

  # ── directories ─────────────────────────────────────────────────────────────
  .require_key(config$directories, c("raw", "intermediate", "processed", "figures",
                                      "logs", "integration"), context = "directories")

  # ── markers ─────────────────────────────────────────────────────────────────
  .check_positive(config$markers$transform_cofactor, "markers$transform_cofactor")

  # ── qc ──────────────────────────────────────────────────────────────────────
  if (isTRUE(config$qc$enabled)) {
    pqc <- config$qc$peacoqc
    if (!is.null(pqc)) {
      .check_range(pqc$IT_limit %||% 0.6, 0, 1, "qc$peacoqc$IT_limit")
      .check_positive(pqc$MAD %||% 6, "qc$peacoqc$MAD")
    }
    if (!is.null(config$qc$max_pct_removed)) {
      .check_range(config$qc$max_pct_removed, 0, 100, "qc$max_pct_removed")
    }
  }

  # ── clustering ──────────────────────────────────────────────────────────────
  clus <- config$clustering
  .check_positive(clus$xdim %||% 10, "clustering$xdim")
  .check_positive(clus$ydim %||% 10, "clustering$ydim")
  .check_positive(clus$n_metaclusters %||% 20, "clustering$n_metaclusters")

  # ── integration ─────────────────────────────────────────────────────────────
  .check_range(config$integration$similarity_threshold %||% 0.7, 0, 1,
               "integration$similarity_threshold")

  # ── scoring ─────────────────────────────────────────────────────────────────
  .check_positive(config$scoring$z_score_threshold %||% 2, "scoring$z_score_threshold")

  # ── n_cores ──────────────────────────────────────────────────────────────────
  n_cores_raw <- config$n_cores %||% 1L
  config$n_cores <- max(1L, as.integer(n_cores_raw))

  # ── random_seed ──────────────────────────────────────────────────────────────
  config$random_seed <- as.integer(config$random_seed %||% 42L)

  # ── YAML integer normalisation for FlowSOM params ───────────────────────────
  config$clustering$xdim          <- as.integer(config$clustering$xdim %||% 10L)
  config$clustering$ydim          <- as.integer(config$clustering$ydim %||% 10L)
  config$clustering$n_metaclusters <- as.integer(config$clustering$n_metaclusters %||% 20L)

  message("[Config] Validation passed.")
  invisible(config)
}
