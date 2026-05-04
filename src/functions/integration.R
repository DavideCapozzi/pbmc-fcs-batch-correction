# src/functions/integration.R
# ==============================================================================
# RANDOM EFFECTS META-ANALYSIS MODULE
# Description: DerSimonian-Laird (DL) pooling of population frequencies across
#              batches. Uses arcsine-square-root transformation for variance
#              stabilization of proportions. The DL implementation was
#              cross-validated against metafor::rma().
# Dependencies: none (base R only)
# ==============================================================================

#' @title Arcsine Square Root Transformation
#' @description Variance-stabilizing transformation for proportions: y = arcsin(sqrt(p)).
#' @param freq_vec Numeric vector of proportions in [0,1].
#' @return Numeric vector of transformed values.
asin_sqrt_transform <- function(freq_vec) {
  asin(sqrt(pmax(0, pmin(1, freq_vec))))
}

#' @title Back-Transform from Arcsine Square Root Scale
#' @param y_vec Numeric vector on the arcsin(sqrt) scale.
#' @return Numeric vector of proportions in [0,1].
asin_sqrt_backtransform <- function(y_vec) {
  sin(y_vec)^2
}

#' @title DerSimonian-Laird Random Effects Meta-Analysis (Single Population)
#' @description Internal function. Pools one population's transformed frequencies.
#' @param y Numeric vector of arcsin-sqrt transformed frequencies (one per sample).
#' @param n Integer vector of total cell counts per sample (used to estimate v_s).
#' @return Named numeric vector: y_pool, se_pool, ci_lower, ci_upper,
#'   tau2, Q, I2, k.
.dl_single_population <- function(y, n) {
  k   <- length(y)
  v_s <- 1 / (4 * n)
  w_s <- 1 / v_s

  y_bar_fe <- sum(w_s * y) / sum(w_s)
  Q        <- sum(w_s * (y - y_bar_fe)^2)
  C        <- sum(w_s) - sum(w_s^2) / sum(w_s)
  tau2     <- max(0, (Q - (k - 1)) / C)

  w_star    <- 1 / (v_s + tau2)
  y_pool    <- sum(w_star * y) / sum(w_star)
  se_pool   <- sqrt(1 / sum(w_star))
  ci_lower  <- y_pool - 1.96 * se_pool
  ci_upper  <- y_pool + 1.96 * se_pool
  I2        <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0

  c(y_pool   = y_pool,
    se_pool  = se_pool,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    tau2     = tau2,
    Q        = Q,
    I2       = I2,
    k        = k)
}

#' @title DerSimonian-Laird Random Effects Integration
#' @description Pools population frequencies across all samples using the DL
#'   estimator after arcsine-square-root transformation. Produces back-transformed
#'   pooled frequencies, 95% analytical CIs, tau2, and I2. Optionally runs a
#'   stratified bootstrap to validate the analytical CIs.
#' @param freq_matrix data.frame with columns: sample_id, batch, total_cells,
#'   pop_1, pop_2, ... (proportions in [0,1]).
#' @param count_matrix data.frame same structure with raw integer counts.
#' @param n_bootstrap Integer. Bootstrap iterations. Default 1000.
#' @param seed Integer. Random seed. Default 42.
#' @return List:
#'   $baseline    data.frame [n_populations x cols]
#'   $per_sample  data.frame with DL weights per sample per population
#'   $bootstrap   data.frame with bootstrap CI bounds per population
run_random_effects_integration <- function(freq_matrix, count_matrix,
                                           n_bootstrap = 1000L, seed = 42L) {

  pop_cols <- grep("^pop_", colnames(freq_matrix), value = TRUE)
  batches  <- sort(unique(freq_matrix$batch))

  # -- Analytical DL ----------------------------------------------------------
  baseline_list <- lapply(pop_cols, function(pop) {
    y <- asin_sqrt_transform(freq_matrix[[pop]])
    n <- as.integer(freq_matrix$total_cells)
    r <- .dl_single_population(y, n)

    data.frame(
      population        = pop,
      freq_pooled       = asin_sqrt_backtransform(r["y_pool"]),
      freq_lower_95ci   = asin_sqrt_backtransform(r["ci_lower"]),
      freq_upper_95ci   = asin_sqrt_backtransform(r["ci_upper"]),
      tau2              = r["tau2"],
      I2_pct            = r["I2"],
      Q                 = r["Q"],
      k                 = as.integer(r["k"]),
      flag_heterogeneity = r["I2"] > 75,
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
  })

  baseline <- do.call(rbind, baseline_list)

  # -- Per-sample DL weights (for reporting) ----------------------------------
  per_sample_list <- lapply(pop_cols, function(pop) {
    y    <- asin_sqrt_transform(freq_matrix[[pop]])
    n    <- as.integer(freq_matrix$total_cells)
    v_s  <- 1 / (4 * n)
    tau2 <- baseline$tau2[baseline$population == pop]
    w_re <- 1 / (v_s + tau2)
    data.frame(
      sample_id  = freq_matrix$sample_id,
      batch      = freq_matrix$batch,
      population = pop,
      y_trans    = y,
      weight_re  = w_re / sum(w_re),
      stringsAsFactors = FALSE
    )
  })
  per_sample <- do.call(rbind, per_sample_list)

  # -- Stratified Bootstrap ---------------------------------------------------
  set.seed(seed)
  boot_pool <- matrix(NA_real_, nrow = n_bootstrap, ncol = length(pop_cols),
                      dimnames = list(NULL, pop_cols))

  for (b in seq_len(n_bootstrap)) {
    boot_idx <- unlist(lapply(batches, function(bat) {
      idx <- which(freq_matrix$batch == bat)
      sample(idx, length(idx), replace = TRUE)
    }))
    fm_b <- freq_matrix[boot_idx, ]
    cm_b <- count_matrix[boot_idx, ]
    for (pop in pop_cols) {
      y_b   <- asin_sqrt_transform(fm_b[[pop]])
      n_b   <- as.integer(cm_b$total_cells)
      r_b   <- .dl_single_population(y_b, n_b)
      boot_pool[b, pop] <- asin_sqrt_backtransform(r_b["y_pool"])
    }
  }

  bootstrap <- do.call(rbind, lapply(pop_cols, function(pop) {
    data.frame(
      population    = pop,
      boot_ci_lower = quantile(boot_pool[, pop], 0.025, na.rm = TRUE),
      boot_ci_upper = quantile(boot_pool[, pop], 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  list(baseline = baseline, per_sample = per_sample, bootstrap = bootstrap)
}

#' @title Format Baseline Report
#' @description Merges DL results with bootstrap CIs into a single summary table.
#' @param integration_result List from run_random_effects_integration().
#' @param match_table data.frame from match_metaclusters().
#' @return data.frame ready for Excel export.
format_baseline_report <- function(integration_result, match_table) {
  baseline  <- integration_result$baseline
  bootstrap <- integration_result$bootstrap

  report <- merge(baseline, bootstrap, by = "population", all.x = TRUE)

  pop_id_map <- match_table[!is.na(match_table$matched_population_id),
                             c("matched_population_id", "cluster_a", "cluster_b")]
  pop_id_map <- unique(pop_id_map)
  pop_id_map$population <- paste0("pop_", pop_id_map$matched_population_id)

  report <- merge(report, pop_id_map[, c("population", "cluster_a", "cluster_b")],
                  by = "population", all.x = TRUE)

  col_order <- c("population", "cluster_a", "cluster_b",
                 "freq_pooled", "freq_lower_95ci", "freq_upper_95ci",
                 "boot_ci_lower", "boot_ci_upper",
                 "tau2", "I2_pct", "Q", "k", "flag_heterogeneity")
  col_order <- col_order[col_order %in% colnames(report)]
  report[, col_order]
}
