#!/usr/bin/env Rscript
# src/steps/06_evaluate_integration.R
# ==============================================================================
# STEP 06: INTEGRATION VALIDATION — PEER-REVIEWER STATISTICS
# Description: Computes all validation metrics for the population-level
#              integration: JSD between matched cluster profiles, Spearman
#              correlation of frequencies across batches, PERMANOVA on the
#              frequency matrix, and I² heterogeneity summary.
# Dependencies: vegan (PERMANOVA)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(vegan)
  library(writexl)
  library(SingleCellExperiment)
})

source("src/functions/cluster_matching.R")
source("src/functions/integration.R")

message("\n=== PIPELINE STEP 6: INTEGRATION VALIDATION ===")

config      <- read_yaml("config/global_params.yml")
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration %||% "results/integration"
if (!dir.exists(int_out_dir)) dir.create(int_out_dir, recursive = TRUE)

required_files <- c(
  freq_file    = file.path(out_dir, "frequency_matrix.rds"),
  baseline_file = file.path(out_dir, "healthy_baseline.rds"),
  match_file   = file.path(out_dir, "cluster_match_table.rds"),
  centroids_file = file.path(out_dir, "cluster_centroids.rds")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("[Step 06] Missing files: ", paste(names(missing), collapse = ", "),
       ". Run Steps 03-05 first.")
}

tryCatch({

  freq_result        <- readRDS(required_files["freq_file"])
  integration_result <- readRDS(required_files["baseline_file"])
  match_table        <- readRDS(required_files["match_file"])
  centroids_list     <- readRDS(required_files["centroids_file"])

  fm       <- freq_result$freq_matrix
  pop_cols <- grep("^pop_", colnames(fm), value = TRUE)
  batches  <- names(centroids_list)
  markers  <- colnames(centroids_list[[batches[1L]]])

  # --------------------------------------------------------------------------
  # METRIC 1: JSD between matched cluster centroid profiles
  # --------------------------------------------------------------------------
  message("[Step 06] Computing Jensen-Shannon Divergence per matched pair...")
  validation <- validate_cluster_matches(match_table, centroids_list, markers)

  cluster_match_report <- validation$summary
  if (nrow(cluster_match_report) > 0) {
    n_pass <- sum(cluster_match_report$jsd_pass, na.rm = TRUE)
    message(sprintf("[Step 06] JSD < 0.1: %d / %d matched populations",
                    n_pass, nrow(cluster_match_report)))
  }

  # --------------------------------------------------------------------------
  # METRIC 2: Cross-batch Spearman correlation of population frequencies
  # --------------------------------------------------------------------------
  message("[Step 06] Computing cross-batch Spearman frequency correlation...")

  batch_a <- batches[1L]
  batch_b <- batches[2L]
  fm_a    <- fm[fm$batch == batch_a, ]
  fm_b    <- fm[fm$batch == batch_b, ]

  cor_rows <- lapply(pop_cols, function(pop) {
    vec_a <- fm_a[[pop]]
    vec_b <- fm_b[[pop]]
    n_a   <- length(vec_a)
    n_b   <- length(vec_b)

    if (n_a < 3 || n_b < 3) {
      return(data.frame(population = pop, rho = NA_real_, p_value = NA_real_,
                        n_batch_a = n_a, n_batch_b = n_b, stringsAsFactors = FALSE))
    }
    ct <- cor.test(vec_a, vec_b, method = "spearman", exact = FALSE)
    data.frame(population = pop, rho = ct$estimate, p_value = ct$p.value,
               n_batch_a = n_a, n_batch_b = n_b, stringsAsFactors = FALSE,
               row.names = NULL)
  })
  cor_report <- do.call(rbind, cor_rows)

  mean_rho <- mean(cor_report$rho, na.rm = TRUE)
  message(sprintf("[Step 06] Mean Spearman rho across populations: %.3f", mean_rho))
  if (mean_rho < 0.80) {
    warning("[Step 06] Mean Spearman rho < 0.80. Cross-batch frequency concordance is low.",
            call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # METRIC 3: PERMANOVA on frequency matrix (batch as grouping variable)
  # --------------------------------------------------------------------------
  message("[Step 06] Running PERMANOVA (batch ~ frequency matrix)...")

  freq_wide <- fm[, pop_cols, drop = FALSE]
  # Apply arcsine-sqrt transformation for approximate normality
  freq_wide_trans <- as.data.frame(
    lapply(freq_wide, asin_sqrt_transform)
  )
  batch_vec <- fm$batch

  n_perm <- 999L
  set.seed(as.integer(config$random_seed %||% 42L))

  permanova_result <- vegan::adonis2(
    freq_wide_trans ~ batch_vec,
    data         = data.frame(batch_vec = batch_vec),
    method       = "euclidean",
    permutations = n_perm
  )

  permanova_df <- as.data.frame(permanova_result)
  permanova_df$term <- rownames(permanova_df)
  rownames(permanova_df) <- NULL

  pval_batch <- permanova_result["batch_vec", "Pr(>F)"]
  if (!is.na(pval_batch)) {
    status <- if (pval_batch > 0.05) "PASS (batch not significant)" else
              "WARN (batch remains significant)"
    message(sprintf("[Step 06] PERMANOVA p-value (batch): %.4f — %s", pval_batch, status))
  }

  # --------------------------------------------------------------------------
  # METRIC 4: Heterogeneity I² summary
  # --------------------------------------------------------------------------
  message("[Step 06] Compiling I² heterogeneity summary...")

  baseline <- integration_result$baseline
  het_report <- data.frame(
    population         = baseline$population,
    I2_pct             = baseline$I2_pct,
    tau2               = baseline$tau2,
    Q                  = baseline$Q,
    k                  = baseline$k,
    flag_heterogeneity = baseline$flag_heterogeneity,
    interpretation     = ifelse(
      baseline$I2_pct < 25, "low",
      ifelse(baseline$I2_pct < 75, "moderate", "high")
    ),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------------
  # EXPORT
  # --------------------------------------------------------------------------
  writexl::write_xlsx(
    list(
      ClusterMatch         = cluster_match_report,
      FrequencyCorrelation = cor_report,
      PERMANOVA            = permanova_df,
      HeterogeneityI2      = het_report
    ),
    path = file.path(int_out_dir, "integration_validation_report.xlsx")
  )

  message(sprintf("[Step 06] Validation report saved to: %s",
                  file.path(int_out_dir, "integration_validation_report.xlsx")))

  message("\n=== STEP 6 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 06 failed:")
  message(e$message)
  stop(paste("[Step 06 Fatal]", e$message), call. = FALSE)
})
