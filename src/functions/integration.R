# src/functions/integration.R
# ==============================================================================
# STRATIFIED ECOLOGICAL BASELINE MODULE
# Description: Computes independent reference distributions for each matched
#              population, stratified by batch. Batches are treated as isolated,
#              non-exchangeable experimental conditions — no cross-batch pooling
#              or normalization is performed. Future patient samples are evaluated
#              by Z-score against the specific baseline of their assay cohort.
# Dependencies: tidyr, dplyr
# ==============================================================================

library(tidyr)
library(dplyr)

#' @title Compute Stratified Ecological Baseline
#' @description For each (batch, population) stratum, computes descriptive
#'   statistics: mean, SD, and non-parametric 95% reference interval (2.5th
#'   and 97.5th percentiles). No transformation is applied — within-assay
#'   frequency distributions do not require variance stabilization.
#' @param freq_matrix data.frame with columns: sample_id, batch, total_cells,
#'   pop_1, pop_2, ... (proportions in [0,1]). Direct output of
#'   compute_population_frequencies()$freq_matrix.
#' @return data.frame (the "Reference Dictionary") with one row per
#'   (batch, population) and columns: batch, population, n_samples,
#'   mean_freq, sd_freq, ci_lower_95, ci_upper_95.
compute_stratified_baseline <- function(freq_matrix) {
  stopifnot(
    is.data.frame(freq_matrix),
    "batch" %in% colnames(freq_matrix),
    any(grepl("^pop_", colnames(freq_matrix)))
  )

  freq_matrix %>%
    dplyr::select(batch, tidyselect::starts_with("pop_")) %>%
    tidyr::pivot_longer(
      cols      = tidyselect::starts_with("pop_"),
      names_to  = "population",
      values_to = "frequency"
    ) %>%
    dplyr::group_by(batch, population) %>%
    dplyr::summarise(
      n_samples   = dplyr::n(),
      mean_freq   = mean(frequency,                    na.rm = TRUE),
      sd_freq     = sd(frequency,                      na.rm = TRUE),
      ci_lower_95 = quantile(frequency, probs = 0.025, na.rm = TRUE),
      ci_upper_95 = quantile(frequency, probs = 0.975, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    as.data.frame()
}

#' @title Format Stratified Baseline Report
#' @description Enriches the baseline dictionary with human-readable cluster
#'   labels from the match table, and optionally with biological population
#'   names from config (population_labels).
#' @param baseline_dict data.frame from compute_stratified_baseline().
#' @param match_table data.frame from match_metaclusters().
#' @param population_labels Named list mapping pop_N → biological name, or NULL.
#' @return data.frame ready for Excel export.
format_stratified_report <- function(baseline_dict, match_table,
                                     population_labels = NULL) {
  pop_id_map <- match_table %>%
    dplyr::filter(is_matched == TRUE, !is.na(matched_population_id)) %>%
    dplyr::select(matched_population_id, cluster_a, cluster_b) %>%
    dplyr::distinct() %>%
    dplyr::mutate(population = paste0("pop_", matched_population_id)) %>%
    dplyr::select(population, cluster_a, cluster_b)

  report <- baseline_dict %>%
    dplyr::left_join(pop_id_map, by = "population") %>%
    dplyr::select(batch, population, cluster_a, cluster_b,
                  n_samples, mean_freq, sd_freq, ci_lower_95, ci_upper_95)

  if (!is.null(population_labels) && length(population_labels) > 0) {
    label_df <- data.frame(
      population      = names(population_labels),
      population_name = unlist(population_labels, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    report <- dplyr::left_join(report, label_df, by = "population")
    col_order <- c("batch", "population", "population_name", "cluster_a", "cluster_b",
                   "n_samples", "mean_freq", "sd_freq", "ci_lower_95", "ci_upper_95")
    report <- dplyr::select(report, dplyr::any_of(col_order))
  }

  as.data.frame(report)
}

#' @title Format Centroid Report
#' @description Builds a tidy table of marker centroids (mean arcsinh expression)
#'   per matched population per batch. Use this sheet to visually identify
#'   each population's immunophenotype and fill in population_labels in config.
#' @param centroids_list Named list of centroid matrices [n_metaclusters x n_markers],
#'   from cluster_centroids.rds.
#' @param match_table data.frame from match_metaclusters().
#' @param markers Character vector of markers (column names of centroid matrices).
#' @return data.frame with columns: population, batch, cluster_id, <marker_1>, ...
format_centroid_report <- function(centroids_list, match_table, markers) {
  matched <- match_table[
    !is.na(match_table$matched_population_id) & match_table$is_matched == TRUE, ]

  rows <- lapply(seq_len(nrow(matched)), function(i) {
    row    <- matched[i, ]
    pop_id <- paste0("pop_", row$matched_population_id)

    c_a_df <- as.data.frame(
      centroids_list[[row$batch_a]][row$cluster_a, markers, drop = FALSE],
      check.names = FALSE)
    row_a  <- cbind(
      data.frame(population = pop_id, batch = row$batch_a,
                 cluster_id = row$cluster_a, stringsAsFactors = FALSE),
      c_a_df, row.names = NULL)

    c_b_df <- as.data.frame(
      centroids_list[[row$batch_b]][row$cluster_b, markers, drop = FALSE],
      check.names = FALSE)
    row_b  <- cbind(
      data.frame(population = pop_id, batch = row$batch_b,
                 cluster_id = row$cluster_b, stringsAsFactors = FALSE),
      c_b_df, row.names = NULL)

    rbind(row_a, row_b)
  })

  result <- do.call(rbind, rows)
  result[, markers] <- round(result[, markers], 4L)
  result
}
