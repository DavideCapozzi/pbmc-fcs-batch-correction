# src/functions/scoring.R
# ==============================================================================
# PATIENT INFERENCE & SCORING MODULE
# Description: Applies the stratified baseline to new patient frequency data to 
#              compute Z-Scores and clinical status flags.
# Dependencies: dplyr, tidyr
# ==============================================================================

library(dplyr)
library(tidyr)

#' @title Compute Patient Z-Scores
#' @description Calculates Z-scores for patient populations using the matched 
#'              batch-specific baseline.
#' @param patient_freq_matrix data.frame containing patient frequencies.
#' @param baseline_dict data.frame containing baseline mean and sd.
#' @return data.frame with Z-scores and clinical flags.
compute_patient_zscores <- function(patient_freq_matrix, baseline_dict, z_threshold = 2.0) {
  
  stopifnot(
    "batch" %in% colnames(patient_freq_matrix),
    "sample_id" %in% colnames(patient_freq_matrix)
  )
  
  # Pivot patient data to long format
  patient_long <- patient_freq_matrix %>%
    dplyr::select(sample_id, batch, tidyselect::starts_with("pop_")) %>%
    tidyr::pivot_longer(
      cols      = tidyselect::starts_with("pop_"),
      names_to  = "population",
      values_to = "patient_freq"
    )
  
  # Join with baseline strictly by batch AND population
  scored_data <- patient_long %>%
    dplyr::inner_join(
      baseline_dict %>% dplyr::select(batch, population, mean_freq, sd_freq),
      by = c("batch", "population"),
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      # Protect against division by zero (undefined variance in baseline)
      z_score = dplyr::if_else(
        sd_freq > 1e-6, 
        (patient_freq - mean_freq) / sd_freq, 
        NA_real_
      ),
      clinical_flag = dplyr::case_when(
        is.na(z_score) ~ "UNDETERMINED",
        z_score >= z_threshold ~ "EXPANDED",
        z_score <= -z_threshold ~ "DEPLETED",
        TRUE ~ "NORMAL"
      )
    )
  
  return(as.data.frame(scored_data))
}

#' @title Format Clinical Report
#' @description Merges biological labels and formats the scoring output for clinicians.
format_clinical_report <- function(scored_data, match_table, population_labels = NULL) {
  
  pop_id_map <- match_table %>%
    dplyr::filter(is_matched == TRUE, !is.na(matched_population_id)) %>%
    dplyr::select(matched_population_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(population = paste0("pop_", matched_population_id))
  
  report <- scored_data %>%
    dplyr::inner_join(pop_id_map, by = "population")
  
  if (!is.null(population_labels) && length(population_labels) > 0) {
    label_df <- data.frame(
      population      = names(population_labels),
      population_name = unlist(population_labels, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    report <- dplyr::left_join(report, label_df, by = "population")
    
    col_order <- c("sample_id", "batch", "population", "population_name", 
                   "patient_freq", "mean_freq", "sd_freq", "z_score", "clinical_flag")
    report <- dplyr::select(report, dplyr::any_of(col_order))
  }
  
  # Ensure clean formatting for Excel
  numeric_cols <- c("patient_freq", "mean_freq", "sd_freq", "z_score")
  report <- report %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(numeric_cols), ~ round(.x, 4L)))
  
  return(as.data.frame(report))
}