# src/functions/frequency.R
# ==============================================================================
# POPULATION FREQUENCY MODULE
# Description: Computes per-sample relative cell frequencies for each matched
#              metacluster population. The denominator is per-sample total cells
#              (all events, matched and unmatched), preserving the true proportion.
# Dependencies: SingleCellExperiment, dplyr
# ==============================================================================

library(SingleCellExperiment)
library(dplyr)

#' @title Compute Per-Sample Population Frequencies
#' @description For each sample and each matched population, computes
#'   freq_{s,p} = count_{s,p} / total_cells_s. Cells in unmatched clusters
#'   (matched_population_id = NA) are counted in the denominator but excluded
#'   from the population numerator.
#' @param sce_list Named list of per-batch SCE objects, each with colData columns:
#'   sample_id, batch, metacluster_id, matched_population_id.
#' @param match_table data.frame from match_metaclusters().
#' @return List with two data.frames:
#'   $freq_matrix  [n_samples x (2 + n_populations)]: sample_id, batch, pop_N columns (proportions)
#'   $count_matrix [n_samples x (2 + n_populations)]: same structure with raw counts
compute_population_frequencies <- function(sce_list, match_table) {

  matched_pop_ids <- sort(unique(match_table$matched_population_id[
    !is.na(match_table$matched_population_id)
  ]))

  pop_cols <- paste0("pop_", matched_pop_ids)

  all_freq   <- list()
  all_counts <- list()

  for (batch_name in names(sce_list)) {
    sce <- sce_list[[batch_name]]

    stopifnot(
      "sample_id"            %in% names(colData(sce)),
      "matched_population_id" %in% names(colData(sce))
    )

    cell_df <- data.frame(
      sample_id             = as.character(sce$sample_id),
      matched_population_id = as.character(sce$matched_population_id),
      stringsAsFactors      = FALSE
    )

    samples <- sort(unique(cell_df$sample_id))

    for (s in samples) {
      cells_s     <- cell_df[cell_df$sample_id == s, ]
      total_cells <- nrow(cells_s)

      count_row <- setNames(
        as.list(rep(0L, length(matched_pop_ids))),
        pop_cols
      )
      for (pid in matched_pop_ids) {
        count_row[[paste0("pop_", pid)]] <- sum(
          cells_s$matched_population_id == as.character(pid),
          na.rm = TRUE
        )
      }

      freq_row <- lapply(count_row, function(cnt) cnt / total_cells)

      all_counts[[length(all_counts) + 1L]] <- c(
        list(sample_id = s, batch = batch_name, total_cells = total_cells),
        count_row
      )
      all_freq[[length(all_freq) + 1L]] <- c(
        list(sample_id = s, batch = batch_name, total_cells = total_cells),
        freq_row
      )
    }
  }

  freq_matrix  <- as.data.frame(do.call(rbind, lapply(all_freq,   as.data.frame)))
  count_matrix <- as.data.frame(do.call(rbind, lapply(all_counts, as.data.frame)))

  for (col in pop_cols) {
    freq_matrix[[col]]  <- as.numeric(freq_matrix[[col]])
    count_matrix[[col]] <- as.integer(count_matrix[[col]])
  }

  list(freq_matrix = freq_matrix, count_matrix = count_matrix)
}
