# src/functions/step_logger.R
# ==============================================================================
# STEP LOGGER: JSON-structured per-step execution audit trail
# Produces one JSON file per pipeline step in results/logs/steps/
# Functions operate on a mutable environment — no re-assignment needed.
# ==============================================================================

library(jsonlite)

#' Initialize a step log environment.
#' @param step_name character(1) e.g. "03_cluster_per_batch"
#' @param step_number integer(1) ordinal position in the pipeline
#' @param input_files character vector of paths consumed by this step
#' @return environment to pass to add_metric / finalize_step_log / write_step_json
init_step_log <- function(step_name, step_number, input_files) {
  now   <- Sys.time()
  log_e <- new.env(parent = emptyenv())

  log_e$step        <- as.character(step_name)
  log_e$step_number <- as.integer(step_number)
  log_e$run_id      <- if (exists("PIPELINE_RUN_ID", envir = globalenv(), inherits = FALSE)) {
    get("PIPELINE_RUN_ID", envir = globalenv())
  } else {
    format(now, "%Y%m%d_%H%M%S")
  }
  log_e$started_at   <- format(now, "%Y-%m-%dT%H:%M:%S%z")
  log_e$input_files  <- as.character(input_files)
  log_e$metrics      <- list()
  log_e$status       <- "RUNNING"
  log_e$.start_time  <- now

  log_e
}

#' Add a named metric to the step log (mutates in-place).
#' @param log_obj environment from init_step_log()
#' @param key character(1) metric name
#' @param value ANY — will be serialized by jsonlite
#' @return invisible(log_obj)
add_metric <- function(log_obj, key, value) {
  stopifnot(is.environment(log_obj), is.character(key), length(key) == 1L)
  log_obj$metrics[[key]] <- value
  invisible(log_obj)
}

#' Finalize timing, outputs, and status (mutates in-place).
#' @param log_obj environment from init_step_log()
#' @param output_files character vector of paths produced by this step
#' @param status "SUCCESS" or "FAILURE"
#' @return invisible(log_obj)
finalize_step_log <- function(log_obj, output_files, status = "SUCCESS") {
  now <- Sys.time()
  log_obj$finished_at      <- format(now, "%Y-%m-%dT%H:%M:%S%z")
  log_obj$duration_seconds <- round(as.numeric(difftime(now, log_obj$.start_time, units = "secs")), 2)
  log_obj$output_files     <- as.character(output_files)
  log_obj$status           <- as.character(status)
  invisible(log_obj)
}

#' Serialize the log to a JSON file in logs_dir/steps/.
#' File name: {NN}_{step_name}_{run_id}.json
#' @param log_obj finalized environment
#' @param logs_dir character(1) root log directory (e.g. "results/logs")
#' @return invisible character(1) path of the written file
write_step_json <- function(log_obj, logs_dir) {
  effective_logs_dir <- if (exists("PIPELINE_LOGS_DIR", envir = globalenv(), inherits = FALSE)) {
    get("PIPELINE_LOGS_DIR", envir = globalenv())
  } else {
    logs_dir
  }
  steps_dir <- file.path(effective_logs_dir, "steps")
  if (!dir.exists(steps_dir)) dir.create(steps_dir, recursive = TRUE)

  fname <- sprintf("%02d_%s_%s.json",
                   log_obj$step_number,
                   log_obj$step,
                   log_obj$run_id)
  out_path <- file.path(steps_dir, fname)

  log_list <- list(
    step             = log_obj$step,
    step_number      = log_obj$step_number,
    run_id           = log_obj$run_id,
    started_at       = log_obj$started_at,
    finished_at      = if (!is.null(log_obj$finished_at)) log_obj$finished_at else NA_character_,
    duration_seconds = if (!is.null(log_obj$duration_seconds)) log_obj$duration_seconds else NA_real_,
    status           = log_obj$status,
    inputs           = log_obj$input_files,
    outputs          = if (!is.null(log_obj$output_files)) log_obj$output_files else character(0L),
    metrics          = log_obj$metrics
  )

  jsonlite::write_json(log_list, path = out_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  message(sprintf("[Logger] Step log written: %s", out_path))
  invisible(out_path)
}
