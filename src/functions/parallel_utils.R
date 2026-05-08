# src/functions/parallel_utils.R
# ==============================================================================
# PARALLEL EXECUTION UTILITIES
# Provides a single pipeline_lapply() wrapper that transparently switches
# between parallel::mclapply (Linux/WSL2) and sequential lapply.
#
# Worker count is resolved by get_n_workers(config) which respects:
#   1. config$parallel$n_workers (explicit override in YAML)
#   2. config$n_cores             (legacy top-level key)
#   3. parallel::detectCores() - 1 (auto-detect upper bound)
#
# Set env var PIPELINE_SEQUENTIAL=1 to force sequential execution for debugging.
# ==============================================================================

#' Determine the number of parallel workers to use.
#' @param config list from read_yaml() (validated by validate_config)
#' @return integer >= 1
get_n_workers <- function(config) {
  requested  <- config$parallel$n_workers %||% config$n_cores %||% 1L
  requested  <- as.integer(requested)
  available  <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  n_workers  <- min(requested, available)
  message(sprintf("[Parallel] Workers: %d (requested %d, available %d)",
                  n_workers, requested, available))
  n_workers
}

#' Parallelised lapply wrapper.
#'
#' Runs FUN over X using parallel::mclapply on Linux/WSL2 when n_workers > 1.
#' Falls back to sequential lapply when:
#'   - n_workers <= 1
#'   - env var PIPELINE_SEQUENTIAL is set to "1"
#'
#' Output order is always preserved (mclapply guarantees this).
#' Worker failures are surfaced as a single stop() after collection.
#'
#' @param X object to iterate over (as in lapply)
#' @param FUN function to apply
#' @param n_workers integer, number of parallel processes
#' @return list of results, same length and order as X
pipeline_lapply <- function(X, FUN, n_workers = 1L) {
  force_sequential <- identical(Sys.getenv("PIPELINE_SEQUENTIAL"), "1")

  if (force_sequential || n_workers <= 1L) {
    return(lapply(X, FUN))
  }

  results <- parallel::mclapply(X, FUN,
                                 mc.cores       = n_workers,
                                 mc.preschedule = FALSE,
                                 mc.set.seed    = TRUE)

  is_error <- vapply(results, inherits, logical(1L), "try-error")
  if (any(is_error)) {
    err_msgs <- vapply(results[is_error], as.character, character(1L))
    stop(sprintf("[Parallel] %d worker(s) failed:\n%s",
                 sum(is_error),
                 paste(head(err_msgs, 3L), collapse = "\n")),
         call. = FALSE)
  }

  results
}
