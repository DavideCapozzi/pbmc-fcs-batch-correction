#!/usr/bin/env Rscript
# tests/smoke_regression.R
# ==============================================================================
# Compares current pipeline output digests against the frozen baseline.
# Run after refactoring to verify no numerical regressions were introduced.
#
# Usage (from project root, after running the pipeline with testing.enabled = true):
#   Rscript tests/smoke_regression.R
# ==============================================================================

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required. Install with: install.packages('digest')")
}

baseline_path <- "tests/baseline_digests.rds"
if (!file.exists(baseline_path)) {
  stop("Baseline not found. Run tests/freeze_baseline.R on the pre-refactor codebase first.")
}

baseline <- readRDS(baseline_path)
config   <- yaml::read_yaml("config/global_params.yml")

search_dirs <- c(
  config$directories$processed,
  config$directories$intermediate,
  config$directories$integration
)

current_files <- unlist(lapply(search_dirs, function(d) {
  if (!dir.exists(d)) return(character(0))
  list.files(d, pattern = "\\.(rds|xlsx|csv)$", full.names = TRUE, recursive = FALSE)
}))

if (length(current_files) == 0L) {
  stop("No current pipeline outputs found. Run Rscript main.R first.")
}

current_digests <- setNames(
  vapply(current_files, digest::digest, character(1L), file = TRUE, algo = "md5"),
  current_files
)

# Compare files that exist in both baseline and current run
common_files <- intersect(names(baseline), names(current_digests))
missing_from_current  <- setdiff(names(baseline), names(current_digests))
new_in_current        <- setdiff(names(current_digests), names(baseline))

n_pass <- 0L
n_fail <- 0L
fails  <- character(0L)

for (f in common_files) {
  if (identical(baseline[[f]], current_digests[[f]])) {
    n_pass <- n_pass + 1L
  } else {
    n_fail <- n_fail + 1L
    fails  <- c(fails, basename(f))
  }
}

message(sprintf("[Smoke] Compared %d files: %d PASS, %d FAIL", length(common_files), n_pass, n_fail))
if (length(missing_from_current) > 0L)
  message(sprintf("[Smoke] %d baseline files missing from current run: %s",
                  length(missing_from_current), paste(basename(missing_from_current), collapse = ", ")))
if (length(new_in_current) > 0L)
  message(sprintf("[Smoke] %d new files not in baseline: %s",
                  length(new_in_current), paste(basename(new_in_current), collapse = ", ")))

if (n_fail > 0L) {
  stop(sprintf("[Smoke] REGRESSION DETECTED in %d file(s): %s",
               n_fail, paste(fails, collapse = ", ")))
} else {
  message("[Smoke] All files match baseline — no regressions detected.")
}
