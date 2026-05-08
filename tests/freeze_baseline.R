#!/usr/bin/env Rscript
# tests/freeze_baseline.R
# ==============================================================================
# Run this script ONCE on the main branch (before any refactoring) to capture
# digest fingerprints of all pipeline outputs. These fingerprints are used by
# tests/smoke_regression.R to verify that refactoring does not change results.
#
# Usage (from project root):
#   Rscript tests/freeze_baseline.R
# ==============================================================================

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required. Install with: install.packages('digest')")
}

config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("Run from project root: Rscript tests/freeze_baseline.R")

config <- yaml::read_yaml(config_path)

search_dirs <- c(
  config$directories$processed,
  config$directories$intermediate,
  config$directories$integration
)

all_files <- unlist(lapply(search_dirs, function(d) {
  if (!dir.exists(d)) return(character(0))
  list.files(d, pattern = "\\.(rds|xlsx|csv)$", full.names = TRUE, recursive = FALSE)
}))

if (length(all_files) == 0L) {
  stop("No pipeline output files found. Run the pipeline first (Rscript main.R).")
}

digests <- setNames(
  vapply(all_files, digest::digest, character(1L), file = TRUE, algo = "md5"),
  all_files
)

out_path <- "tests/baseline_digests.rds"
saveRDS(digests, out_path)
message(sprintf("[Freeze] Baseline captured: %d files → %s", length(digests), out_path))
for (f in head(names(digests), 10)) {
  message(sprintf("  %s  %s", digests[[f]], basename(f)))
}
if (length(digests) > 10) message(sprintf("  ... and %d more", length(digests) - 10))
