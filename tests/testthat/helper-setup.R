# tests/testthat/helper-setup.R
# Defines .proj_root (absolute path to project root) and .src() helper.
# Loaded by testthat before all test files.

# This file lives at <root>/tests/testthat/ so two levels up is the root.
.this_file <- normalizePath(
  if (!is.null(sys.frame(1L)$ofile)) sys.frame(1L)$ofile else "helper-setup.R"
)
.proj_root <- normalizePath(file.path(dirname(.this_file), "../.."))

# Source a file relative to the project root, in the global env so functions
# are accessible in subsequent test files.
.src <- function(...) source(file.path(.proj_root, ...), local = FALSE)
