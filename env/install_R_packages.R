# Force standard libcurl for secure downloads via Conda certificates
options(download.file.method = "libcurl")

# Bypass potential GitHub API rate limits or user-agent blocks
options(HTTPUserAgent = sprintf(
  "R/%s R (%s)", 
  getRversion(), 
  paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
))

# Install cyCombine strictly without updating dependencies
message("Starting strict installation of cyCombine...")

remotes::install_github(
  repo = "biosurf/cyCombine",
  dependencies = FALSE,
  upgrade = "never",
  build_vignettes = FALSE,
  force = TRUE
)

# Verification step
if (requireNamespace("cyCombine", quietly = TRUE)) {
  message("SUCCESS: cyCombine has been correctly installed and loaded.")
  library(cyCombine)
} else {
  stop("CRITICAL ERROR: cyCombine installation failed.")
}