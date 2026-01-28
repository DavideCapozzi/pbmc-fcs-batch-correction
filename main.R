# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR (RStudio Optimized)
# Description: Executes pipeline steps sequentially within the current session.
#              Captures real-time logs to file and console.
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
# Clean global environment to simulate a fresh start (Best effort isolation)
rm(list = ls()) 
graphics.off()

# Load essential libraries for the orchestrator
suppressPackageStartupMessages({
  library(yaml)
  library(tools)
})

# 2. Configuration & Logging Setup
# ------------------------------------------------------------------------------
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("[Main] FATAL: Config file not found at ", config_path)

config <- read_yaml(config_path)

# Define Output Structure
out_root <- config$directories$logs
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# Define Log File
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path(out_root, paste0("pipeline_execution_", timestamp, ".log"))

# 3. Helper Functions
# ------------------------------------------------------------------------------

#' @title Run Pipeline Step (Source Mode)
#' @description Sources an R script capturing output and handling errors.
run_pipeline_step <- function(script_rel_path) {
  
  if (!file.exists(script_rel_path)) {
    stop(paste("[Main] Script missing:", script_rel_path))
  }
  
  script_name <- basename(script_rel_path)
  
  # Log Header
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat(sprintf(">>> STARTING STEP: %s\n", script_name))
  cat(paste0(strrep("=", 60), "\n"))
  
  start_time <- Sys.time()
  
  # Execute Script via source()
  # 'local = FALSE' runs in GlobalEnv (needed for RStudio variable inspection)
  # 'echo = FALSE' prevents printing the code itself, keeps log clean
  tryCatch({
    source(script_rel_path, local = FALSE, echo = FALSE)
  }, error = function(e) {
    # Re-throw error with context
    stop(sprintf("Error in '%s': %s", script_name, e$message), call. = FALSE)
  })
  
  # Log Footer & Timing
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "secs"), 1)
  
  cat(sprintf("\n>>> SUCCESS: %s completed in %s seconds.\n", script_name, duration))
}

# 4. Execution Block
# ------------------------------------------------------------------------------
message(sprintf("[System] Pipeline started at %s", Sys.time()))
message(sprintf("[System] Logging to: %s", log_file))

# Start Logging (Split output to Console AND File)
con <- file(log_file, open = "wt")
sink(con, split = TRUE)      # Capture Standard Output
sink(con, type = "message")  # Capture Messages/Warnings (IMPORTANT for RStudio)

# Use tryCatch to ensure sink() is closed even if the pipeline crashes
tryCatch({
  
  # --- STEP 1: LOAD & TRANSFORM ---
  run_pipeline_step("src/steps/01_load.R")
  
  # --- STEP 2: BATCH CORRECTION ---
  run_pipeline_step("src/steps/02_correct.R")
  
  # --- STEP 3: CLUSTERING ---
  run_pipeline_step("src/steps/03_cluster.R")
  
  # --- STEP 4: CLUSTERING ---
  run_pipeline_step("src/steps/04_dimred.R")
  
  # Final Success Message
  cat("\n==================================================\n")
  cat(sprintf("PIPELINE FINISHED SUCCESSFULLY at %s\n", Sys.time()))
  cat("==================================================\n")
  
}, error = function(e) {
  
  # Handle Critical Failures
  cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("[FATAL ERROR] Pipeline stopped unexpectedly.\n")
  cat("Error Message: ", e$message, "\n")
  cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  
  # We don't stop() here because we want the 'finally' block to run cleanly
  
}, finally = {
  
  # 5. Safe Cleanup
  # ------------------------------------------------------------------------------
  # Restore output to console only
  sink(type = "message")
  sink()
  close(con)
  
  message("[System] Log file closed. Execution environment preserved for debugging.")
})
