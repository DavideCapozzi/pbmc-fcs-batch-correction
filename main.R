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

# Load essential libraries for the orchestrator
suppressPackageStartupMessages({
  library(yaml)
  library(tools)
  library(parallel)
})

# 2. Configuration & Logging Setup
# ------------------------------------------------------------------------------
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("[Main] FATAL: Config file not found at ", config_path)

# Source shared utilities before anything else — defines %||%, safe_div, etc.
source("src/functions/utils.R")
source("src/functions/config_validator.R")
source("src/functions/parallel_utils.R")

config <- validate_config(read_yaml(config_path))

# Define Output Structure
out_root <- config$directories$logs
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# Define run-specific log directory — isolates each execution for clinical audit
timestamp         <- format(Sys.time(), "%Y%m%d_%H%M%S")
PIPELINE_RUN_ID   <- timestamp
PIPELINE_LOGS_DIR <- file.path(out_root, paste0("run_", timestamp))
if (!dir.exists(PIPELINE_LOGS_DIR)) dir.create(PIPELINE_LOGS_DIR, recursive = TRUE)
log_file          <- file.path(PIPELINE_LOGS_DIR, paste0("pipeline_execution_", timestamp, ".log"))

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
  message(paste0("\n", strrep("=", 60)))
  message(sprintf(">>> STARTING STEP: %s", script_name))
  message(strrep("=", 60))
  
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
  
  message(sprintf("\n>>> SUCCESS: %s completed in %s seconds.", script_name, duration))
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

  # --- STEP 1: QC FILTERING ---
  run_pipeline_step("src/steps/01_qc.R")

  # --- STEP 2: LOAD & TRANSFORM (applies QC filters from Step 1) ---
  run_pipeline_step("src/steps/02_load.R")

  # --- CONDITIONAL BRANCH: MFI CORRECTION vs. FREQUENCY INTEGRATION ---
  run_per_batch <- isTRUE(config$clustering$run_per_batch)

  if (!run_per_batch) {
    message("\n[Main] clustering.run_per_batch = FALSE: running legacy MFI correction path.")

    # --- STEP 2: BATCH CORRECTION (LEGACY) ---
    run_pipeline_step("src/legacy/steps/02_correct.R")

    # --- STEP 2b: CORRECTION EVALUATION (LEGACY) ---
    run_pipeline_step("src/legacy/steps/02b_evaluate.R")

    # --- STEP 3: CLUSTERING ON COMBINED CORRECTED DATA (LEGACY) ---
    run_pipeline_step("src/legacy/steps/03_cluster.R")

  } else {
    message("\n[Main] clustering.run_per_batch = TRUE: running population frequency integration path.")

    # --- STEP 3: INDEPENDENT PER-BATCH CLUSTERING + CLUSTER MATCHING ---
    run_pipeline_step("src/steps/03_cluster_per_batch.R")

    # --- STEP 4: PER-SAMPLE POPULATION FREQUENCIES ---
    run_pipeline_step("src/steps/04_frequencies.R")

    # --- STEP 5: STRATIFIED ECOLOGICAL BASELINE ---
    run_pipeline_step("src/steps/05_integrate.R")

    # --- STEP 6: INTEGRATION VALIDATION (PEER-REVIEWER STATISTICS) ---
    run_pipeline_step("src/steps/06_evaluate_integration.R")
    
    # --- STEP 7: CLINICAL INFERENCE (Z-SCORE) ---
    if (!is.null(config$scoring) && isTRUE(config$scoring$run_patient_scoring)) {
      run_pipeline_step("src/steps/07_patient_scoring.R")
    }
  }

  # --- STEP 8: DIMRED: UMAP VISUALIZATION (BOTH PATHS) ---
  run_pipeline_step("src/steps/08_dimred.R")

  # Final Success Message
  message("\n==================================================")
  message(sprintf("PIPELINE FINISHED SUCCESSFULLY at %s", Sys.time()))
  message("==================================================")
  
}, error = function(e) {
  
  # Handle Critical Failures
  message("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("[FATAL ERROR] Pipeline stopped unexpectedly.")
  message("Error Message: ", e$message)
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  
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
