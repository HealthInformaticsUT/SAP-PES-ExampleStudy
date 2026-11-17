# Connect to your OMOP database first and name the CDMConnector reference cdm.
# Optional: supply nihss_table with person_id, measurement_date and nihss_score,
# or place nihss.csv beside this file. Without scores, NIHSS analyses are skipped.
# Run this script with source("/path/to/SAP-PES-ExampleStudy/code2Run.R").

study_result <- local({
  script <- tail(unlist(lapply(sys.frames(), function(frame) frame$ofile)), 1L)
  if (!length(script)) stop("Run this file using source().")
  study_dir <- dirname(normalizePath(script, mustWork = TRUE))
  if (!exists("cdm", inherits = TRUE)) stop("Create your connected CDMConnector reference named cdm first.")
  scores <- if (exists("nihss_table", inherits = TRUE)) get("nihss_table", inherits = TRUE) else {
    path <- file.path(study_dir, "nihss.csv")
    if (file.exists(path)) path else NULL
  }
  for (file in c("assumptions.R", file.path("R",
    c("data.R", "trajectories.R", "observed.R", "survival.R", "economics.R", "export.R", "run.R")))) {
    sys.source(file.path(study_dir, file), envir = environment())
  }
  run_study(cdm = cdm, nihss = scores, config = study_settings(study_dir))
})

# Protected tables are also available in study_result$tables.
# study_result$archive is the ZIP to return to the study team.
