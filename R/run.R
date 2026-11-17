# Read this function first: it follows the study from cohorts to result files.
run_study <- function(cdm, nihss = NULL, config) {
  for (package in c("Cohort2Trajectory", "TrajectoryMarkovAnalysis", "CDMConnector",
                    "CirceR", "SqlRender", "dplyr", "dbplyr", "jsonlite", "survival", "zip")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " before running this study.")
  }
  if (utils::packageVersion("Cohort2Trajectory") < "1.3.0" ||
      utils::packageVersion("TrajectoryMarkovAnalysis") < "1.1.0") {
    stop("This study requires Cohort2Trajectory >= 1.3.0 and TrajectoryMarkovAnalysis >= 1.1.0.")
  }
  if (!inherits(cdm, "cdm_reference")) stop("Supply a connected CDMConnector cdm reference.")
  nihss <- read_nihss(nihss)

  # C2T creates temporary patient files. They are deleted when this run ends.
  work_dir <- tempfile("tma_study_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  message("1/6 Generating cohorts and weekly trajectories.")
  cohorts <- generate_study_cohorts(cdm, config, basename(work_dir))
  cdm <- cohorts$cdm
  trajectories <- create_weekly_trajectories(cohorts$rows, work_dir)
  patients <- add_baseline_information(cdm, get_patient_dates(trajectories), trajectories$eligibility_cohort)
  markers <- add_dysphagia_markers(cdm, patients, config)
  cox_inputs <- prepare_cox_inputs(cdm, markers$patients, nihss)
  patients <- cox_inputs$patients
  if (!nrow(nihss)) message("No NIHSS table supplied. NIHSS models will be reported as unavailable.")

  message("2/6 Estimating transitions and costs.")
  cost_patients <- patients$subject_id[patients$eligible | patients$broad_eligible]
  cost_histories <- TrajectoryMarkovAnalysis::prepareTrajectories(
    trajectories$trajectories[trajectories$trajectories$subject_id %in% cost_patients, ],
    cycleLength = 7L,
    states = c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory"),
    absorbingStates = c("DeathNoHistory", "DeathHistory"))
  cost_records <- TrajectoryMarkovAnalysis::extractCosts(cdm, trajectories = cost_histories,
    measure = "total_charge", domains = "Visit Detail", allocation = "duration")
  cost_records$subject_id <- as.character(cost_records$subject_id)
  eligible_ids <- patients$subject_id[patients$eligible]
  eligible_cost_records <- cost_records[cost_records$subject_id %in% eligible_ids, ]
  primary_result <- estimate_model_inputs(trajectories, patients, config, eligible_cost_records)

  message("3/6 Re-estimating the three cohort and chronology sensitivities.")
  definitions <- c(
    strict_chronology = "Same eligible patients. Same-day first pneumonia/death is not prior pneumonia history. Other pneumonia dates stay unchanged.",
    broad_marker = "Swallowing/care marker roots and descendants recorded during Days 0-7. No pre-index carry-forward.",
    day8_landmark = "Fixed dysphagia group, alive and pneumonia-free through Day 7 and observed on Day 8. Follow-up and four-week effect begin Day 8."
  )
  variants <- list()
  for (scenario in names(definitions)) {
    scenario_patients <- apply_sensitivity_definition(patients, scenario)
    if (!any(scenario_patients$eligible)) next
    directory <- file.path(work_dir, scenario)
    dir.create(directory)
    scenario_cohorts <- cohort_rows_from_patients(scenario_patients)
    scenario_trajectories <- create_weekly_trajectories(scenario_cohorts, directory)
    eligible_ids <- scenario_patients$subject_id[scenario_patients$eligible]
    scenario_costs <- cost_records[cost_records$subject_id %in% eligible_ids, ]
    scenario_result <- estimate_model_inputs(scenario_trajectories, scenario_patients, config, scenario_costs)
    scenario_result$variant_note <- unname(definitions[scenario])
    scenario_result$costs$fit <- NULL
    variants[[scenario]] <- scenario_result
  }

  message("4/6 Summarising observed care, survival validation and Cox models.")
  observed <- summarise_observed_care(trajectories, patients, config, primary_result,
    config$validation_replicates, config$seed)
  survival <- analyse_survival(patients, landmark_eligible = patients$landmark_eligible,
    landmark_policy = "Fixed dysphagia group. Day 8 survivors without prior pneumonia.")

  message("5/6 Projecting economic scenarios and simulation intervals.")
  economics <- analyse_cost_effectiveness(primary_result, config, variants,
    n_simulations = config$economic_replicates, seed = config$seed + 1L)
  model_inputs <- study_input_reports(primary_result, markers, config)
  results <- c(observed, cox_inputs$tables, survival, economics, model_inputs)

  message("6/6 Exporting protected result tables and figure data.")
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir <- tempfile(paste0("study_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = config$output_dir)
  exported <- export_results(results, config, file.path(run_dir, "shareable"))
  archive <- file.path(normalizePath(run_dir), "shareable_results.zip")
  zip::zipr(archive, list.files(exported$export_dir, recursive = TRUE),
    root = exported$export_dir, include_directories = FALSE, mode = "mirror")
  message("Complete. Share only: ", archive)
  invisible(list(archive = archive, export_dir = exported$export_dir, tables = exported$tables))
}

study_input_reports <- function(primary_result, markers, config) {
  schedule <- primary_result$schedule
  schedule$time_years <- schedule$cycle * 7 / config$economic$days_per_year
  for (probability in c("background_probability", "p_death_no_pneumonia", "p_death_pneumonia")) {
    schedule[[paste0(probability, "_annual_equivalent")]] <-
      -expm1(log1p(-schedule[[probability]]) * config$economic$days_per_year / 7)
  }
  list(eligibility_definitions = markers$counts,
       mortality_calibration = primary_result$mortality$calibration,
       mortality_schedule = schedule)
}
