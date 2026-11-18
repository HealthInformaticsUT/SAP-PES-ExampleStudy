# Cohort definitions, CDM inputs and weekly C2T construction.

# Optional NIHSS input stays in memory. CSV identifiers are always read as text.
# Score, date and linkage validation belongs to cox_link_nihss().
read_nihss <- function(x = NULL) {
  if (is.null(x)) {
    return(data.frame(person_id = character(), measurement_date = as.Date(character()),
      nihss_score = numeric(), stringsAsFactors = FALSE))
  }
  if (inherits(x, "data.frame")) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    if (!file.exists(x)) stop("NIHSS CSV does not exist: ", x)
    x <- utils::read.csv(x, colClasses = "character", stringsAsFactors = FALSE)
  } else {
    stop("Supply NIHSS as a data frame, a CSV path, or NULL.")
  }
  required <- c("person_id", "measurement_date", "nihss_score")
  if (!all(required %in% names(x))) {
    stop("NIHSS input needs person_id, measurement_date and nihss_score columns.")
  }
  x$person_id <- as.character(x$person_id)
  x
}

collect_table <- function(table) {
  as.data.frame(dplyr::collect(table))
}

# Include the stated roots and their valid OMOP descendants, as ATLAS does.
expand_concepts <- function(cdm, roots) {
  valid_concepts <- cdm$concept |>
    dplyr::filter(is.na(.data$invalid_reason)) |>
    dplyr::select("concept_id")
  descendants <- cdm$concept_ancestor |>
    dplyr::filter(.data$ancestor_concept_id %in% !!roots) |>
    dplyr::semi_join(valid_concepts, by = c("descendant_concept_id" = "concept_id")) |>
    dplyr::select("descendant_concept_id") |>
    collect_table()
  unique(c(as.integer(roots), as.integer(descendants$descendant_concept_id)))
}

generate_study_cohorts <- function(cdm, config, run_id) {
  definitions <- CDMConnector::readCohortSet(config$cohort_dir)
  table_name <- paste0(config$cohort_table_prefix, "_", run_id)
  cdm <- CDMConnector::generateCohortSet(
    cdm, definitions, name = table_name, computeAttrition = FALSE, overwrite = FALSE
  )
  cohorts <- cdm[[table_name]]
  target_ids <- cohorts |>
    dplyr::filter(.data$cohort_definition_id == 1L) |>
    dplyr::distinct(.data$subject_id)

  # Keep all target/eligibility rows so episode mismatches remain detectable.
  # Pneumonia and death are needed only for members of the target cohort.
  relevant <- dplyr::union_all(
    cohorts |> dplyr::filter(.data$cohort_definition_id %in% c(1L, 2L)),
    cohorts |> dplyr::filter(.data$cohort_definition_id %in% c(3L, 4L)) |>
      dplyr::semi_join(target_ids, by = "subject_id")
  )
  rows <- relevant |>
    dplyr::transmute(
      cohort_definition_id = as.integer(.data$cohort_definition_id),
      subject_id = as.character(.data$subject_id),
      cohort_start_date = .data$cohort_start_date,
      cohort_end_date = .data$cohort_end_date
    ) |>
    collect_table() |>
    dplyr::distinct()
  rows$cohort_start_date <- as.Date(rows$cohort_start_date)
  rows$cohort_end_date <- as.Date(rows$cohort_end_date)
  list(cdm = cdm, rows = rows, table_name = table_name)
}

define_state_intervals <- function(cohort_rows) {
  cohorts <- unique(as.data.frame(cohort_rows))
  cohorts$subject_id <- as.character(cohorts$subject_id)
  cohorts$cohort_start_date <- as.Date(cohorts$cohort_start_date)
  cohorts$cohort_end_date <- as.Date(cohorts$cohort_end_date)
  if (anyNA(cohorts) || any(cohorts$cohort_end_date < cohorts$cohort_start_date)) {
    stop("Generated cohorts contain missing values or reversed dates.")
  }
  patients <- cohorts |>
    dplyr::filter(.data$cohort_definition_id == 1L) |>
    dplyr::transmute(subject_id = .data$subject_id, index_date = .data$cohort_start_date,
      observation_end = .data$cohort_end_date)
  if (!nrow(patients) || anyDuplicated(patients$subject_id)) {
    stop("Expected one stroke target episode per patient.")
  }

  # Use events overlapping the target episode. A pre-index event that carries
  # into this episode is an error, rather than an onset on the index date.
  events <- cohorts |>
    dplyr::filter(.data$cohort_definition_id %in% c(3L, 4L)) |>
    dplyr::inner_join(patients, by = "subject_id") |>
    dplyr::filter(.data$cohort_end_date >= .data$index_date,
      .data$cohort_start_date <= .data$observation_end)
  if (any(events$cohort_start_date < events$index_date)) {
    stop("A pneumonia or death event starts before the target episode and carries into follow-up.")
  }
  pneumonia <- events |>
    dplyr::filter(.data$cohort_definition_id == 3L) |>
    dplyr::arrange(.data$cohort_start_date) |>
    dplyr::distinct(.data$subject_id, .keep_all = TRUE) |>
    dplyr::transmute(subject_id = .data$subject_id, pneumonia_date = .data$cohort_start_date)
  deaths <- events |>
    dplyr::filter(.data$cohort_definition_id == 4L) |>
    dplyr::distinct(.data$subject_id, death_date = .data$cohort_start_date)
  if (anyDuplicated(deaths$subject_id)) stop("Conflicting recorded death dates within target follow-up.")
  patients <- patients |>
    dplyr::left_join(pneumonia, by = "subject_id") |>
    dplyr::left_join(deaths, by = "subject_id")
  if (any(patients$pneumonia_date > patients$death_date, na.rm = TRUE)) {
    stop("Pneumonia entry occurs after recorded death.")
  }

  eligibility <- cohorts[cohorts$cohort_definition_id == 2L, ]
  target_position <- match(eligibility$subject_id, patients$subject_id)
  if (anyNA(target_position) || anyDuplicated(eligibility$subject_id) ||
      any(eligibility$cohort_start_date != patients$index_date[target_position]) ||
      any(eligibility$cohort_end_date != patients$observation_end[target_position])) {
    stop("Generated dysphagia eligibility does not match the stroke target episodes.")
  }

  target <- dplyr::transmute(patients, subject_id = .data$subject_id,
    cohort_definition_id = "0", cohort_start_date = .data$index_date,
    cohort_end_date = .data$observation_end)
  no_pneumonia <- target
  no_pneumonia$cohort_definition_id <- "AliveNoPneum"
  pneumonia <- patients |>
    dplyr::filter(!is.na(.data$pneumonia_date)) |>
    dplyr::transmute(subject_id = .data$subject_id, cohort_definition_id = "AlivePneum",
      cohort_start_date = .data$pneumonia_date,
      cohort_end_date = pmin(.data$observation_end, .data$death_date, na.rm = TRUE))
  deaths <- patients |>
    dplyr::filter(!is.na(.data$death_date)) |>
    dplyr::transmute(subject_id = .data$subject_id,
      cohort_definition_id = dplyr::if_else(is.na(.data$pneumonia_date), "DeathNoHistory", "DeathHistory"),
      cohort_start_date = .data$death_date, cohort_end_date = .data$observation_end)
  intervals <- dplyr::bind_rows(target, no_pneumonia, pneumonia,
    deaths[deaths$cohort_definition_id == "DeathNoHistory", ],
    deaths[deaths$cohort_definition_id == "DeathHistory", ])
  list(intervals = intervals, eligibility = eligibility, patients = patients)
}

create_weekly_trajectories <- function(cohort_rows, run_dir, config) {
  facts <- define_state_intervals(cohort_rows)
  states <- c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory")
  study <- Cohort2Trajectory::cohort2TrajectoryConfiguration(
    studyName = "c2t", pathToStudy = run_dir,
    atlasTargetCohort = "0", useCDM = FALSE, runSavedStudy = FALSE,
    trajectoryDataObject = facts$intervals,
    stateCohortLabels = states, stateCohortPriorityOrder = rev(states),
    stateCohortAbsorbing = states[3:4], stateSelectionType = "Priority",
    trajectoryType = "Discrete", cycleLength = config$economic$cycle_days,
    outOfCohortAllowed = FALSE, mergeStates = FALSE
  )
  Cohort2Trajectory::getDataForStudy(studyEnv = study)
  trajectories <- Cohort2Trajectory::createTrajectories(studyEnv = study, useCDM = FALSE)
  trajectories$subject_id <- as.character(trajectories$subject_id)

  # Exact dates are retained because one weekly label cannot show an onset
  # followed by death in the same week.
  list(trajectories = trajectories, source_intervals = facts$intervals,
    eligibility_cohort = facts$eligibility, patients = facts$patients)
}

get_patient_dates <- function(c2t) {
  patients <- c2t$patients
  patients$pneumonia_day <- as.integer(patients$pneumonia_date - patients$index_date)
  patients$death_day <- as.integer(patients$death_date - patients$index_date)
  patients
}

add_baseline_information <- function(cdm, patients, eligibility_cohort) {
  ids <- patients$subject_id
  people <- cdm$person |>
    dplyr::filter(.data$person_id %in% !!ids) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      year_of_birth = .data$year_of_birth, gender_concept_id = .data$gender_concept_id) |>
    collect_table()
  if (anyDuplicated(people$subject_id) || !setequal(people$subject_id, ids)) {
    stop("The person table must supply exactly one record for every target patient.")
  }
  patients <- dplyr::left_join(patients, people, by = "subject_id")
  patients$age <- as.integer(format(patients$index_date, "%Y")) - as.integer(patients$year_of_birth)
  patients$eligible <- patients$subject_id %in% eligibility_cohort$subject_id
  patients
}

# Read the frozen broad-marker concepts across the original clinical domains.
# These are raw recorded occurrences. No pre-index history is carried forward.
add_dysphagia_markers <- function(cdm, patients, config) {
  ids <- as.character(patients$subject_id)
  roots <- utils::read.csv(config$broad_marker_file)$source_concept_id
  concepts <- expand_concepts(cdm, roots)
  direct <- expand_concepts(cdm, 31317L)
  condition_records <- cdm$condition_occurrence |>
    dplyr::filter(.data$person_id %in% !!ids, .data$condition_concept_id %in% !!concepts) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      concept_id = .data$condition_concept_id, record_date = .data$condition_start_date,
      visit_id = as.character(.data$visit_occurrence_id)) |>
    collect_table()
  procedure_records <- cdm$procedure_occurrence |>
    dplyr::filter(.data$person_id %in% !!ids, .data$procedure_concept_id %in% !!concepts) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      concept_id = .data$procedure_concept_id, record_date = .data$procedure_date,
      visit_id = as.character(.data$visit_occurrence_id)) |>
    collect_table()
  observations <- cdm$observation |>
    dplyr::filter(.data$person_id %in% !!ids, .data$observation_concept_id %in% !!concepts) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      concept_id = .data$observation_concept_id, record_date = .data$observation_date,
      visit_id = as.character(.data$visit_occurrence_id)) |>
    collect_table()
  measurements <- cdm$measurement |>
    dplyr::filter(.data$person_id %in% !!ids, .data$measurement_concept_id %in% !!concepts) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      concept_id = .data$measurement_concept_id, record_date = .data$measurement_date,
      visit_id = as.character(.data$visit_occurrence_id)) |>
    collect_table()
  records <- dplyr::bind_rows(
    condition_occurrence = condition_records, procedure_occurrence = procedure_records,
    observation = observations, measurement = measurements, .id = "domain"
  )
  records$record_date <- as.Date(records$record_date)
  records$direct_dysphagia <- records$domain == "condition_occurrence" & records$concept_id %in% direct
  visits <- cdm$visit_occurrence |>
    dplyr::filter(.data$person_id %in% !!ids) |>
    dplyr::transmute(subject_id = as.character(.data$person_id), visit_id = as.character(.data$visit_occurrence_id),
      visit_concept_id = .data$visit_concept_id, visit_start = .data$visit_start_date,
      visit_end = .data$visit_end_date) |> collect_table()
  visits$visit_start <- as.Date(visits$visit_start)
  visits$visit_end <- as.Date(visits$visit_end)
  if (anyDuplicated(visits[c("subject_id", "visit_id")])) stop("Duplicate patient/visit keys in visit_occurrence.")
  # Resolve the same stroke visit selected by the target definition.
  definition <- jsonlite::fromJSON(file.path(config$cohort_dir, "target.json"), simplifyVector = FALSE)
  source_roots <- vapply(definition$ConceptSets[[1]]$expression$items, function(z) z$concept$CONCEPT_ID, numeric(1))
  source_codes <- expand_concepts(cdm, source_roots)
  acute <- expand_concepts(cdm, c(262L, 9201L, 9203L))
  stroke_links <- cdm$condition_occurrence |>
    dplyr::filter(.data$person_id %in% !!ids, .data$condition_source_concept_id %in% !!source_codes) |>
    dplyr::transmute(subject_id = as.character(.data$person_id), visit_id = as.character(.data$visit_occurrence_id)) |>
    dplyr::distinct() |> collect_table()
  index_visits <- visits |>
    dplyr::inner_join(stroke_links, by = c("subject_id", "visit_id")) |>
    dplyr::inner_join(patients[c("subject_id", "index_date")], by = "subject_id") |>
    dplyr::filter(.data$visit_concept_id %in% acute, .data$visit_start == .data$index_date) |>
    dplyr::add_count(.data$subject_id, name = "matching_visits") |>
    dplyr::filter(.data$matching_visits == 1L, !is.na(.data$visit_end),
      .data$visit_end >= .data$visit_start) |>
    dplyr::transmute(subject_id = .data$subject_id,
      stroke_visit_start = .data$visit_start, stroke_visit_end = .data$visit_end)
  linked <- records |>
    dplyr::left_join(visits[c("subject_id", "visit_id", "visit_start", "visit_end")], by = c("subject_id", "visit_id")) |>
    dplyr::left_join(patients[c("subject_id", "index_date", "observation_end")], by = "subject_id") |>
    dplyr::left_join(index_visits[c("subject_id", "stroke_visit_start", "stroke_visit_end")], by = "subject_id")
  linked$qualifying_visit <- with(linked, !is.na(visit_start) & !is.na(visit_end) &
    !is.na(stroke_visit_end) & visit_end >= visit_start &
    visit_end >= stroke_visit_start & visit_start <= stroke_visit_end + 7L)
  linked$qualifying_visit[is.na(linked$qualifying_visit)] <- FALSE
  # This cross-check prevents a second, subtly different eligibility definition.
  reproduced <- unique(linked$subject_id[linked$direct_dysphagia & linked$qualifying_visit])
  if (!setequal(reproduced, patients$subject_id[patients$eligible])) {
    stop("Visit-linked dysphagia reconstruction disagrees with the generated cohort. No landmark or broad-visit analysis was run.")
  }
  recorded_by_day7 <- !is.na(linked$record_date) & linked$record_date <= linked$index_date + 7L
  known <- unique(linked$subject_id[
    linked$direct_dysphagia & linked$qualifying_visit & recorded_by_day7
  ])
  patients$landmark_eligible <- patients$eligible
  broad <- recorded_by_day7 & linked$record_date >= linked$index_date &
    linked$record_date <= linked$observation_end
  patients$broad_eligible <- patients$subject_id %in% unique(linked$subject_id[broad])
  any_direct <- linked$direct_dysphagia & !is.na(linked$record_date) &
    linked$record_date >= linked$index_date & linked$record_date <= linked$observation_end
  patients$direct_dysphagia_any_followup <- patients$subject_id %in% unique(linked$subject_id[any_direct])
  list(patients = patients,
    counts = data.frame(target_patients = nrow(patients), primary_eligible = sum(patients$eligible),
      qualifying_dysphagia_recorded_by_day7 = length(known), broad_eligible = sum(patients$broad_eligible),
      direct_dysphagia_any_followup = sum(patients$direct_dysphagia_any_followup)))
}

cohort_rows_from_patients <- function(patients) {
  target <- dplyr::transmute(patients,
    cohort_definition_id = 1L, subject_id = .data$subject_id,
    cohort_start_date = .data$index_date, cohort_end_date = .data$observation_end)
  eligibility <- target[patients$eligible, ]
  eligibility$cohort_definition_id <- rep(2L, nrow(eligibility))
  pneumonia <- patients |>
    dplyr::filter(!is.na(.data$pneumonia_date)) |>
    dplyr::transmute(cohort_definition_id = 3L, subject_id = .data$subject_id,
      cohort_start_date = .data$pneumonia_date,
      cohort_end_date = pmin(.data$observation_end, .data$death_date, na.rm = TRUE))
  deaths <- patients |>
    dplyr::filter(!is.na(.data$death_date)) |>
    dplyr::transmute(cohort_definition_id = 4L, subject_id = .data$subject_id,
      cohort_start_date = .data$death_date, cohort_end_date = .data$death_date)
  dplyr::bind_rows(target, eligibility, pneumonia, deaths)
}

# Rebuild variants from dates, never relabel existing weekly rows by guesswork.
apply_sensitivity_definition <- function(patients, variant) {
  if (variant == "strict_chronology") {
    same <- !is.na(patients$pneumonia_date) & !is.na(patients$death_date) &
      patients$pneumonia_date >= patients$death_date
    patients$pneumonia_date[same] <- as.Date(NA)
  } else if (variant == "broad_marker") {
    patients$eligible <- patients$broad_eligible
  } else if (variant == "day8_landmark") {
    cutoff <- patients$index_date + 7L
    keep <- patients$observation_end >= cutoff + 1L &
      (is.na(patients$death_date) | patients$death_date > cutoff) &
      (is.na(patients$pneumonia_date) | patients$pneumonia_date > cutoff)
    patients <- patients[keep, , drop = FALSE]
    patients$eligible <- patients$landmark_eligible
    patients$age <- patients$age + 8 / 365.25
    patients$index_date <- patients$index_date + 8L
  } else {
    stop("Unknown variant: ", variant)
  }
  patients$pneumonia_day <- as.integer(patients$pneumonia_date - patients$index_date)
  patients$death_day <- as.integer(patients$death_date - patients$index_date)
  patients
}
