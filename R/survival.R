# Cox analyses use exact source dates, rather than weekly state labels.

cox_concepts <- function() {
  data.frame(
    concept_id = c(132797L, 257315L, 1344143L, 3046279L, 4052413L,
                   4052415L, 4116186L, 4085730L, 4095048L, 4189544L,
                   4260440L, 4301474L),
    domain = c("Condition", "Condition", "Drug", "Measurement", "Procedure",
               "Procedure", "Procedure", "Procedure", "Measurement",
               "Measurement", "Measurement", "Observation"),
    variable = c("prior_sepsis", "prior_bacterial_pneumonia", "prior_albumin",
                 "prior_procalcitonin", rep("prior_central_venous_access", 3),
                 "prior_antibiotic_therapy", "prior_acid_base_evaluation",
                 "prior_microbial_culture", "prior_immunochromatographic_test",
                 "prior_invasive_arterial_monitoring"),
    stringsAsFactors = FALSE
  )
}

cox_collect <- function(records) {
  if (!inherits(records, "data.frame")) records <- dplyr::collect(records)
  as.data.frame(records, stringsAsFactors = FALSE)
}

check_survival_inputs <- function(patients) {
  patients <- as.data.frame(patients, stringsAsFactors = FALSE)
  required <- c("subject_id", "index_date", "observation_end", "pneumonia_date",
                "death_date", "age", "gender_concept_id", "eligible")
  if (!all(required %in% names(patients))) stop("Cox patient data lack: ",
    paste(setdiff(required, names(patients)), collapse = ", "))
  patients$subject_id <- as.character(patients$subject_id)
  if (!nrow(patients) || anyNA(patients$subject_id) || anyDuplicated(patients$subject_id)) {
    stop("Cox analysis requires one target episode per patient, with nonmissing IDs.")
  }
  for (column in c("index_date", "observation_end", "pneumonia_date", "death_date")) {
    patients[[column]] <- as.Date(patients[[column]])
  }
  if (anyNA(patients$index_date) || anyNA(patients$observation_end) || anyNA(patients$eligible) ||
      any(patients$observation_end < patients$index_date)) stop("Invalid Cox observation or eligibility data.")
  for (column in c("pneumonia_date", "death_date")) {
    if (any(patients[[column]] < patients$index_date | patients[[column]] > patients$observation_end, na.rm = TRUE)) {
      stop(column, " lies outside the target observation interval.")
    }
  }
  if (any(patients$pneumonia_date > patients$death_date, na.rm = TRUE)) {
    stop("Pneumonia after recorded death must be resolved before Cox analysis.")
  }
  patients
}

# Exact same-index-date NIHSS. Identical repeated scores count once.
# Conflicting or invalid same-day values exclude the patient from the subset.
# Records without patient IDs cannot be linked and are counted separately.
cox_link_nihss <- function(patients, nihss_table = NULL) {
  patients$nihss_score <- NA_real_
  patients$nihss_status <- "No same-day score"
  audit <- data.frame(metric = c("patients", "usable_same_day", "conflicting_same_day",
    "invalid_same_day", "no_same_day_score", "exact_duplicate_rows",
    "input_records", "records_missing_patient_id", "records_with_patient_id"),
    count = c(nrow(patients), 0L, 0L, 0L, nrow(patients), 0L, 0L, 0L, 0L))
  if (is.null(nihss_table)) {
    return(list(patients = patients, audit = audit, status = "NIHSS_NOT_SUPPLIED"))
  }
  score_records <- as.data.frame(nihss_table, stringsAsFactors = FALSE)
  patient_id <- if ("person_id" %in% names(score_records)) "person_id" else "subject_id"
  required <- c(patient_id, "measurement_date", "nihss_score")
  if (!all(required %in% names(score_records))) {
    stop("nihss_table must contain person_id (or subject_id), measurement_date and nihss_score.")
  }
  if (!nrow(score_records)) return(list(patients = patients, audit = audit, status = "NIHSS_NOT_SUPPLIED"))
  input_records <- nrow(score_records)
  raw_id <- as.character(score_records[[patient_id]])
  missing_id <- is.na(raw_id) | trimws(raw_id) %in% c("", "NA", "NULL")
  records_missing_id <- sum(missing_id)
  score_records <- score_records[!missing_id, , drop = FALSE]
  score_records <- data.frame(subject_id = as.character(score_records[[patient_id]]),
    measurement_date = as.Date(score_records$measurement_date),
    nihss_score = suppressWarnings(as.numeric(as.character(score_records$nihss_score))))
  if (anyNA(score_records$measurement_date)) {
    stop("NIHSS records with a patient ID have missing dates. Check the source values and date format.")
  }
  duplicates <- nrow(score_records) - nrow(unique(score_records))
  score_records <- unique(score_records)
  position <- match(score_records$subject_id, patients$subject_id)
  same <- !is.na(position) & score_records$measurement_date == patients$index_date[position]
  score_records <- score_records[which(same), , drop = FALSE]
  for (patient_id in unique(score_records$subject_id)) {
    scores <- score_records$nihss_score[score_records$subject_id == patient_id]
    valid <- is.finite(scores) & scores >= 0 & scores <= 42 & scores == round(scores)
    patient_row <- match(patient_id, patients$subject_id)
    if (!all(valid)) patients$nihss_status[patient_row] <- "Invalid same-day score"
    else if (length(unique(scores)) > 1L) patients$nihss_status[patient_row] <- "Conflicting same-day scores"
    else {
      patients$nihss_status[patient_row] <- "Usable same-day score"
      patients$nihss_score[patient_row] <- scores[1L]
    }
  }
  audit$count <- c(nrow(patients), sum(!is.na(patients$nihss_score)),
    sum(patients$nihss_status == "Conflicting same-day scores"),
    sum(patients$nihss_status == "Invalid same-day score"),
    sum(patients$nihss_status == "No same-day score"), duplicates,
    input_records, records_missing_id, input_records - records_missing_id)
  list(patients = patients, audit = audit, status = "SUPPLIED")
}

# Historical baseline variables are reproduced exactly, without descendant
# expansion. All flags refer to the 365 days before the target visit starts.
prepare_cox_inputs <- function(cdm, patients, nihss_table = NULL) {
  patients <- check_survival_inputs(patients)
  concepts <- cox_concepts()
  domains <- data.frame(domain = c("Condition", "Drug", "Measurement", "Procedure", "Observation"),
    table = c("condition_occurrence", "drug_exposure", "measurement", "procedure_occurrence", "observation"),
    concept = c("condition_concept_id", "drug_concept_id", "measurement_concept_id", "procedure_concept_id", "observation_concept_id"),
    date = c("condition_start_date", "drug_exposure_start_date", "measurement_date", "procedure_date", "observation_date"))
  ids <- patients$subject_id
  min_date <- min(patients$index_date) - 365L
  max_date <- max(patients$index_date) - 1L
  for (variable in unique(concepts$variable)) patients[[variable]] <- 0L
  for (domain_number in seq_len(nrow(domains))) {
    domain <- domains[domain_number, ]
    if (!domain$table %in% names(cdm)) stop("Cox baseline extraction requires cdm$", domain$table, ".")
    records <- cdm[[domain$table]]
    if (!all(c("person_id", domain$concept, domain$date) %in% colnames(records))) {
      stop("Missing baseline columns in cdm$", domain$table, ".")
    }
    concept_ids <- concepts$concept_id[concepts$domain == domain$domain]
    records <- records |>
      dplyr::filter(.data$person_id %in% !!ids,
        .data[[domain$concept]] %in% !!concept_ids,
        .data[[domain$date]] >= !!min_date, .data[[domain$date]] <= !!max_date) |>
      dplyr::transmute(subject_id = as.character(.data$person_id),
        concept_id = .data[[domain$concept]], record_date = .data[[domain$date]]) |>
      dplyr::distinct()
    records <- cox_collect(records)
    position <- match(records$subject_id, patients$subject_id)
    dates <- as.Date(records$record_date)
    records <- records[which(dates >= patients$index_date[position] - 365L & dates < patients$index_date[position]), , drop = FALSE]
    for (variable in unique(concepts$variable[concepts$domain == domain$domain])) {
      hits <- records$subject_id[records$concept_id %in% concepts$concept_id[concepts$variable == variable]]
      patients[[variable]][patients$subject_id %in% hits] <- 1L
    }
  }
  if (!"visit_occurrence" %in% names(cdm)) stop("Cox baseline extraction requires cdm$visit_occurrence.")
  acute <- cdm$visit_occurrence |>
    dplyr::filter(.data$person_id %in% !!ids, .data$visit_concept_id %in% c(262L, 9201L, 9203L),
      .data$visit_start_date >= !!min_date, .data$visit_start_date <= !!max_date) |>
    dplyr::transmute(subject_id = as.character(.data$person_id),
      visit_id = as.character(.data$visit_occurrence_id), record_date = .data$visit_start_date) |>
    dplyr::distinct()
  acute <- cox_collect(acute)
  position <- match(acute$subject_id, patients$subject_id)
  acute <- acute[which(as.Date(acute$record_date) >= patients$index_date[position] - 365L &
                         as.Date(acute$record_date) < patients$index_date[position]), , drop = FALSE]
  counts <- table(acute$subject_id)
  patients$prior_acute_visit_count <- as.integer(counts[match(patients$subject_id, names(counts))])
  patients$prior_acute_visit_count[is.na(patients$prior_acute_visit_count)] <- 0L
  patients$age_per_10 <- patients$age / 10
  patients$male <- ifelse(patients$gender_concept_id == 8507, 1,
    ifelse(patients$gender_concept_id == 8532, 0, NA_real_))
  patients$index_year <- as.integer(format(patients$index_date, "%Y"))
  patients$index_year_centered <- patients$index_year - stats::median(patients$index_year)
  patients$log1p_prior_acute_visits <- log1p(patients$prior_acute_visit_count)
  linked <- cox_link_nihss(patients, nihss_table)
  score_rows <- list()
  for (group in c("All target", "Eligible", "Not eligible")) {
    group_patients <- linked$patients
    if (group == "Eligible") group_patients <- group_patients[group_patients$eligible, , drop = FALSE]
    if (group == "Not eligible") group_patients <- group_patients[!group_patients$eligible, , drop = FALSE]
    available <- group_patients$nihss_score[!is.na(group_patients$nihss_score)]
    quantiles <- rep(NA_real_, 3)
    if (length(available)) quantiles <- stats::quantile(available, c(.25, .5, .75), names = FALSE)
    patient_count <- nrow(group_patients)
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      group = group, patients = patient_count, available = length(available),
      unavailable = patient_count - length(available),
      conflicting = sum(group_patients$nihss_status == "Conflicting same-day scores"),
      invalid = sum(group_patients$nihss_status == "Invalid same-day score"),
      percent_available = if (patient_count) 100 * length(available) / patient_count else NA_real_,
      median = quantiles[2L], q25 = quantiles[1L], q75 = quantiles[3L],
      status = if (length(available)) "AVAILABLE" else "NIHSS_UNAVAILABLE",
      .support_n = length(available))
  }
  score_summary <- dplyr::bind_rows(score_rows)
  list(patients = linked$patients, tables = list(nihss_linkage = linked$audit,
    nihss_summary = score_summary, cox_covariate_concepts = concepts))
}

cox_bands <- function() c("unexposed", "0-6 days", "7-27 days", "28-89 days", "90-364 days", "365+ days")

# Each calendar day contributes risk through its end. A Day-d death has tstop
# d+1, so a one-day exposure lag leaves a same-day diagnosis/death unexposed.
# Band elapsed time starts when lagged exposure starts, consistent with the exposure definition.
build_cox_intervals <- function(patients, horizon_days = Inf, lag_days = 1L, banded = FALSE) {
  if (!lag_days %in% c(0L, 1L)) stop("This analysis defines lag_days as 0 or 1.")
  band_boundaries <- c(0, 7, 28, 90, 365, Inf)
  band_labels <- cox_bands()
  if (!nrow(patients)) {
    patients$tstart <- patients$tstop <- numeric()
    patients$event <- patients$pneumonia_history <- integer()
    patients$exposure <- factor(character(), levels = band_labels)
    return(patients)
  }
  patient_intervals <- list()
  for (patient_number in seq_len(nrow(patients))) {
    patient <- patients[patient_number, , drop = FALSE]
    followup_end <- patient$observation_end
    if (is.finite(horizon_days)) followup_end <- min(followup_end, patient$index_date + horizon_days - 1L)
    died <- !is.na(patient$death_date) && patient$death_date <= followup_end
    if (died) followup_end <- patient$death_date
    stop_day <- as.integer(followup_end - patient$index_date) + 1L
    if (stop_day < 1L) next
    exposure_day <- Inf
    if (!is.na(patient$pneumonia_date)) {
      exposure_day <- as.integer(patient$pneumonia_date - patient$index_date) + lag_days
    }
    starts <- 0
    if (exposure_day > 0 && exposure_day < stop_day) starts <- c(starts, exposure_day)
    if (banded && exposure_day < stop_day) starts <- c(starts, exposure_day + band_boundaries[2:5])
    starts <- sort(unique(starts[starts >= 0 & starts < stop_day]))
    ends <- c(starts[-1L], stop_day)
    intervals <- patient[rep(1L, length(starts)), , drop = FALSE]
    intervals$tstart <- starts
    intervals$tstop <- ends
    intervals$event <- as.integer(died & ends == stop_day)
    intervals$pneumonia_history <- as.integer(starts >= exposure_day)
    intervals$exposure <- "unexposed"
    if (banded) {
      exposed <- starts >= exposure_day
      band_number <- findInterval(starts[exposed] - exposure_day, band_boundaries) + 1L
      intervals$exposure[exposed] <- band_labels[band_number]
    }
    patient_intervals[[patient_number]] <- intervals
  }
  intervals <- dplyr::bind_rows(patient_intervals)
  if (nrow(intervals)) intervals$exposure <- factor(intervals$exposure, levels = band_labels)
  intervals
}

cox_support <- function(long, banded = FALSE) {
  exposure_levels <- c("unexposed", "exposed")
  exposure <- ifelse(long$pneumonia_history == 1L, "exposed", "unexposed")
  if (banded) {
    exposure_levels <- cox_bands()
    exposure <- as.character(long$exposure)
  }
  rows <- list()
  for (level in exposure_levels) {
    intervals <- long[exposure == level, , drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(exposure = level,
      patients = length(unique(intervals$subject_id)), deaths = sum(intervals$event),
      person_days = sum(intervals$tstop - intervals$tstart), stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}

fit_cox_model <- function(long, rhs, robust = FALSE) {
  if (!nrow(long) || sum(long$event) < 2L) {
    return(list(fit = NULL, status = "INSUFFICIENT_EVENTS", warnings = character()))
  }
  warnings <- character()
  fit <- tryCatch(withCallingHandlers({
    formula <- stats::as.formula(paste("survival::Surv(tstart, tstop, event) ~", rhs))
    args <- list(formula = formula, data = long, ties = "efron", x = TRUE,
      model = TRUE, control = survival::coxph.control(iter.max = 100L))
    if (robust) {
      args$cluster <- long$subject_id
      args$robust <- TRUE
    }
    do.call(survival::coxph, args)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    warnings <<- c(warnings, conditionMessage(e))
    NULL
  })
  list(fit = fit, status = if (is.null(fit)) "FIT_FAILED" else if (length(warnings)) "FIT_WARNING" else "FIT",
       warnings = unique(warnings))
}

cox_coefficient <- function(result, term, model, support_n, support_events) {
  row <- data.frame(model = model, term = term, hazard_ratio = NA_real_, ci95_low = NA_real_,
    ci95_high = NA_real_, p_value = NA_real_, status = result$status,
    .support_n = support_n, .support_events = support_events, stringsAsFactors = FALSE)
  if (result$status != "FIT") return(row)
  fit <- result$fit
  coefficient <- stats::coef(fit)
  if (!term %in% names(coefficient) || !is.finite(coefficient[[term]])) {
    row$status <- "TERM_NOT_ESTIMABLE"
    return(row)
  }
  standard_error <- sqrt(stats::vcov(fit)[term, term])
  values <- c(exp(coefficient[[term]]), exp(coefficient[[term]] - stats::qnorm(.975) * standard_error),
    exp(coefficient[[term]] + stats::qnorm(.975) * standard_error), 2 * stats::pnorm(-abs(coefficient[[term]] / standard_error)))
  if (!all(is.finite(values)) || any(values[1:3] <= 0)) {
    row$status <- "TERM_NOT_ESTIMABLE"
    return(row)
  }
  row[c("hazard_ratio", "ci95_low", "ci95_high", "p_value")] <- as.list(values)
  row
}

cox_ph_tests <- function(result, model, support_n, support_events) {
  empty <- data.frame(model = model, term = "GLOBAL", chisq = NA_real_, df = NA_real_,
    p_value = NA_real_, status = result$status, .support_n = support_n,
    .support_events = support_events)
  if (result$status != "FIT") return(empty)
  ph_test <- tryCatch(survival::cox.zph(result$fit, transform = "km", terms = TRUE)$table,
    error = function(e) NULL)
  if (is.null(ph_test)) {
    empty$status <- "PH_TEST_FAILED"
    return(empty)
  }
  data.frame(model = model, term = rownames(ph_test), chisq = ph_test[, "chisq"],
    df = ph_test[, "df"], p_value = ph_test[, "p"], status = "FIT", .support_n = support_n,
    .support_events = support_events, row.names = NULL)
}

compare_cox_time_effects <- function(reduced, full, model, support_n, support_events) {
  out <- data.frame(model = model, test = "banded_vs_single_pneumonia_effect", statistic = NA_real_,
    df = NA_real_, p_value = NA_real_, status = "NOT_ESTIMABLE",
    .support_n = support_n, .support_events = support_events)
  if (reduced$status != "FIT" || full$status != "FIT") return(out)
  # Both fits use the exact same rows, risk sets and adjustment variables.
  degrees_of_freedom <- sum(is.finite(stats::coef(full$fit))) - sum(is.finite(stats::coef(reduced$fit)))
  statistic <- 2 * (tail(full$fit$loglik, 1) - tail(reduced$fit$loglik, 1))
  if (degrees_of_freedom <= 0 || !is.finite(statistic) || statistic < -1e-7) return(out)
  out$statistic <- max(0, statistic)
  out$df <- degrees_of_freedom
  out$p_value <- stats::pchisq(out$statistic, degrees_of_freedom, lower.tail = FALSE)
  out$status <- "FIT"
  out
}

# landmark_eligible is supplied explicitly because a visit-interval phenotype
# may only become available after Day 7. Never infer prospective availability
# from the retrospective eligibility flag.
analyse_survival <- function(patients, landmark_eligible = NULL,
                              landmark_policy = "NOT_SUPPLIED") {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival for the Cox analyses.")
  patients <- check_survival_inputs(patients)
  required <- c("age_per_10", "male", "index_year_centered", "log1p_prior_acute_visits",
                unique(cox_concepts()$variable))
  if (!all(required %in% names(patients))) stop("Run prepare_cox_inputs() before analyse_survival().")
  if (!"nihss_score" %in% names(patients)) patients$nihss_score <- NA_real_
  patients$nihss5 <- patients$nihss_score / 5
  patients$landmark_eligible <- if (is.null(landmark_eligible)) NA else landmark_eligible
  if (!is.null(landmark_eligible) && (!is.logical(landmark_eligible) ||
      length(landmark_eligible) != nrow(patients) || anyNA(landmark_eligible))) {
    stop("landmark_eligible must supply one nonmissing logical flag per patient.")
  }
  varying <- missing_columns <- character()
  for (variable in required) {
    observed_values <- patients[[variable]][!is.na(patients[[variable]])]
    if (!length(observed_values)) missing_columns <- c(missing_columns, variable)
    if (length(unique(observed_values)) > 1L) varying <- c(varying, variable)
  }
  # Exclude missing adjustment values before fitting. A completely missing
  # variable must not disappear from the requested adjustment set.
  if (length(missing_columns)) stop("No usable values for adjustment variable(s): ", paste(missing_columns, collapse = ", "))
  complete_adjustment <- stats::complete.cases(patients[required])
  adjusted_patients <- patients[complete_adjustment, , drop = FALSE]
  estimates <- ph <- supports <- summaries <- lrts <- list()
  diagnostics <- list()
  adjustment <- data.frame(variable = required,
    status = ifelse(required %in% varying, "INCLUDED", "OMITTED_CONSTANT"))
  for (lag in c(0L, 1L)) {
    model_name <- if (lag == 0L) "adjusted_365_day" else "adjusted_365_day_lag1"
    long <- build_cox_intervals(adjusted_patients, 365L, lag, FALSE)
    support <- cox_support(long)
    rhs <- paste(c("pneumonia_history", varying), collapse = " + ")
    fitted <- fit_cox_model(long, rhs, robust = TRUE)
    support_patients <- min(support$patients)
    support_deaths <- min(support$deaths)
    estimates[[model_name]] <- cox_coefficient(fitted, "pneumonia_history", model_name, support_patients, support_deaths)
    ph[[model_name]] <- cox_ph_tests(fitted, model_name, support_patients, support_deaths)
    support$model <- model_name
    supports[[model_name]] <- support
    summaries[[model_name]] <- data.frame(model = model_name, supplied_patients = nrow(patients),
      model_patients = nrow(adjusted_patients), excluded_missing_adjustment = sum(!complete_adjustment),
      deaths = sum(long$event), recorded_pneumonia = sum(!is.na(adjusted_patients$pneumonia_date) &
        adjusted_patients$pneumonia_date <= adjusted_patients$index_date + 364L),
      same_day_pneumonia_deaths = sum(adjusted_patients$pneumonia_date == adjusted_patients$death_date &
        adjusted_patients$death_date <= adjusted_patients$index_date + 364L, na.rm = TRUE))
    diagnostics[[model_name]] <- data.frame(model = model_name, status = fitted$status,
      warning = paste(fitted$warnings, collapse = " | "))
  }
  for (kind in c("full", "nihss_unadjusted", "nihss_adjusted")) {
    model_name <- paste0("bands_lag1_", kind)
    model_patients <- if (kind == "full") patients else patients[!is.na(patients$nihss_score), , drop = FALSE]
    long <- build_cox_intervals(model_patients, Inf, 1L, TRUE)
    support <- cox_support(long, TRUE)
    suffix <- if (kind == "nihss_adjusted") " + nihss5" else ""
    fitted <- fit_cox_model(long, paste0("exposure", suffix))
    binary <- fit_cox_model(long, paste0("pneumonia_history", suffix))
    if (!nrow(model_patients)) fitted$status <- binary$status <- "NIHSS_UNAVAILABLE"
    support_patients <- min(support$patients)
    support_deaths <- min(support$deaths)
    band_rows <- list()
    for (band in cox_bands()[-1L]) {
      compared_groups <- support$exposure %in% c("unexposed", band)
      row <- cox_coefficient(fitted, paste0("exposure", band), model_name,
        min(support$patients[compared_groups]), min(support$deaths[compared_groups]))
      row$band <- band
      band_rows[[length(band_rows) + 1L]] <- row
    }
    estimates[[model_name]] <- dplyr::bind_rows(band_rows)
    if (kind == "nihss_adjusted") estimates[[paste0(model_name, "_nihss")]] <-
      cox_coefficient(fitted, "nihss5", model_name, support_patients, support_deaths)
    ph[[model_name]] <- cox_ph_tests(fitted, model_name, support_patients, support_deaths)
    lrts[[model_name]] <- compare_cox_time_effects(binary, fitted, model_name, support_patients, support_deaths)
    if (!nrow(model_patients)) lrts[[model_name]]$status <- "NIHSS_UNAVAILABLE"
    support$model <- model_name
    supports[[model_name]] <- support
    summaries[[model_name]] <- data.frame(model = model_name, supplied_patients = nrow(patients),
      model_patients = nrow(model_patients), deaths = sum(long$event))
    diagnostics[[model_name]] <- data.frame(model = model_name, status = fitted$status,
      warning = paste(fitted$warnings, collapse = " | "))
  }
  landmark <- cox_landmark(patients, supplied = !is.null(landmark_eligible), policy = landmark_policy)
  tables <- list(cox_estimates = dplyr::bind_rows(estimates), cox_model_support = dplyr::bind_rows(supports),
    cox_model_summary = dplyr::bind_rows(summaries), cox_ph_tests = dplyr::bind_rows(ph),
    cox_band_comparison = dplyr::bind_rows(lrts), cox_model_diagnostics = dplyr::bind_rows(diagnostics),
    cox_adjustment_variables = adjustment)
  tables <- c(tables, landmark)
  tables
}

cox_landmark <- function(patients, supplied, policy) {
  if (!supplied) return(list(landmark_cox_status = data.frame(
    status = "LANDMARK_ELIGIBILITY_NOT_SUPPLIED", policy = policy)))
  keep <- patients$observation_end >= patients$index_date + 8L &
    (is.na(patients$death_date) | patients$death_date > patients$index_date + 7L) &
    (is.na(patients$pneumonia_date) | patients$pneumonia_date > patients$index_date + 7L)
  landmark_patients <- patients[keep, , drop = FALSE]
  flow <- data.frame(group = c("All target", "Alive and pneumonia-free with follow-up on Day 8",
    "Landmark eligible", "NIHSS complete", "NIHSS complete and landmark eligible"),
    patients = c(nrow(patients), nrow(landmark_patients), sum(landmark_patients$landmark_eligible), sum(!is.na(landmark_patients$nihss_score)),
                 sum(!is.na(landmark_patients$nihss_score) & landmark_patients$landmark_eligible)))
  rows <- supports <- ph <- diagnostics <- list()
  for (outcome in c("pneumonia", "death")) {
    for (adjustment in c("full_unadjusted", "nihss_unadjusted", "nihss_adjusted")) {
      model_name <- paste("day8", outcome, adjustment, sep = "_")
      model_patients <- if (adjustment == "full_unadjusted") landmark_patients else landmark_patients[!is.na(landmark_patients$nihss_score), , drop = FALSE]
      start <- model_patients$index_date + 8L
      end <- model_patients$observation_end
      died <- !is.na(model_patients$death_date) & model_patients$death_date <= end
      end[died] <- model_patients$death_date[died]
      event_date <- if (outcome == "pneumonia") model_patients$pneumonia_date else model_patients$death_date
      event <- !is.na(event_date) & event_date >= start & event_date <= end
      end[event] <- event_date[event]
      model_patients$tstart <- rep(0, nrow(model_patients))
      model_patients$tstop <- as.integer(end - start) + 1L
      model_patients$event <- as.integer(event)
      model_patients$direct <- as.integer(model_patients$landmark_eligible)
      # For pneumonia, deaths before a pneumonia record censor follow-up.
      # Same-day pneumonia/death is counted as a pneumonia event.
      support <- data.frame(model = model_name, outcome = outcome,
        exposure = c("reference", "eligible"),
        patients = c(sum(model_patients$direct == 0), sum(model_patients$direct == 1)),
        events = c(sum(model_patients$event[model_patients$direct == 0]), sum(model_patients$event[model_patients$direct == 1])))
      fit <- fit_cox_model(model_patients, if (adjustment == "nihss_adjusted") "direct + nihss5" else "direct")
      if (!nrow(model_patients)) fit$status <- if (nrow(landmark_patients)) "NIHSS_UNAVAILABLE" else "NO_LANDMARK_FOLLOWUP"
      diagnostics[[model_name]] <- data.frame(model = model_name, status = fit$status,
        warning = paste(fit$warnings, collapse = " | "))
      support_patients <- min(support$patients)
      support_events <- min(support$events)
      rows[[model_name]] <- cox_coefficient(fit, "direct", model_name, support_patients, support_events)
      if (adjustment == "nihss_adjusted") rows[[paste0(model_name, "_nihss")]] <-
        cox_coefficient(fit, "nihss5", model_name, support_patients, support_events)
      supports[[model_name]] <- support
      ph[[model_name]] <- cox_ph_tests(fit, model_name, support_patients, support_events)
    }
  }
  estimates <- dplyr::bind_rows(rows)
  support_table <- dplyr::bind_rows(supports)
  status <- if (!nrow(landmark_patients)) "NO_LANDMARK_FOLLOWUP" else if (!nrow(estimates) ||
    !any(estimates$status == "FIT")) "NO_ESTIMABLE_MODELS" else "RUN"
  list(landmark_flow = flow, landmark_cox_estimates = estimates,
    landmark_cox_support = support_table, landmark_cox_ph = dplyr::bind_rows(ph),
    landmark_cox_diagnostics = dplyr::bind_rows(diagnostics),
    landmark_cox_status = data.frame(status = status, policy = policy))
}
