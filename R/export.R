# Disclosure rules live here so the analysis code only calculates results.
study_release_rules <- function(tables) {
  counts <- list(
    demographics = c("patients", "age_available", "age_missing", "female", "male", "sex_other_or_missing"),
    observed_population = c("patients", "source_state_rows", "observed_patient_cycles", "shortened_patient_cycles",
      "pneumonia_full_followup", "deaths_full_followup", "censored_alive", "observed_beyond_365_days",
      "pneumonia_by_day365", "deaths_by_day365", "same_day_pneumonia_death", "represented_cycle_days", "observed_alive_days"),
    transition_estimates = c("events", "non_events", "denominator", "source_patients"),
    observed_occupancy = c("patients", "baseline_patients"),
    survival_validation = c("patients", "observed_deaths"),
    weekly_charges = c("patient_cycles", "observed_days"),
    nihss_linkage = "count",
    nihss_summary = c("patients", "available", "unavailable", "conflicting", "invalid"),
    cox_model_support = c("patients", "deaths", "person_days"),
    cox_model_summary = c("supplied_patients", "model_patients", "excluded_missing_adjustment",
      "deaths", "recorded_pneumonia", "same_day_pneumonia_deaths"),
    landmark_flow = "patients",
    landmark_cox_support = c("patients", "events"),
    eligibility_definitions = names(tables$eligibility_definitions)
  )
  rules <- list()
  for (name in names(tables)) {
    group <- "cohort_characterisation"
    if (name %in% c("transition_estimates", "fitted_occupancy", "survival_validation", "mortality_schedule")) group <- "empirical_model"
    if (name == "observed_occupancy") group <- "observed_occupancy"
    if (name == "weekly_charges") group <- "costs"
    if (name %in% c("cox_covariate_concepts", "cox_adjustment_variables")) group <- "public_specification"
    if (startsWith(name, "economic_") || name == "discussion_comparisons") {
      group <- "economics"
      counts[[name]] <- intersect("eligible_n", names(tables[[name]]))
    }
    rules[[name]] <- list(group = group, counts = counts[[name]],
      support = grep("^\\.support", names(tables[[name]]), value = TRUE), extra = numeric())
  }

  # Whole-cohort and NIHSS-subset tables can reveal the omitted subgroup by subtraction.
  support <- tables$cox_model_support
  whole <- support[support$model == "bands_lag1_full", ]
  subset <- support[support$model == "bands_lag1_nihss_unadjusted", ]
  if (nrow(whole) && nrow(subset)) {
    subset <- subset[match(whole$exposure, subset$exposure), ]
    rules$cox_model_support$extra <- unlist(
      whole[c("patients", "deaths", "person_days")] - subset[c("patients", "deaths", "person_days")])
  }
  flow <- tables$landmark_flow$patients
  if (length(flow) == 5L) {
    total <- flow[1]
    survivors <- flow[2]
    eligible <- flow[3]
    with_nihss <- flow[4]
    eligible_with_nihss <- flow[5]
    rules$landmark_flow$extra <- c(total - survivors, survivors - eligible,
      survivors - with_nihss, eligible - eligible_with_nihss, with_nihss - eligible_with_nihss,
      (survivors - eligible) - (with_nihss - eligible_with_nihss))
  }
  support <- tables$landmark_cox_support
  for (outcome in c("pneumonia", "death")) {
    whole <- support[support$model == paste0("day8_", outcome, "_full_unadjusted"), ]
    subset <- support[support$model == paste0("day8_", outcome, "_nihss_unadjusted"), ]
    if (nrow(whole) && nrow(subset)) {
      subset <- subset[match(whole$exposure, subset$exposure), ]
      rules$landmark_cox_support$extra <- c(rules$landmark_cox_support$extra,
        unlist(whole[c("patients", "events")] - subset[c("patients", "events")]))
    }
  }
  rules
}

number_text <- function(values) {
  if (!is.numeric(values)) return(as.character(values))
  text <- rep("NA", length(values))
  known <- is.finite(values)
  text[known] <- format(values[known], digits = 15, scientific = FALSE, trim = TRUE)
  text
}

suppress_small_cells <- function(tables, rules) {
  is_small <- function(values) is.finite(values) & values > 0 & values < 5
  forbidden <- c("subject_id", "person_id", "simulation_id", "index_date", "pneumonia_date",
    "death_date", "observation_end", "observation_end_date", "measurement_date", "visit_id",
    "visit_occurrence_id", "state_start_date", "state_end_date", "record_date")
  unsafe_groups <- character()
  small_population <- FALSE

  # First identify groups of linked results that contain a protected count.
  for (name in names(tables)) {
    table <- tables[[name]]
    rule <- rules[[name]]
    if (!is.data.frame(table) || any(names(table) %in% forbidden) ||
        any(vapply(table, function(column) is.list(column) || inherits(column, "Date") || inherits(column, "POSIXt"), logical(1)))) {
      stop("Only flat aggregate tables without patient identifiers or calendar dates may be exported.")
    }
    columns <- c(rule$counts, rule$support)
    if (nrow(table) && !all(columns %in% names(table))) stop("Missing count/support column in ", name)
    columns <- intersect(columns, names(table))
    if (!all(vapply(table[columns], is.numeric, logical(1)))) stop("Counts and support must be numeric.")
    values <- c(unlist(table[columns]), rule$extra)
    for (pair in list(c("events", "cycles"), c("deaths", "cycles"), c("events", "denominator"),
                      c("deaths", "denominator"), c("events", "patients"), c("deaths", "patients"),
                      c("deaths", "model_patients"), c("recorded_pneumonia", "model_patients"))) {
      if (all(pair %in% columns)) values <- c(values, table[[pair[2]]] - table[[pair[1]]])
    }
    if (any(is_small(values))) unsafe_groups <- union(unsafe_groups, rule$group)
    population_columns <- intersect(c("eligible_n", "eligible_patients", "primary_eligible"), names(table))
    small_population <- small_population || any(is_small(unlist(table[population_columns])))
  }
  if (any(c("empirical_model", "costs") %in% unsafe_groups)) unsafe_groups <- union(unsafe_groups, "economics")
  if (small_population) unsafe_groups <- union(unsafe_groups, c("economics", "cohort_characterisation", "empirical_model"))

  # Then replace small counts by <5 and withhold their linked numeric results.
  output <- list()
  for (name in names(tables)) {
    table <- tables[[name]]
    rule <- rules[[name]]
    table <- table[setdiff(names(table), rule$support)]
    for (column in names(table)) {
      values <- table[[column]]
      if (is.numeric(values) && rule$group %in% unsafe_groups) {
        table[[column]] <- rep("SUPPRESSED", length(values))
        if (column %in% rule$counts) table[[column]][is_small(values)] <- "<5"
      } else if (is.numeric(values)) {
        table[[column]] <- ifelse(is.na(values), "NA", format(values, digits = 12, trim = TRUE, scientific = FALSE))
      } else {
        table[[column]] <- as.character(values)
      }
    }
    output[[name]] <- table
  }
  output
}

# Figures 5.3/5.4 follow the agreed direct observed-cell rule. This deliberately
# does not add complementary or cross-time suppression to observed occupancy.
protect_observed_counts <- function(occupancy) {
  groups <- c("Eligible", "Not_eligible")
  states <- c("AliveNoPneum", "AlivePneum", "Dead", "ObservationEnded")
  columns <- c("group", "cycle", "day", "state", "patients", "baseline_patients", "proportion")
  occupancy <- occupancy[occupancy$group %in% groups, columns]
  if (anyNA(occupancy) || any(!occupancy$state %in% states) || anyDuplicated(occupancy[c("group", "cycle", "state")]) ||
      any(occupancy$patients < 0 | occupancy$patients != floor(occupancy$patients)) ||
      any(occupancy$baseline_patients <= 0 | occupancy$baseline_patients != floor(occupancy$baseline_patients)) ||
      nrow(occupancy) != length(unique(occupancy$cycle)) * length(groups) * length(states)) {
    stop("Observed occupancy needs four valid state counts per group and weekly boundary.")
  }
  for (group in groups) {
    group_rows <- occupancy[occupancy$group == group, ]
    baseline_count <- unique(group_rows$baseline_patients)
    if (length(baseline_count) != 1L || any(tapply(group_rows$patients, group_rows$cycle, sum) != baseline_count)) stop("Observed occupancy does not sum to baseline.")
  }
  if (any(abs(occupancy$proportion - occupancy$patients / occupancy$baseline_patients) > 1e-10)) stop("Observed counts and proportions disagree.")
  small <- occupancy$patients > 0 & occupancy$patients < 5
  small_denominator <- occupancy$baseline_patients > 0 & occupancy$baseline_patients < 5
  occupancy$disclosure <- ifelse(small, "Count <5; percentage withheld", "Released")
  occupancy$disclosure[small_denominator & !small] <- "Denominator <5; percentage withheld"
  occupancy$patients <- number_text(occupancy$patients)
  occupancy$patients[small] <- "<5"
  occupancy$baseline_patients <- number_text(occupancy$baseline_patients)
  occupancy$baseline_patients[small_denominator] <- "<5"
  occupancy$proportion <- number_text(occupancy$proportion)
  occupancy$proportion[small | small_denominator] <- NA_character_
  occupancy[order(match(occupancy$group, groups), occupancy$cycle, match(occupancy$state, states)), ]
}

# Every figure coordinate comes from the protected tables above.
figure_data <- function(public) {
  take <- function(name, columns) public[[name]][columns]
  percent_text <- function(values, missing = "NA") {
    values <- as.character(values)
    numeric_values <- suppressWarnings(as.numeric(values))
    percent <- rep(missing, length(values))
    known <- is.finite(numeric_values)
    percent[known] <- format(100 * numeric_values[known], digits = 12, trim = TRUE, scientific = FALSE)
    percent[!is.na(values) & values == "SUPPRESSED"] <- "SUPPRESSED"
    percent
  }
  labels <- c(AliveNoPneum = "Alive, no recorded pneumonia",
    AlivePneum = "Alive, pneumonia history", Dead = "Dead", ObservationEnded = "Observation ended")
  occupancy <- take("observed_occupancy", c("group", "cycle", "day", "state", "patients", "baseline_patients", "proportion"))
  occupancy$state_label <- unname(labels[occupancy$state])
  occupancy$percent <- percent_text(occupancy$proportion, NA_character_)
  occupancy$disclosure <- public$observed_occupancy$disclosure
  observed <- occupancy[occupancy$group == "Eligible", ]
  observed$model <- "Observed"
  names(observed)[names(observed) == "patients"] <- "observed_or_expected_patients"
  fitted <- take("fitted_occupancy", c("model", "cycle", "day", "state", "expected_patients", "proportion"))
  names(fitted)[names(fitted) == "expected_patients"] <- "observed_or_expected_patients"
  fitted$state_label <- unname(labels[fitted$state])
  fitted$percent <- percent_text(fitted$proportion, NA_character_)
  fitted$disclosure <- ifelse(fitted$proportion == "SUPPRESSED", "Protected model estimate", "Modelled estimate")
  comparison_columns <- c("model", "cycle", "day", "state", "state_label", "observed_or_expected_patients", "proportion", "percent", "disclosure")

  transitions <- take("transition_estimates", c("group", "period", "stage", "events", "non_events", "denominator", "source_patients", "probability"))
  periods <- c("Week 0", "Weeks 1-4", "Weeks 5-12", "Weeks 13-26", "Weeks 27-52", "Year 2", "Years 3-5", "Years 6-10")
  order <- match(transitions$period, periods)
  if (anyNA(order)) stop("Unknown empirical transition period.")
  bounds <- c(0L, 7L, 35L, 91L, 189L, 371L, 735L, 1827L, 3647L)
  transitions$period_order <- as.character(order)
  transitions$start_day <- as.character(head(bounds, -1L)[order])
  transitions$end_day_exclusive <- as.character(tail(bounds, -1L)[order])
  transitions$percent <- percent_text(transitions$probability)

  mortality <- public$mortality_schedule
  mortality_columns <- c(background = "background_probability", no_pneumonia = "p_death_no_pneumonia", pneumonia = "p_death_pneumonia")
  mortality_long <- do.call(rbind, lapply(names(mortality_columns), function(series) {
    column <- mortality_columns[[series]]
    source <- if (series == "background") rep("External life table", nrow(mortality)) else mortality[[paste0("source_period_", series)]]
    data.frame(series = series, cycle = mortality$cycle, time_years = mortality$time_years,
      period = mortality$period, source_period = source, excess_multiplier = mortality$excess_multiplier,
      weekly_probability = mortality[[column]], annual_equivalent_probability = mortality[[paste0(column, "_annual_equivalent")]],
      annual_equivalent_percent = percent_text(mortality[[paste0(column, "_annual_equivalent")]]),
      stringsAsFactors = FALSE)
  }))
  markers <- data.frame(series = character(), marker = character(), cycle = character(),
    time_years = character(), annual_equivalent_probability = character(), annual_equivalent_percent = character())
  for (series in names(mortality_columns)) {
    series_rows <- mortality_long[mortality_long$series == series, ]
    time_years <- suppressWarnings(as.numeric(series_rows$time_years))
    death_probability <- suppressWarnings(as.numeric(series_rows$weekly_probability))
    background_probability <- suppressWarnings(as.numeric(mortality$background_probability))
    excess_multiplier <- suppressWarnings(as.numeric(series_rows$excess_multiplier))
    if (!all(is.finite(time_years)) || !all(is.finite(death_probability)) || !all(is.finite(background_probability)) || !all(is.finite(excess_multiplier))) next
    add <- function(i, label) {
      if (!length(i)) return(NULL)
      data.frame(series = series, marker = label, cycle = series_rows$cycle[i], time_years = series_rows$time_years[i],
        annual_equivalent_probability = series_rows$annual_equivalent_probability[i],
        annual_equivalent_percent = series_rows$annual_equivalent_percent[i])
    }
    changes <- which(c(FALSE, series_rows$source_period[-1L] != head(series_rows$source_period, -1L)))
    markers <- rbind(markers, add(changes, "Source period changes"))
    if (series != "background") {
      fading <- which(excess_multiplier < 1)
      if (length(fading)) markers <- rbind(markers, add(max(1L, fading[1L] - 1L), "End of direct local support"))
      faded <- which(excess_multiplier == 0)
      if (length(faded)) markers <- rbind(markers, add(faded[1L], "End of assumed excess-hazard decline"))
      non_background <- which(abs(death_probability - background_probability) > 1e-11)
      convergence <- if (!length(non_background)) 1L else max(non_background) + 1L
      if (convergence <= nrow(series_rows)) markers <- rbind(markers, add(convergence, "Converges to background"))
    }
  }
  horizons <- take("economic_horizons", c("reduction", "horizon_years", "eligible_n", "effect_cycles",
    "discount_rate", "incremental_cost", "incremental_qalys", "icer_per_qaly", "result_classification"))
  horizons$reduction_percent <- percent_text(horizons$reduction)
  icer <- suppressWarnings(as.numeric(horizons$icer_per_qaly))
  horizons$positive_icer_for_log_axis <- ifelse(horizons$icer_per_qaly == "SUPPRESSED", "SUPPRESSED",
    ifelse(is.finite(icer) & icer > 0, "TRUE", "FALSE"))
  list(figure_5_3_data = occupancy,
    figure_5_4_data = rbind(observed[comparison_columns], fitted[comparison_columns]),
    figure_5_6_data = transitions, figure_C_1_data = mortality_long,
    figure_C_1_markers = markers, figure_C_2_data = horizons)
}

export_results <- function(tables, config, output_dir) {
  keep <- c("demographics", "observed_population", "transition_estimates",
    "observed_occupancy", "fitted_occupancy", "survival_validation", "weekly_charges",
    "nihss_linkage", "nihss_summary", "cox_covariate_concepts", "cox_estimates",
    "cox_model_support", "cox_model_summary", "cox_ph_tests", "cox_band_comparison",
    "cox_model_diagnostics", "cox_adjustment_variables", "landmark_flow",
    "landmark_cox_estimates", "landmark_cox_support", "landmark_cox_ph",
    "landmark_cox_diagnostics", "landmark_cox_status", "economic_primary",
    "economic_horizons", "economic_duration", "economic_implementation",
    "economic_downstream", "economic_record_variants", "economic_variant_status",
    "economic_resource_inputs", "economic_simulation_check", "economic_simulation_icer",
    "eligibility_definitions", "mortality_schedule", "discussion_comparisons")
  missing <- setdiff(c(keep, "mortality_calibration"), names(tables))
  if (length(missing)) stop("Required result missing: ", paste(missing, collapse = ", "))
  tables$nihss_linkage <- tables$nihss_linkage[
    tables$nihss_linkage$metric != "exact_duplicate_rows", , drop = FALSE]
  calibration <- tables$mortality_calibration
  if (!is.data.frame(calibration) || !all(c("cycles", "deaths", "used_for_projection") %in% names(calibration))) {
    stop("Mortality calibration support is required before export.")
  }
  calibration <- calibration[calibration$used_for_projection %in% TRUE, ]
  tables <- tables[keep]
  rules <- study_release_rules(tables)
  tables$transition_estimates <- tables$transition_estimates[
    tables$transition_estimates$group == "Eligible" &
      tables$transition_estimates$period != tail(config$period_labels, 1L), ]
  mortality_support <- c(calibration$cycles, calibration$deaths, calibration$cycles - calibration$deaths)
  rules$transition_estimates$extra <- mortality_support
  # Figure C.1 has mortality coordinates only; repeated charge inputs stay out.
  tables$mortality_schedule <- tables$mortality_schedule[setdiff(names(tables$mortality_schedule),
    c("cost_no_pneumonia", "cost_pneumonia", "p_onset", "externally_surviving"))]
  rules$mortality_schedule$extra <- mortality_support
  public <- suppress_small_cells(tables, rules)
  public$observed_occupancy <- protect_observed_counts(tables$observed_occupancy)
  public[] <- lapply(public, function(x) { x[] <- lapply(x, number_text); x })
  figures <- figure_data(public)
  main <- public[setdiff(names(public), c("observed_occupancy", "fitted_occupancy", "mortality_schedule", "economic_horizons"))]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) stop("Choose an empty results directory.")
  for (section in c("tables", "figures")) {
    dir.create(file.path(output_dir, section))
    contents <- if (section == "tables") main else figures
    for (name in names(contents)) utils::write.csv(contents[[name]],
      file.path(output_dir, section, paste0(name, ".csv")), row.names = FALSE,
      na = if (name %in% c("figure_5_3_data", "figure_5_4_data")) "" else "NA")
  }
  saveRDS(c(main, figures), file.path(output_dir, "thesis_results.rds"))
  writeLines(c(
    "Thesis results: share this directory only.",
    "tables/ contains manuscript and appendix results, model status and diagnostics.",
    "figures/ contains numeric inputs for Figures 5.3, 5.4, 5.6, C.1 and C.2; C.1 markers use the same weekly coordinates. No images are generated.",
    "thesis_results.rds contains the same flat tables for reading in R.",
    "Counts 1-4 are <5; SUPPRESSED marks related values withheld for disclosure protection. Zero is retained. Do not reconstruct or interpolate masked values.",
    "Observed Figures 5.3/5.4 use direct cell masking: counts 1-4 and their proportions/percentages are hidden. Totals can permit subtraction; no complementary or cross-time suppression is applied to these observed curves.",
    "Observed occupancy retains ObservationEnded in the original denominator. Fitted models are uncensored. Plot only the provided coordinates; do not smooth mortality jumps.",
    "Mortality annual equivalents are model inputs, not observed annual death rates. The one-year economic result covers 52 undiscounted weeks; longer horizons use the stated discount rate.",
    sprintf("Configured assumptions: utility %s; annual discount rate %s; intervention EUR %s; charge price year %s.",
      config$economic$utility, config$economic$discount_rate, config$economic$intervention_cost, config$cost$price_year),
    "FIT means the model/test was computed; inspect its warnings and proportional-hazards tests. Missing or non-estimable results remain missing. Simulation path intervals are not parameter confidence intervals.",
    "Day 8 analyses retain the fixed dysphagia group, restricted to patients alive, observed and without pneumonia through Day 7. CohortContrast and published multi-database results are outside this local rerun."
  ), file.path(output_dir, "README.txt"))
  invisible(list(export_dir = normalizePath(output_dir), tables = c(main, figures)))
}
