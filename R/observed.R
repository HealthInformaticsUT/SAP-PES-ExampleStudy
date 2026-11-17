# Scientific analyses: observed follow-up, exact-date Cox models and economics.
# Report functions return aggregate tables without retaining fitted model objects.

# Observed summaries, period-specific occupancy and fixed-input validation.
cohort_groups <- function(patients) {
  list(All = rep(TRUE, nrow(patients)), Eligible = patients$eligible,
       Not_eligible = !patients$eligible)
}

summarise_demographics <- function(patients) {
  groups <- cohort_groups(patients)
  rows <- list()
  for (group_number in seq_along(groups)) {
    group <- names(groups)[group_number]
    group_patients <- patients[groups[[group_number]], , drop = FALSE]
    patient_count <- nrow(group_patients)
    age <- group_patients$age[is.finite(group_patients$age)]
    index_year <- as.integer(format(group_patients$index_date, "%Y"))
    index_year <- index_year[!is.na(index_year)]
    year_quartiles <- rep(NA_real_, 3)
    if (length(index_year)) {
      year_quartiles <- stats::quantile(index_year, c(.25, .5, .75), names = FALSE)
    }
    sex <- as.character(group_patients$gender_concept_id)
    women <- sum(sex == "8532", na.rm = TRUE)
    men <- sum(sex == "8507", na.rm = TRUE)
    rows[[group_number]] <- data.frame(
      group = group, patients = patient_count,
      age_available = length(age), age_missing = patient_count - length(age),
      age_mean = if (length(age)) mean(age) else NA_real_,
      age_sd = if (length(age) > 1L) stats::sd(age) else NA_real_,
      female = women, male = men,
      sex_other_or_missing = sum(is.na(sex) | !sex %in% c("8507", "8532")),
      female_percent = if (patient_count) 100 * women / patient_count else NA_real_,
      index_year_median = year_quartiles[2L], index_year_q1 = year_quartiles[1L],
      index_year_q3 = year_quartiles[3L], .support_n = patient_count,
      .support_age = length(age), .support_index_year = length(index_year))
  }
  do.call(rbind, rows)
}

summarise_followup <- function(c2t, patients, histories) {
  groups <- cohort_groups(patients)
  source_intervals <- c2t$source_intervals
  if (is.null(source_intervals)) source_intervals <- c2t$cleaned_intervals
  rows <- list()
  for (group_number in seq_along(groups)) {
    group <- names(groups)[group_number]
    group_patients <- patients[groups[[group_number]], , drop = FALSE]
    patient_cycles <- histories$data[
      histories$data$subject_id %in% group_patients$subject_id, , drop = FALSE]
    observation_days <- as.integer(group_patients$observation_end - group_patients$index_date)
    died <- !is.na(group_patients$death_day)
    had_pneumonia <- !is.na(group_patients$pneumonia_day)
    # First-year events include Day 0 through Day 364, inclusive.
    first_year_pneumonia <- sum(had_pneumonia & group_patients$pneumonia_day >= 0L &
      group_patients$pneumonia_day < 365L)
    first_year_deaths <- sum(died & group_patients$death_day >= 0L & group_patients$death_day < 365L)
    patient_count <- nrow(group_patients)
    rows[[group_number]] <- data.frame(
      group = group, patients = patient_count,
      source_state_rows = if (is.null(source_intervals)) NA_integer_ else
        sum(source_intervals$subject_id %in% group_patients$subject_id),
      observed_patient_cycles = nrow(patient_cycles),
      shortened_patient_cycles = sum(patient_cycles$observed_days < 7L),
      pneumonia_full_followup = sum(had_pneumonia), deaths_full_followup = sum(died),
      censored_alive = sum(!died),
      observed_beyond_365_days = sum(observation_days > 365L & (!died | group_patients$death_day > 365L)),
      pneumonia_by_day365 = first_year_pneumonia, deaths_by_day365 = first_year_deaths,
      deaths_by_day365_percent = if (patient_count) 100 * first_year_deaths / patient_count else NA_real_,
      same_day_pneumonia_death = sum(had_pneumonia & died &
        group_patients$pneumonia_day == group_patients$death_day),
      represented_cycle_days = sum(patient_cycles$observed_days),
      observed_alive_days = sum(pmin(observation_days,
        ifelse(died, group_patients$death_day, Inf)) + 1L),
      .support_n = patient_count)
  }
  do.call(rbind, rows)
}

summarise_transitions <- function(cycles, patients, config) {
  groups <- cohort_groups(patients)
  stages <- c("Pneumonia onset", "Mortality without pneumonia", "Mortality after pneumonia")
  rows <- list()
  for (group_number in seq_along(groups)) {
    group <- names(groups)[group_number]
    patient_ids <- patients$subject_id[groups[[group_number]]]
    group_cycles <- cycles[cycles$subject_id %in% patient_ids, , drop = FALSE]
    if (!nrow(group_cycles)) {
      empty_periods <- expand.grid(stage = stages, period = config$period_labels,
        stringsAsFactors = FALSE)
      rows[[group_number]] <- data.frame(group = group, period = empty_periods$period,
        stage = empty_periods$stage, events = 0, non_events = 0, denominator = 0,
        source_patients = 0, probability = NA_real_, .support_n = 0)
      next
    }
    estimates <- estimate_transitions(group_cycles, config)$estimates
    onset <- estimates$stage == "onset" & estimates$from_state == "AliveNoPneum" &
      estimates$to_state == "AlivePneum"
    death_without_pneumonia <- estimates$stage == "mortality" &
      estimates$from_state == "AliveNoPneum" & estimates$to_state == "DeathNoHistory"
    death_after_pneumonia <- estimates$stage == "mortality" &
      estimates$from_state == "AlivePneum" & estimates$to_state == "DeathHistory"
    estimates <- estimates[onset | death_without_pneumonia | death_after_pneumonia, , drop = FALSE]
    stage_label <- rep(stages[3L], nrow(estimates))
    stage_label[estimates$from_state == "AliveNoPneum"] <- stages[2L]
    stage_label[estimates$stage == "onset"] <- stages[1L]
    group_results <- data.frame(
      group = group, period = estimates$period, stage = stage_label,
      events = estimates$count, non_events = estimates$denominator - estimates$count,
      denominator = estimates$denominator, source_patients = estimates$source_patients,
      probability = estimates$probability, .support_n = estimates$source_patients)
    row_order <- order(match(group_results$period, config$period_labels),
      match(group_results$stage, stages))
    rows[[group_number]] <- group_results[row_order, ]
  }
  do.call(rbind, rows)
}

# Weekly boundary occupancy is evaluated just before that boundary. Thus a
# record on Day 0 appears in the distribution at Day 7, after the index cycle.
observed_state_occupancy <- function(patients, max_cycle = 52L) {
  states <- c("AliveNoPneum", "AlivePneum", "Dead", "ObservationEnded")
  groups <- cohort_groups(patients)
  rows <- list()
  for (group in names(groups)) {
    group_patients <- patients[groups[[group]], , drop = FALSE]
    observation_days <- as.integer(group_patients$observation_end - group_patients$index_date)
    for (cycle in 0:max_cycle) {
      day <- cycle * 7L
      state <- ifelse(!is.na(group_patients$pneumonia_day) & group_patients$pneumonia_day < day,
        "AlivePneum", "AliveNoPneum")
      state <- ifelse(observation_days < day, "ObservationEnded", state)
      state <- ifelse(!is.na(group_patients$death_day) & group_patients$death_day < day, "Dead", state)
      counts <- tabulate(match(state, states), length(states))
      rows[[length(rows) + 1L]] <- data.frame(group = group, cycle = cycle,
        day = day, state = states, patients = counts, baseline_patients = nrow(group_patients),
        proportion = if (nrow(group_patients)) counts / nrow(group_patients) else NA_real_, .support_n = nrow(group_patients))
    }
  }
  do.call(rbind, rows)
}

period_transition_probabilities <- function(transitions, config, pooled = FALSE) {
  rates <- transitions[transitions$group == "Eligible" &
    transitions$period != tail(config$period_labels, 1L), , drop = FALSE]
  periods <- head(config$period_labels, -1L)
  if (pooled) {
    # Figure 5.4 pools the five first-year source periods, through Weeks 27-52.
    first_year <- config$period_labels[head(config$period_breaks, -1L) < 365]
    rates <- rates[rates$period %in% first_year, , drop = FALSE]
    rates <- stats::aggregate(cbind(events, denominator) ~ stage, data = rates, FUN = sum)
    rates$period <- "Pooled"
    periods <- "Pooled"
  }
  probabilities <- data.frame(period = periods)
  stage_names <- c(onset = "Pneumonia onset", death_no = "Mortality without pneumonia",
    death_pn = "Mortality after pneumonia")
  for (column in names(stage_names)) {
    stage_rates <- rates[rates$stage == stage_names[[column]], , drop = FALSE]
    probability <- stage_rates$events / stage_rates$denominator
    probabilities[[column]] <- probability[match(periods, stage_rates$period)]
  }
  probabilities
}

project_state_occupancy <- function(transitions, config, patients, max_cycle = 52L) {
  eligible_count <- sum(patients$eligible)
  states <- c("AliveNoPneum", "AlivePneum", "Dead")
  identity <- diag(length(states))
  dimnames(identity) <- list(states, states)
  results <- list()
  for (pooled in c(FALSE, TRUE)) {
    rates <- period_transition_probabilities(transitions, config, pooled)
    period_labels <- config$period_labels
    period_breaks <- config$period_breaks
    model_name <- "Period specific"
    if (pooled) {
      period_labels <- "Pooled"
      period_breaks <- c(0, Inf)
      model_name <- "Pooled"
    }
    onset <- mortality <- list()
    for (period in period_labels) {
      rate <- rates[match(period, rates$period), , drop = FALSE]
      onset_matrix <- mortality_matrix <- identity
      # Unsupported rows stay missing; the projection rejects a missing row
      # only if patients reach it. The extrapolated period is not used here.
      onset_matrix["AliveNoPneum", ] <- NA_real_
      mortality_matrix[c("AliveNoPneum", "AlivePneum"), ] <- NA_real_
      if (is.finite(rate$onset)) {
        onset_matrix["AliveNoPneum", ] <- c(1 - rate$onset, rate$onset, 0)
      }
      if (is.finite(rate$death_no)) {
        mortality_matrix["AliveNoPneum", ] <- c(1 - rate$death_no, 0, rate$death_no)
      }
      if (is.finite(rate$death_pn)) {
        mortality_matrix["AlivePneum", ] <- c(0, 1 - rate$death_pn, rate$death_pn)
      }
      onset[[period]] <- onset_matrix
      mortality[[period]] <- mortality_matrix
    }
    model <- TrajectoryMarkovAnalysis::composeTransitionStages(
      list(onset = onset, mortality = mortality), cycleLength = config$economic$cycle_days,
      absorbingStates = "Dead", periodBreaks = period_breaks, periodLabels = period_labels)
    occupancy <- TrajectoryMarkovAnalysis::projectCohort(model, nCycles = max_cycle,
      initial = c(AliveNoPneum = 1), cohortSize = eligible_count)$occupancy
    results[[length(results) + 1L]] <- data.frame(
      model = model_name, cycle = occupancy$cycle, day = occupancy$elapsed_day,
      state = occupancy$state, expected_patients = occupancy$expected_count,
      proportion = occupancy$probability, .support_n = eligible_count,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, results)
}

# Simulate first exits with geometric waiting times within each fixed-rate
# period. This is equivalent to weekly draws. Each patient keeps their recorded
# observation end, which may coincide with death: this is an internal fit check.
# Transition inputs are fixed; no parameter or causal-effect uncertainty is drawn.
simulate_validation_deaths <- function(patients, stage, config, replicates, horizon_days) {
  patient_count <- nrow(patients)
  patient_number <- rep(seq_len(patient_count), times = replicates)
  censor_week <- pmin(as.integer(patients$observation_end - patients$index_date) + 1L, horizon_days)[patient_number] / 7
  simulated_state <- integer(length(patient_number)) # 0: no pneumonia, 1: pneumonia, 2: dead
  death_week <- rep(NA_integer_, length(patient_number))
  period_starts <- head(config$period_breaks, -2L) / 7
  period_stops <- config$period_breaks[seq.int(2L, length(config$period_breaks) - 1L)] / 7
  waiting_time <- function(number, probability) {
    if (probability <= 0) rep(Inf, number) else stats::rgeom(number, probability) + 1
  }
  for (period_number in seq_along(period_starts)) {
    period_start <- period_starts[[period_number]]
    period_end <- min(period_stops[[period_number]], horizon_days / 7)
    if (period_end <= period_start) next
    active_patients <- which(simulated_state != 2L & censor_week > period_start)
    if (!length(active_patients)) break
    observed_weeks <- pmax(0, pmin(censor_week[active_patients], period_end) - period_start)
    complete_weeks <- floor(observed_weeks + 1e-10)
    partial_week <- observed_weeks - complete_weeks
    probability <- stage[period_number, ]
    if (any(!is.finite(unlist(probability[c("onset", "death_no", "death_pn")]))))
      stop("Validation lacks supported stage probabilities in ", probability$period, ".")
    already_has_pneumonia <- simulated_state[active_patients] == 1L
    # Patients who entered the period with pneumonia.
    within_period <- which(already_has_pneumonia)
    if (length(within_period)) {
      wait <- waiting_time(length(within_period), probability$death_pn)
      died <- wait <= complete_weeks[within_period]
      death_week[active_patients[within_period[died]]] <- period_start + wait[died]
      simulated_state[active_patients[within_period[died]]] <- 2L
    }
    # First exit from no pneumonia, allowing same-cycle pneumonia and death.
    within_period <- which(!already_has_pneumonia)
    if (length(within_period)) {
      exit_probability <- 1 - (1 - probability$onset) * (1 - probability$death_no)
      wait <- waiting_time(length(within_period), exit_probability)
      exits <- which(wait <= complete_weeks[within_period])
      if (length(exits)) {
        exit_positions <- within_period[exits]
        patients_with_event <- active_patients[exit_positions]
        survived_pneumonia <- stats::runif(length(exit_positions)) < probability$onset * (1 - probability$death_pn) / exit_probability
        deaths <- patients_with_event[!survived_pneumonia]
        death_week[deaths] <- period_start + wait[exits][!survived_pneumonia]
        simulated_state[deaths] <- 2L
        pneumonia_survivors <- patients_with_event[survived_pneumonia]
        simulated_state[pneumonia_survivors] <- 1L
        if (length(pneumonia_survivors)) {
          later_death_wait <- waiting_time(length(pneumonia_survivors), probability$death_pn)
          remaining_weeks <- complete_weeks[exit_positions[survived_pneumonia]] - wait[exits][survived_pneumonia]
          died <- later_death_wait <= remaining_weeks
          death_week[pneumonia_survivors[died]] <- period_start + wait[exits][survived_pneumonia][died] + later_death_wait[died]
          simulated_state[pneumonia_survivors[died]] <- 2L
        }
      }
    }
    # For a shortened final cycle, convert fixed weekly probabilities using a
    # constant hazard over its observed fraction. The log-rank groups remain
    # weekly, so death and censoring within that fraction share a weekly bin.
    within_period <- which(partial_week > 1e-10 & simulated_state[active_patients] != 2L)
    if (length(within_period)) {
      patients_with_event <- active_patients[within_period]
      observed_fraction <- partial_week[within_period]
      onset_probability <- -expm1(log1p(-probability$onset) * observed_fraction)
      new_pneumonia <- simulated_state[patients_with_event] == 0L & stats::runif(length(within_period)) < onset_probability
      simulated_state[patients_with_event[new_pneumonia]] <- 1L
      death_probability <- ifelse(simulated_state[patients_with_event] == 1L, probability$death_pn, probability$death_no)
      deaths <- stats::runif(length(within_period)) < -expm1(log1p(-death_probability) * observed_fraction)
      death_week[patients_with_event[deaths]] <- period_start + complete_weeks[within_period[deaths]] + 1L
      simulated_state[patients_with_event[deaths]] <- 2L
    }
  }
  matrix(death_week, nrow = patient_count, ncol = replicates)
}

validate_survival <- function(patients, transitions, config, replicates = 10000L,
                                    seed = 20261111L, batch_size = 250L) {
  if (!is.numeric(replicates) || length(replicates) != 1L || replicates < 1L ||
      replicates != as.integer(replicates)) stop("validation_replicates must be a positive integer.")
  eligible_patients <- patients[patients$eligible, , drop = FALSE]
  empirical_end <- tail(config$period_breaks[is.finite(config$period_breaks)], 1L)
  horizon <- min(empirical_end, max(as.integer(eligible_patients$observation_end - eligible_patients$index_date) + 1L))
  week_count <- ceiling(horizon / 7)
  censor_week <- ceiling(pmin(as.integer(eligible_patients$observation_end - eligible_patients$index_date) + 1L, horizon) / 7)
  death_week <- ifelse(!is.na(eligible_patients$death_day) & eligible_patients$death_day < horizon, floor(eligible_patients$death_day / 7) + 1L, NA_integer_)
  death_week[!is.na(death_week) & death_week > censor_week] <- NA_integer_
  followup_end <- ifelse(is.na(death_week), censor_week, death_week)
  observed_risk <- rev(cumsum(rev(tabulate(followup_end, week_count))))
  observed_death <- tabulate(death_week[!is.na(death_week)], week_count)
  life_table <- data.frame(cycle = seq_len(week_count) - 1L, at_risk = observed_risk,
    deaths = observed_death, .support_n = observed_risk)
  stage <- period_transition_probabilities(transitions, config)
  # Only rates reached by observed windows are needed.
  reached <- which(head(config$period_breaks, -2L) < horizon)
  if (any(!is.finite(as.matrix(stage[reached, c("onset", "death_no", "death_pn")])))) {
    return(list(summary = data.frame(replicates = 0L, valid_comparisons = 0L,
      median_p_value = NA_real_, comparisons_below_0_05 = NA_integer_,
      proportion_below_0_05 = NA_real_, monte_carlo_se = NA_real_,
      horizon_days = horizon, patients = nrow(eligible_patients), observed_deaths = sum(observed_death),
      status = "Not estimated: an empirical period/state lacks support", .support_n = nrow(eligible_patients)),
      life_table = life_table, p_values = numeric()))
  }
  old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old_seed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  pvalues <- numeric(replicates)
  for (first in seq.int(1L, replicates, by = batch_size)) {
    size <- min(batch_size, replicates - first + 1L)
    simulated <- simulate_validation_deaths(eligible_patients, stage, config, size, horizon)
    for (replicate_number in seq_len(size)) {
      simulated_death <- simulated[, replicate_number]
      endpoint <- ifelse(is.na(simulated_death), censor_week, simulated_death)
      simulated_risk <- rev(cumsum(rev(tabulate(endpoint, week_count))))
      simulated_deaths <- tabulate(simulated_death[!is.na(simulated_death)], week_count)
      combined_risk <- observed_risk + simulated_risk
      combined_deaths <- observed_death + simulated_deaths
      # Divide before multiplying so integer patient counts cannot overflow.
      # These are the same expected deaths and hypergeometric variance.
      expected <- ifelse(combined_risk > 0, (observed_risk / combined_risk) * combined_deaths, 0)
      variance <- ifelse(combined_risk > 1, (observed_risk / combined_risk) * (simulated_risk / combined_risk) * combined_deaths *
        ((combined_risk - combined_deaths) / (combined_risk - 1)), 0)
      score <- sum(observed_death - expected)
      total_variance <- sum(variance)
      pvalues[first + replicate_number - 1L] <- if (total_variance > 0) stats::pchisq(score^2 / total_variance, 1L, lower.tail = FALSE) else NA_real_
    }
  }
  valid <- sum(is.finite(pvalues))
  below <- sum(pvalues < .05, na.rm = TRUE)
  rejection_fraction <- if (valid) below / valid else NA_real_
  list(summary = data.frame(replicates = replicates, valid_comparisons = valid,
    median_p_value = if (valid) stats::median(pvalues, na.rm = TRUE) else NA_real_,
    comparisons_below_0_05 = below, proportion_below_0_05 = rejection_fraction,
    monte_carlo_se = if (valid) sqrt(rejection_fraction * (1 - rejection_fraction) / valid) else NA_real_,
    horizon_days = horizon, patients = nrow(eligible_patients), observed_deaths = sum(observed_death),
    status = if (valid) "Estimated" else "No informative death comparisons", .support_n = nrow(eligible_patients)),
    life_table = life_table, p_values = pvalues)
}

summarise_observed_care <- function(c2t, patients, config, primary_result = NULL,
                                  validation_replicates = 10000L, seed = 20261111L) {
  required <- c("subject_id", "index_date", "observation_end", "pneumonia_day", "death_day",
                "age", "gender_concept_id", "eligible")
  if (!all(required %in% names(patients)) || anyNA(patients$eligible) || anyDuplicated(patients$subject_id))
    stop("Observed reports need one valid baseline row per target patient and an eligibility flag.")
  histories <- TrajectoryMarkovAnalysis::prepareTrajectories(c2t$trajectories, cycleLength = 7L,
    states = c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory"),
    absorbingStates = c("DeathNoHistory", "DeathHistory"))
  cycles <- build_weekly_analysis_rows(patients, histories, config)
  transitions <- summarise_transitions(cycles, patients, config)
  validation <- validate_survival(patients, transitions, config, validation_replicates, seed)
  tables <- list(demographics = summarise_demographics(patients),
    observed_population = summarise_followup(c2t, patients, histories),
    transition_estimates = transitions,
    observed_occupancy = observed_state_occupancy(patients),
    fitted_occupancy = project_state_occupancy(transitions, config, patients),
    survival_validation = validation$summary,
    validation_life_table = validation$life_table)
  if (!is.null(primary_result)) {
    charges <- primary_result$costs$rewards
    intervals <- primary_result$costs$fit$intervalCosts
    charges$.support_n <- integer(nrow(charges))
    for (row in seq_len(nrow(charges))) {
      if (isTRUE(charges$extrapolated[row])) {
        contributing_cycles <- intervals[intervals$period == primary_result$costs$last_supported_period, , drop = FALSE]
      } else {
        contributing_cycles <- intervals[intervals$period == charges$period[row] &
          intervals$cost_state == charges$state[row], , drop = FALSE]
      }
      charges$.support_n[row] <- length(unique(contributing_cycles$subject_id))
    }
    tables$weekly_charges <- charges
  }
  tables
}
