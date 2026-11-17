# Economic outcomes and sensitivities from fitted source-cycle schedules.

intervention_cost_scenarios <- function(config) {
  index <- utils::read.csv(config$cost$price_index_file, stringsAsFactors = FALSE)
  if (!all(c("year", "index") %in% names(index))) {
    stop("The supplied CPI file needs year and index columns.", call. = FALSE)
  }
  get_index <- function(year) {
    value <- index$index[index$year == year]
    if (length(value) != 1L || !is.finite(value) || value <= 0) stop("Missing CPI for ", year)
    value
  }
  equipment_factor <- get_index(config$cost$price_year) / get_index(2023L)
  staff_factor <- get_index(config$cost$price_year) / get_index(2025L)
  catheter <- 1000 / 0.85740 * equipment_factor
  session <- 110 / 1.4690 * staff_factor
  station_per_patient <- 10000 / 0.85740 * equipment_factor / 25
  data.frame(
    cost_scenario = c("catheter_only", "three_sessions", "six_sessions", "base_station_25"),
    catheter_cost = catheter,
    staff_sessions = c(0L, 3L, 6L, 3L),
    staff_cost = session * c(0L, 3L, 6L, 3L),
    base_station_allocation = c(0, 0, 0, station_per_patient),
    cost_per_patient = round(catheter + session * c(0L, 3L, 6L, 3L) +
                               c(0, 0, 0, station_per_patient)),
    price_year = config$cost$price_year,
    source = "Archived thesis resource and exchange-rate assumptions with the bundled CPI",
    stringsAsFactors = FALSE
  )
}

classify_economic_result <- function(cost, qaly) {
  if (abs(qaly) < 1e-10) {
    if (cost > 0) return("equal_qalys_higher_cost")
    if (cost < 0) return("equal_qalys_lower_cost")
    return("equal_cost_and_qalys")
  }
  if (qaly > 0) {
    if (cost <= 0) return("intervention_dominant")
    return("more_effective_more_costly")
  }
  if (cost >= 0) return("intervention_dominated")
  "less_effective_less_costly"
}

summarise_economic_outcomes <- function(projection, config, horizon_years,
                                  cost_per_patient = config$economic$intervention_cost) {
  requested_cycles <- ceiling(horizon_years * config$economic$days_per_year / config$economic$cycle_days)
  discount_rate <- config$economic$discount_rate
  if (horizon_years == 1) {
    requested_cycles <- 52L
    discount_rate <- 0
  }
  if (requested_cycles > projection$n_cycles) {
    stop("Requested horizon exceeds the fitted schedule.")
  }
  eligible_count <- projection$eligible_n
  program_cost <- eligible_count * cost_per_patient
  summaries <- arm_totals <- list()
  for (arm_name in c("standard", "intervention")) {
    arm <- projection$arms[[arm_name]]
    rewards <- arm$rewards
    rewards$costDiscountRate <- rewards$healthDiscountRate <- discount_rate
    rewards$upfrontCosts <- c(program = if (arm_name == "intervention") cost_per_patient else 0)
    summary <- TrajectoryMarkovAnalysis::summariseOutcomes(arm$projection, rewards,
      nCycles = min(requested_cycles, nrow(arm$by_cycle)))
    summaries[[arm_name]] <- summary
    totals <- stats::setNames(summary$totals$total, summary$totals$metric)
    cycles <- arm$by_cycle[arm$by_cycle$cycle < requested_cycles, , drop = FALSE]
    care_cost <- unname(totals["discounted_cost"])
    if (arm_name == "intervention") care_cost <- care_cost - program_cost
    arm_totals[[arm_name]] <- list(
      cost = care_cost, life_years = unname(totals["discounted_life_years"]),
      qalys = unname(totals["discounted_qalys"]),
      pneumonia = sum(cycles$pneumonia_onsets), deaths = sum(cycles$deaths),
      remaining_alive = tail(cycles$alive_end, 1L), cycles_run = nrow(cycles))
  }
  standard <- arm_totals$standard
  intervention <- arm_totals$intervention
  delta_cost <- intervention$cost + program_cost - standard$cost
  delta_qaly <- intervention$qalys - standard$qalys
  delta_ly <- intervention$life_years - standard$life_years
  if (abs(delta_qaly) < 1e-10) delta_qaly <- 0
  if (abs(delta_ly) < 1e-10) delta_ly <- 0

  # TMA compares equal-length summaries. Early stopping can give the two arms
  # different lengths; then compare their accumulated totals directly.
  comparison <- NULL
  if (summaries$standard$metadata$nCycles == summaries$intervention$metadata$nCycles) {
    comparison <- TrajectoryMarkovAnalysis::compareScenarios(summaries$standard, summaries$intervention)
  }
  icer_ly <- icer_qaly <- NA_real_
  if (delta_ly > 1e-10) {
    icer_ly <- delta_cost / delta_ly
    if (!is.null(comparison)) icer_ly <- comparison$costEffectiveness$icer[
      comparison$costEffectiveness$outcome == "life_years"]
  }
  if (delta_qaly > 1e-10) {
    icer_qaly <- delta_cost / delta_qaly
    if (!is.null(comparison)) icer_qaly <- comparison$costEffectiveness$icer[
      comparison$costEffectiveness$outcome == "qalys"]
  }
  data.frame(
    eligible_n = eligible_count, horizon_years = horizon_years, requested_cycles = requested_cycles,
    standard_cycles_run = standard$cycles_run, intervention_cycles_run = intervention$cycles_run,
    discount_rate = discount_rate, cost_per_patient = cost_per_patient, program_cost = program_cost,
    standard_cost = standard$cost, intervention_healthcare_cost = intervention$cost,
    intervention_cost = intervention$cost + program_cost, incremental_cost = delta_cost,
    healthcare_cost_change = intervention$cost - standard$cost,
    standard_life_years = standard$life_years, intervention_life_years = intervention$life_years,
    incremental_life_years = delta_ly,
    standard_qalys = standard$qalys, intervention_qalys = intervention$qalys,
    incremental_qalys = delta_qaly, icer_per_life_year = icer_ly, icer_per_qaly = icer_qaly,
    standard_pneumonia = standard$pneumonia, intervention_pneumonia = intervention$pneumonia,
    pneumonia_cases_avoided = standard$pneumonia - intervention$pneumonia,
    standard_deaths = standard$deaths, intervention_deaths = intervention$deaths,
    deaths_avoided = standard$deaths - intervention$deaths,
    standard_remaining_alive = standard$remaining_alive,
    intervention_remaining_alive = intervention$remaining_alive,
    result_classification = classify_economic_result(delta_cost, delta_qaly),
    .support_n = eligible_count, stringsAsFactors = FALSE)
}

# The latent prevented state preserves the pneumonia-associated mortality and
# charge schedule to the extent specified by the reassignment fraction. It can
# still enter actual pneumonia history later. At fraction one, aggregation of
# this state with AliveNoPneum equals the public TMA primary projection.
project_downstream_scenario <- function(schedule, eligible_n, config, fraction) {
  if (length(fraction) != 1L || !is.finite(fraction) || fraction < 0 || fraction > 1) {
    stop("Reassignment fraction must be between zero and one.")
  }
  economic <- config$economic
  states <- c("no_pneumonia", "pneumonia", "prevented", "Dead")
  periods <- paste0("source_cycle_", schedule$cycle)
  breaks <- c(schedule$cycle * economic$cycle_days, Inf)
  identity <- diag(4)
  dimnames(identity) <- list(states, states)

  # Costs and utilities are the same in both arms. Only pneumonia onset differs.
  prevented_cost <- fraction * schedule$cost_no_pneumonia +
    (1 - fraction) * schedule$cost_pneumonia
  costs <- data.frame(period = rep(periods, each = 4), state = rep(states, nrow(schedule)),
    cost = as.vector(t(cbind(schedule$cost_no_pneumonia, schedule$cost_pneumonia, prevented_cost, 0))))
  rewards <- TrajectoryMarkovAnalysis::createRewardModel(
    stateCosts = costs, upfrontCosts = c(program = 0), aliveStates = states[1:3],
    utilities = stats::setNames(c(rep(economic$utility, 3), 0), states),
    costTiming = "onset", healthTiming = "start", discountTiming = "midpoint",
    costDiscountRate = economic$discount_rate, healthDiscountRate = economic$discount_rate,
    currency = "EUR", priceYear = config$cost$price_year, daysPerYear = economic$days_per_year,
    periodBreaks = breaks, periodLabels = periods)

  arms <- list()
  for (arm_name in c("standard", "intervention")) {
    reduction <- if (arm_name == "intervention") economic$reduction else 0
    onset <- mortality <- stats::setNames(vector("list", nrow(schedule)), periods)
    for (cycle_number in seq_len(nrow(schedule))) {
      cycle <- schedule[cycle_number, ]
      active_reduction <- reduction * (cycle$cycle < economic$effect_cycles)
      onset_matrix <- mortality_matrix <- identity
      onset_matrix["no_pneumonia", ] <- c(1 - cycle$p_onset,
        cycle$p_onset * (1 - active_reduction), cycle$p_onset * active_reduction, 0)
      onset_matrix["prevented", ] <- c(0, cycle$p_onset * (1 - active_reduction),
        1 - cycle$p_onset * (1 - active_reduction), 0)
      death_probability <- c(cycle$p_death_no_pneumonia, cycle$p_death_pneumonia,
        fraction * cycle$p_death_no_pneumonia + (1 - fraction) * cycle$p_death_pneumonia)
      for (living_state in 1:3) {
        mortality_matrix[living_state, living_state] <- 1 - death_probability[living_state]
        mortality_matrix[living_state, 4L] <- death_probability[living_state]
      }
      onset[[cycle_number]] <- onset_matrix
      mortality[[cycle_number]] <- mortality_matrix
    }
    model <- TrajectoryMarkovAnalysis::composeTransitionStages(
      list(onset = onset, mortality = mortality), cycleLength = economic$cycle_days,
      absorbingStates = "Dead", periodBreaks = breaks, periodLabels = periods)
    projected <- TrajectoryMarkovAnalysis::projectCohort(model, nCycles = nrow(schedule),
      initial = c(no_pneumonia = 1), cohortSize = eligible_n)
    flows <- projected$stageFlows
    pneumonia <- flows$stage == "onset" & flows$from_state %in% c("no_pneumonia", "prevented") &
      flows$to_state == "pneumonia"
    deaths <- flows$stage == "mortality" & flows$from_state != "Dead" & flows$to_state == "Dead"
    event_counts <- list(pneumonia_onsets = numeric(nrow(schedule)), deaths = numeric(nrow(schedule)))
    for (event_name in names(event_counts)) {
      event_flows <- flows[if (event_name == "deaths") deaths else pneumonia, ]
      if (nrow(event_flows)) {
        counts <- rowsum(event_flows$expected_count, event_flows$source_cycle, reorder = FALSE)
        event_counts[[event_name]][as.integer(rownames(counts)) + 1L] <- counts[, 1L]
      }
    }
    living <- projected$occupancy[projected$occupancy$state != "Dead", ]
    alive <- rowsum(living$expected_count, living$cycle, reorder = FALSE)[
      as.character(0:nrow(schedule)), 1L]
    by_cycle <- data.frame(cycle = schedule$cycle, alive_end = unname(tail(alive, -1)),
      pneumonia_onsets = event_counts$pneumonia_onsets, deaths = event_counts$deaths)
    if (!is.null(economic$stop_alive_threshold)) {
      stopping_cycle <- which(by_cycle$alive_end <= economic$stop_alive_threshold)
      if (length(stopping_cycle)) by_cycle <- by_cycle[seq_len(stopping_cycle[1L]), ]
    }
    arms[[arm_name]] <- list(by_cycle = by_cycle, projection = projected, rewards = rewards)
  }
  list(arms = arms, eligible_n = eligible_n, n_cycles = nrow(schedule))
}

# Simulate standard care and intervention together. Counting patients in each
# pair of states gives the same distribution as drawing every patient separately.
# States 1, 2 and 3 mean no pneumonia, pneumonia history and dead.
# These runs describe path variation with fixed model inputs.
draw_paired_events <- function(counts, p_standard, p_intervention, destinations) {
  next_counts <- matrix(0, nrow(counts), 9L)
  # Nine columns represent every standard/intervention state pair. Each row is
  # one simulation. Patients in the same pair share the same event probabilities.
  for (standard_state in 1:3) {
    for (intervention_state in 1:3) {
      current_pair <- (standard_state - 1L) * 3L + intervention_state
      patients <- counts[, current_pair]
      if (!any(patients > 0)) next
      smaller_probability <- min(p_standard[standard_state], p_intervention[intervention_state])
      larger_probability <- max(p_standard[standard_state], p_intervention[intervention_state])

      # Shared uniforms imply: both events below the smaller probability,
      # one event between the two probabilities, neither above the larger one.
      both <- numeric(length(patients))
      if (smaller_probability == 1) {
        both <- patients
      } else if (smaller_probability > 0) {
        both <- stats::rbinom(length(patients), patients, smaller_probability)
      }
      remaining <- patients - both
      one <- numeric(length(patients))
      if (larger_probability != smaller_probability && smaller_probability != 1) {
        if (larger_probability == 1) {
          one <- remaining
        } else {
          conditional_probability <- (larger_probability - smaller_probability) / (1 - smaller_probability)
          one <- stats::rbinom(length(patients), remaining, conditional_probability)
        }
      }
      neither <- remaining - one
      both_destination <- (destinations[standard_state] - 1L) * 3L + destinations[intervention_state]
      if (p_standard[standard_state] > p_intervention[intervention_state]) {
        one_destination <- (destinations[standard_state] - 1L) * 3L + intervention_state
      } else {
        one_destination <- (standard_state - 1L) * 3L + destinations[intervention_state]
      }
      next_counts[, both_destination] <- next_counts[, both_destination] + both
      next_counts[, one_destination] <- next_counts[, one_destination] + one
      next_counts[, current_pair] <- next_counts[, current_pair] + neither
    }
  }
  next_counts
}

simulate_economic_outcomes <- function(schedule, eligible_n, config, projection, deterministic,
                                     n_simulations = 10000L, seed = 120917L) {
  if (length(n_simulations) != 1L || !is.finite(n_simulations) ||
      n_simulations < 2 || n_simulations != floor(n_simulations)) {
    stop("n_simulations must be an integer of at least two.")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  economic <- config$economic
  paired_counts <- matrix(0, n_simulations, 9L)
  paired_counts[, 1L] <- eligible_n
  standard_state <- rep(1:3, each = 3L)
  intervention_state <- rep(1:3, times = 3L)
  totals <- matrix(0, n_simulations, 4L,
    dimnames = list(NULL, c("standard_cost", "intervention_cost", "standard_qalys", "intervention_qalys")))
  standard_end <- nrow(projection$arms$standard$by_cycle)
  intervention_end <- nrow(projection$arms$intervention$by_cycle)
  for (cycle_number in seq_len(max(standard_end, intervention_end))) {
    cycle <- schedule[cycle_number, ]
    standard_active <- cycle_number <= standard_end
    intervention_active <- cycle_number <= intervention_end
    discount <- (1 + economic$discount_rate)^(
      -(cycle$cycle + 0.5) * economic$cycle_days / economic$days_per_year)
    qaly_factor <- economic$utility * economic$cycle_days / economic$days_per_year * discount

    # Health is accrued by patients alive at the start of the cycle.
    if (standard_active) {
      totals[, "standard_qalys"] <- totals[, "standard_qalys"] +
        as.vector(paired_counts %*% as.numeric(standard_state != 3L)) * qaly_factor
    }
    if (intervention_active) {
      totals[, "intervention_qalys"] <- totals[, "intervention_qalys"] +
        as.vector(paired_counts %*% as.numeric(intervention_state != 3L)) * qaly_factor
    }

    # Pneumonia occurs first; healthcare charges use the resulting state.
    standard_onset <- c(if (standard_active) cycle$p_onset else 0, 0, 0)
    intervention_probability <- 0
    if (intervention_active) {
      intervention_probability <- cycle$p_onset *
        (1 - economic$reduction * (cycle$cycle < economic$effect_cycles))
    }
    intervention_onset <- c(intervention_probability, 0, 0)
    paired_counts <- draw_paired_events(paired_counts, standard_onset,
      intervention_onset, c(2L, 2L, 3L))
    rewards <- c(cycle$cost_no_pneumonia, cycle$cost_pneumonia, 0)
    if (standard_active) {
      totals[, "standard_cost"] <- totals[, "standard_cost"] +
        as.vector(paired_counts %*% rewards[standard_state]) * discount
    }
    if (intervention_active) {
      totals[, "intervention_cost"] <- totals[, "intervention_cost"] +
        as.vector(paired_counts %*% rewards[intervention_state]) * discount
    }

    # Mortality ends the cycle. Stopped arms receive no further events or rewards.
    death_probability <- c(cycle$p_death_no_pneumonia, cycle$p_death_pneumonia, 0)
    paired_counts <- draw_paired_events(paired_counts,
      death_probability * standard_active, death_probability * intervention_active, c(3L, 3L, 3L))
  }
  totals[, "intervention_cost"] <- totals[, "intervention_cost"] + eligible_n * economic$intervention_cost
  draws <- data.frame(totals, incremental_cost = totals[, "intervention_cost"] - totals[, "standard_cost"],
    incremental_qalys = totals[, "intervention_qalys"] - totals[, "standard_qalys"])
  summary_rows <- list()
  for (metric in names(draws)) {
    values <- draws[[metric]]
    target <- deterministic[[metric]][[1L]]
    simulated_mean <- mean(values)
    simulated_sd <- stats::sd(values)
    mcse <- simulated_sd / sqrt(n_simulations)
    difference <- simulated_mean - target
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      metric = metric, deterministic = target, simulated_mean = simulated_mean,
      simulated_sd = simulated_sd, mcse = mcse,
      central_95_lower = unname(stats::quantile(values, .025)),
      central_95_upper = unname(stats::quantile(values, .975)),
      mean_minus_deterministic = difference,
      reconciliation = if (abs(difference) <= max(4 * mcse, 1e-8)) "PASS" else "REVIEW",
      simulations = n_simulations, seed = seed, .support_n = eligible_n)
  }
  summaries <- do.call(rbind, summary_rows)
  list(summary = summaries,
    ratio_of_means = data.frame(
      deterministic_icer = deterministic$icer_per_qaly,
      simulated_icer_from_mean_increments = if (mean(draws$incremental_qalys) > 0)
        mean(draws$incremental_cost) / mean(draws$incremental_qalys) else NA_real_,
      simulations = n_simulations, .support_n = eligible_n))
}

analyse_cost_effectiveness <- function(primary_result, config, variants = list(),
                                   n_simulations = 10000L, seed = 120917L) {
  reductions <- c(.10, .20, .25, .40)
  horizons <- c(1, 5, 10, 20, 40)
  eligible_count <- primary_result$population$eligible_patients[[1L]]
  schedule <- primary_result$schedule
  expected_cycles <- ceiling(40 * config$economic$days_per_year / config$economic$cycle_days)
  if (nrow(schedule) < expected_cycles) stop("Full thesis outputs require a schedule through 40 years.")
  implementation <- intervention_cost_scenarios(config)
  primary_cost <- implementation$cost_per_patient[implementation$cost_scenario == "three_sessions"]
  if (primary_cost != config$economic$intervention_cost) {
    stop("Configured primary intervention price differs from the archived resource assumptions.")
  }
  central_config <- config
  central_config$economic$reduction <- .25
  central_config$economic$effect_cycles <- 4L

  # Main effect assumptions and follow-up horizons.
  horizon_rows <- list()
  for (reduction in reductions) {
    scenario_config <- central_config
    scenario_config$economic$reduction <- reduction
    projection <- project_standard_and_intervention(schedule, eligible_count, scenario_config)
    if (reduction == .25) central_projection <- projection
    rows <- list()
    for (horizon_number in seq_along(horizons)) {
      outcomes <- summarise_economic_outcomes(projection, config, horizons[horizon_number])
      rows[[horizon_number]] <- cbind(reduction = reduction, effect_cycles = 4L, outcomes)
    }
    horizon_rows[[as.character(reduction)]] <- do.call(rbind, rows)
  }
  horizon_results <- do.call(rbind, horizon_rows)
  central_outcomes <- summarise_economic_outcomes(central_projection, config, 40)

  # Effect duration and intervention resource assumptions.
  durations <- c(1L, 4L, 12L, ceiling(10 * config$economic$days_per_year / config$economic$cycle_days))
  duration_rows <- list()
  for (duration_number in seq_along(durations)) {
    effect_cycles <- durations[duration_number]
    scenario_config <- central_config
    scenario_config$economic$effect_cycles <- effect_cycles
    projection <- central_projection
    if (effect_cycles != 4L) projection <- project_standard_and_intervention(schedule, eligible_count, scenario_config)
    duration_rows[[duration_number]] <- cbind(reduction = .25, effect_cycles = effect_cycles,
      summarise_economic_outcomes(projection, config, 40))
  }
  cost_rows <- list()
  for (scenario_number in seq_len(nrow(implementation))) {
    cost_rows[[scenario_number]] <- cbind(cost_scenario = implementation$cost_scenario[scenario_number],
      summarise_economic_outcomes(central_projection, config, 40, implementation$cost_per_patient[scenario_number]))
  }

  # How much of pneumonia-associated mortality and cost follows prevented cases?
  downstream_rows <- list()
  for (fraction in c(0, .5, 1)) {
    projection <- project_downstream_scenario(schedule, eligible_count, central_config, fraction)
    outcomes <- summarise_economic_outcomes(projection, config, 40)
    if (fraction == 1) {
      columns <- c("incremental_cost", "incremental_qalys", "pneumonia_cases_avoided")
      if (!isTRUE(all.equal(outcomes[columns], central_outcomes[columns], tolerance = 1e-8))) {
        stop("The expanded downstream model does not reconcile with the TMA primary projection.")
      }
    }
    downstream_rows[[length(downstream_rows) + 1L]] <- cbind(downstream_reassignment = fraction, outcomes)
  }

  # Repeat the economic calculations using the separately re-estimated cohorts.
  variant_rows <- list()
  variant_status <- data.frame(variant = c("strict_chronology", "broad_marker", "day8_landmark"),
    status = "NOT_SUPPLIED", detail = "No variant fit supplied", stringsAsFactors = FALSE)
  for (variant_number in seq_len(nrow(variant_status))) {
    variant_name <- variant_status$variant[variant_number]
    fit <- variants[[variant_name]]
    if (is.null(fit)) next
    variant_status$status[variant_number] <- "COMPLETE"
    variant_status$detail[variant_number] <- if (!is.null(fit$variant_note)) fit$variant_note else "Re-estimated inputs supplied"
    variant_reductions <- if (variant_name == "broad_marker") reductions else .25
    for (reduction in variant_reductions) {
      scenario_config <- central_config
      scenario_config$economic$reduction <- reduction
      projection <- project_standard_and_intervention(fit$schedule, fit$population$eligible_patients[[1L]], scenario_config)
      outcomes <- summarise_economic_outcomes(projection, config, 40)
      variant_rows[[paste(variant_name, reduction)]] <- cbind(variant = variant_name, reduction = reduction,
        definition = variant_status$detail[variant_number], outcomes)
    }
  }
  record_variants <- if (length(variant_rows)) do.call(rbind, variant_rows) else data.frame()
  simulation <- simulate_economic_outcomes(schedule, eligible_count, central_config,
    central_projection, central_outcomes, n_simulations, seed)

  discussion <- record_variants[record_variants$reduction == .25, ]
  discussion$primary_icer <- rep(central_outcomes$icer_per_qaly, nrow(discussion))
  discussion$icer_percent_change <- 100 * (discussion$icer_per_qaly / central_outcomes$icer_per_qaly - 1)
  cpi <- utils::read.csv(config$cost$price_index_file)
  contextual_reference <- 40000 * cpi$index[cpi$year == 2026] / cpi$index[cpi$year == 2019]
  discussion$contextual_reference_eur_per_qaly <- rep(contextual_reference, nrow(discussion))
  discussion$below_contextual_reference <- discussion$icer_per_qaly < discussion$contextual_reference_eur_per_qaly

  list(economic_primary = horizon_results[horizon_results$horizon_years %in% c(1, 40), ],
    economic_horizons = horizon_results,
    economic_duration = do.call(rbind, duration_rows),
    economic_implementation = do.call(rbind, cost_rows),
    economic_downstream = do.call(rbind, downstream_rows),
    economic_record_variants = record_variants, economic_variant_status = variant_status,
    economic_resource_inputs = implementation, economic_simulation_check = simulation$summary,
    economic_simulation_icer = simulation$ratio_of_means, discussion_comparisons = discussion)
}
