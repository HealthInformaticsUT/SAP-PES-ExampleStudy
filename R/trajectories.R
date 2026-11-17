# Weekly stages, costs and deterministic projections for the study.

build_weekly_analysis_rows <- function(patients, histories, config) {
  weekly <- histories$data
  baseline <- patients[match(weekly$subject_id, patients$subject_id), ]
  start_day <- as.integer(weekly$state_start_date - baseline$index_date)
  end_day <- as.integer(weekly$state_end_date - baseline$index_date)
  if (any(start_day %% 7 != 0)) stop("Weekly histories must start on an I63-index weekly boundary.")
  death <- !is.na(baseline$death_day) & baseline$death_day >= start_day & baseline$death_day <= end_day
  prior_pneumonia <- !is.na(baseline$pneumonia_day) & baseline$pneumonia_day < start_day
  onset <- !prior_pneumonia & !is.na(baseline$pneumonia_day) &
    baseline$pneumonia_day >= start_day & baseline$pneumonia_day <= end_day &
    (is.na(baseline$death_day) | baseline$pneumonia_day <= baseline$death_day)
  mortality_state <- ifelse(prior_pneumonia | onset, "AlivePneum", "AliveNoPneum")
  expected <- ifelse(death,
    ifelse(mortality_state == "AlivePneum", "DeathHistory", "DeathNoHistory"), mortality_state)
  if (any(expected != weekly$state_label) || any(!is.na(baseline$death_day) & baseline$death_day < start_day)) {
    stop("Weekly labels disagree with the exact pneumonia/death dates.")
  }
  next_observed <- paste(weekly$subject_id, start_day + 7) %in% paste(weekly$subject_id, start_day)
  data.frame(subject_id = weekly$subject_id, cycle = start_day / 7,
    period = config$period_labels[findInterval(start_day, config$period_breaks)],
    start_state = ifelse(prior_pneumonia, "AlivePneum", "AliveNoPneum"),
    pneumonia_onset = onset, mortality_state = mortality_state,
    death_in_cycle = death, complete_cycle = next_observed | death)
}

# Onset happens before mortality within the same week. A shortened censoring
# week has no observed destination; a death supplies a complete destination.

stage_observations <- function(cycles, stages = c("onset", "mortality")) {
  rows <- lapply(stages, function(stage) {
    source <- if (stage == "onset") cycles$start_state else cycles$mortality_state
    destination <- if (stage == "onset") {
      ifelse(cycles$pneumonia_onset, "AlivePneum", cycles$start_state)
    } else ifelse(cycles$death_in_cycle,
      ifelse(cycles$mortality_state == "AlivePneum", "DeathHistory", "DeathNoHistory"),
      cycles$mortality_state)
    destination[!cycles$complete_cycle] <- NA_character_
    data.frame(subject_id = cycles$subject_id, source_day = cycles$cycle * 7L,
      stage = stage, from_state = source, to_state = destination,
      observed = as.logical(cycles$complete_cycle), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

estimate_transitions <- function(cycles, config, stages = c("onset", "mortality"),
                               observations = NULL) {
  if (is.null(observations)) observations <- stage_observations(cycles, stages)
  TrajectoryMarkovAnalysis::estimateStageTransitions(observations,
    stageOrder = stages, cycleLength = config$economic$cycle_days,
    states = c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory"),
    absorbingStates = c("DeathNoHistory", "DeathHistory"),
    identityStates = if ("onset" %in% stages) list(onset = "AlivePneum") else list(),
    periodBreaks = config$period_breaks, periodLabels = config$period_labels)
}

estimate_pneumonia_onset <- function(cycles, config) {
  estimated <- estimate_transitions(cycles, config, stages = "onset")$estimates
  onset <- estimated[estimated$stage == "onset" & estimated$from_state == "AliveNoPneum" &
    estimated$to_state == "AlivePneum" & estimated$denominator > 0 &
    estimated$period != tail(config$period_labels, 1L), , drop = FALSE]
  fit <- data.frame(period = onset$period, cycles = onset$denominator, events = onset$count,
    probability = onset$probability, stringsAsFactors = FALSE)
  fit <- fit[order(match(fit$period, config$period_labels)), ]
  if (!nrow(fit)) stop("No complete cycles support onset estimation.")
  last <- max(match(fit$period, config$period_labels))
  if (!all(config$period_labels[seq_len(last)] %in% fit$period)) stop("An internal onset period lacks support.")
  cycle <- seq_len(ceiling(config$economic$horizon_years * config$economic$days_per_year / 7)) - 1L
  period <- config$period_labels[findInterval(cycle * 7, config$period_breaks)]
  lookup <- pmin(match(period, config$period_labels), last)
  list(fit = fit, schedule = data.frame(cycle = cycle,
    p_onset = fit$probability[match(config$period_labels[lookup], fit$period)]))
}


# The primary thesis uses a shared background-mortality schedule for the
# eligible cohort, plus a state- and period-specific excess hazard. This adapter
# prepares those supplied probabilities. TMA performs the cohort propagation.
estimate_mortality <- function(cycles, patients, life_table_file, config) {
  days <- config$economic$cycle_days
  year_days <- config$economic$days_per_year
  max_age <- config$mortality$max_age
  labels <- config$period_labels
  states <- c("AliveNoPneum", "AlivePneum")
  n_cycles <- ceiling(config$economic$horizon_years * year_days / days)
  period_for <- function(cycle) labels[findInterval(cycle * days, config$period_breaks)]
  if (any(!is.finite(patients$age)) || any(patients$age < 0)) {
    stop("Every eligible patient needs a non-negative baseline age; ages are not imputed.")
  }
  if (anyNA(patients$gender_concept_id) || any(!patients$gender_concept_id %in% c(8507, 8532))) {
    stop("The sex-specific life table requires gender_concept_id 8507 or 8532 for every eligible patient.")
  }
  patients$sex <- ifelse(patients$gender_concept_id == 8532, "Female", "Male")

  # The bundled 2025 table supplies exact ages 0-100. Extend its last 21 ages
  # with a log-hazard line, anchored monotonically to the last published value.
  published <- utils::read.csv(life_table_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  life_table <- data.frame(sex = ifelse(published$Sex == "Females", "Female", "Male"),
    age = published$Age, annual_probability = pmin(published[["2025"]], 0.999))
  life_table <- life_table[order(life_table$sex, life_table$age), ]
  for (sex in c("Female", "Male")) {
    observed <- life_table[life_table$sex == sex, ]
    last_age <- max(observed$age)
    fit_data <- observed[observed$age >= max(last_age - 20L, 0L) & observed$annual_probability > 0, ]
    fit_data$log_hazard <- log(-log1p(-fit_data$annual_probability))
    fitted_hazard <- stats::lm(log_hazard ~ age, data = fit_data)
    new_ages <- seq.int(last_age + 1L, max_age)
    predicted <- pmin(-expm1(-exp(stats::predict(fitted_hazard, data.frame(age = new_ages)))), 0.999)
    anchored <- cummax(c(tail(observed$annual_probability, 1L), predicted))[-1L]
    life_table <- rbind(life_table, data.frame(sex = sex, age = new_ages, annual_probability = anchored))
  }
  life_table <- life_table[order(life_table$sex, life_table$age), ]
  life_table_keys <- paste(life_table$sex, life_table$age)
  weekly_probability <- function(age, sex) {
    attained_age <- pmin(pmax(floor(age), 0), max_age)
    annual <- life_table$annual_probability[match(paste(sex, attained_age), life_table_keys)]
    if (anyNA(annual)) stop("An attained age/sex is missing from the life table.")
    1 - (1 - annual)^(days / year_days)
  }

  # TMA matches the background rate to each observed source cycle before
  # estimating the period/state excess hazard. Negative excess is set to zero.
  stage_rows <- stage_observations(cycles, stages = "mortality")
  patient_row <- match(stage_rows$subject_id, patients$subject_id)
  stage_rows$background_age <- pmin(pmax(floor(patients$age[patient_row] + stage_rows$source_day / year_days), 0), max_age)
  stage_rows$sex <- patients$sex[patient_row]
  fitted_stages <- estimate_transitions(cycles, config, stages = "mortality", observations = stage_rows)
  background <- TrajectoryMarkovAnalysis::createModelBundle(
    TrajectoryMarkovAnalysis::createMortalityTransitions(
      data.frame(background_age = life_table$age, sex = life_table$sex,
        death_probability = weekly_probability(life_table$age, life_table$sex)),
      states = c("BackgroundAlive", "BackgroundDeath"), deathState = "BackgroundDeath",
      cycleLength = days, probabilityDays = days, strata = c("background_age", "sex"),
      source = "Supplied Statistics Estonia age/sex probabilities with the study's explicit old-age extension"))
  calibrated <- withCallingHandlers(
    TrajectoryMarkovAnalysis::calibrateExcessMortality(fitted_stages, background,
      deathStates = c("DeathNoHistory", "DeathHistory"), stage = "mortality",
      negativeExcess = "zero", probabilityCeiling = 0.999999)$estimates,
    tma_excess_mortality_floor = function(warning) {
      message(conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  calibrated <- calibrated[calibrated$cycles > 0, ]
  calibration <- data.frame(state = calibrated$from_state, period = calibrated$period,
    period_index = match(calibrated$period, labels), cycles = calibrated$cycles,
    deaths = calibrated$deaths, observed_probability = calibrated$observed_probability,
    background_probability = calibrated$background_probability,
    excess_hazard_per_cycle = calibrated$excess_hazard_per_day * days,
    observed_below_background = calibrated$observed_probability < calibrated$background_probability,
    reconstructed_probability = calibrated$reconstructed_probability,
    used_for_projection = calibrated$period != tail(labels, 1L))
  calibration <- calibration[order(calibration$state, calibration$period_index), ]
  rownames(calibration) <- NULL
  empirical <- cycles[cycles$complete_cycle & cycles$period != tail(labels, 1L), ]
  if (!nrow(empirical)) stop("No complete mortality source cycles exist before the extrapolation period.")
  last_observed_cycle <- max(empirical$cycle)

  # Baseline ages are rounded for this shared demographic projection. Its
  # weights evolve under background mortality alone, identically in both arms.
  demographics <- data.frame(baseline_age = pmin(round(patients$age), max_age), sex = patients$sex, patients = 1L)
  demographics <- stats::aggregate(patients ~ baseline_age + sex, demographics, sum)
  survivors <- demographics$patients
  schedule <- data.frame(cycle = seq_len(n_cycles) - 1L,
    background_probability = NA_real_, externally_surviving = NA_real_)
  for (week in seq_len(n_cycles)) {
    age <- demographics$baseline_age + schedule$cycle[week] * days / year_days
    mortality <- weekly_probability(age, demographics$sex)
    schedule$externally_surviving[week] <- sum(survivors)
    schedule$background_probability[week] <- if (sum(survivors) > 0) sum(survivors * mortality) / sum(survivors) else 1
    survivors <- survivors * (1 - mortality)
  }
  schedule$period <- period_for(schedule$cycle)
  requested_period <- period_for(pmin(schedule$cycle, last_observed_cycle))
  requested_index <- match(requested_period, labels)
  taper_cycles <- ceiling(config$mortality$taper_years * year_days / days)
  elapsed <- pmax(schedule$cycle - last_observed_cycle, 0)
  schedule$excess_multiplier <- pmax(1 - elapsed / taper_cycles, 0)
  for (state in states) {
    available <- calibration[calibration$state == state & calibration$used_for_projection, ]
    available <- available[order(available$period_index), ]
    selected <- findInterval(requested_index, available$period_index)
    if (!nrow(available) || any(selected == 0L)) {
      stop("No earlier observed mortality support for ", state, "; later periods are not borrowed.")
    }
    excess <- available$excess_hazard_per_cycle[selected]
    total_hazard <- -log1p(-pmin(schedule$background_probability, 0.999999)) + excess * schedule$excess_multiplier
    suffix <- if (state == "AliveNoPneum") "no_pneumonia" else "pneumonia"
    schedule[[paste0("p_death_", suffix)]] <- pmin(-expm1(-total_hazard), 0.999999)
    schedule[[paste0("source_period_", suffix)]] <- available$period[selected]
  }
  list(schedule = schedule, calibration = calibration)
}

# TMA allocates the original service duration, including unobserved shares.
# The thesis reward is mean charge per observed weekly row, not a daily rate.

estimate_weekly_costs <- function(histories, config, cost_data) {
  if (!nrow(cost_data) || any(!is.finite(cost_data$amount))) {
    stop("Nonempty Visit Detail total_charge records with finite charges are required; missing charges are not zero.")
  }
  if (anyNA(cost_data$currency_concept_id) || any(cost_data$currency_concept_id != "44818568")) {
    stop("All charges must be EUR (currency_concept_id 44818568); no currency conversion is applied.")
  }
  price_index <- utils::read.csv(config$cost$price_index_file, stringsAsFactors = FALSE)
  cost_data$measure <- rep("total_charge", nrow(cost_data))
  adjusted <- TrajectoryMarkovAnalysis::adjustCostPrices(
    costData = cost_data, priceIndex = price_index, priceYear = config$cost$price_year,
    sourceYear = "event_date", source = config$cost$source)
  fit <- TrajectoryMarkovAnalysis::estimateStateCosts(
    trajectories = histories, costData = adjusted,
    periodBreaks = config$period_breaks, periodLabels = config$period_labels,
    terminalCostMap = c(DeathNoHistory = "AliveNoPneum", DeathHistory = "AlivePneum"), allocation = "duration")
  intervals <- fit$intervalCosts
  rewards <- as.data.frame(dplyr::summarise(
    dplyr::group_by(intervals, .data$period, state = .data$cost_state),
    patient_cycles = dplyr::n(), observed_days = sum(.data$observed_days),
    mean_cycle_charge = mean(.data$total_cost), .groups = "drop"))
  living <- c("AliveNoPneum", "AlivePneum")
  coverage <- table(factor(rewards$period, levels = config$period_labels), factor(rewards$state, levels = living))
  empirical <- seq_len(nrow(coverage) - 1L)
  supported <- empirical[rowSums(coverage[empirical, , drop = FALSE] > 0) == 2L]
  if (!length(supported)) stop("No period has observed weekly cost support for both living states.")
  last <- max(supported)
  if (any(coverage[seq_len(last), , drop = FALSE] == 0L)) stop("Missing internal period/state cost cells; these are not imputed as zero.")
  later <- empirical[empirical > last]
  if (length(later) && any(rowSums(coverage[later, , drop = FALSE] > 0) == 1L)) {
    stop("A later empirical cost period supports only one living state; resolve this support gap before extrapolating.")
  }
  # Both states use the patient-cycle-weighted pooled mean from the last
  # jointly supported period. The final configured period is always extrapolated.
  tail_rows <- rewards$period == config$period_labels[last]
  tail_charge <- stats::weighted.mean(rewards$mean_cycle_charge[tail_rows], rewards$patient_cycles[tail_rows])
  rewards <- rewards[match(rewards$period, config$period_labels) <= last, ]
  rewards$extrapolated <- FALSE
  tail_rewards <- expand.grid(period = config$period_labels[seq.int(last + 1L, length(config$period_labels))],
    state = living, stringsAsFactors = FALSE)
  tail_rewards$patient_cycles <- 0L
  tail_rewards$observed_days <- 0L
  tail_rewards$mean_cycle_charge <- tail_charge
  tail_rewards$extrapolated <- TRUE
  rewards <- rbind(rewards, tail_rewards[names(rewards)])
  rewards <- rewards[order(match(rewards$period, config$period_labels), match(rewards$state, living)), ]
  rownames(rewards) <- NULL
  balance <- data.frame(adjusted_input_charge = sum(adjusted$amount),
    allocated_charge = sum(fit$allocatedCosts$amount), unallocated_charge = sum(fit$unallocatedCosts$amount))
  balance$residual <- balance$adjusted_input_charge - balance$allocated_charge - balance$unallocated_charge
  if (abs(balance$residual) > 1e-8 * max(1, sum(abs(adjusted$amount)))) stop("Allocated and unallocated charges do not reconcile.")
  list(fit = fit, rewards = rewards, balance = balance,
    last_supported_period = config$period_labels[last], pooled_tail_charge = tail_charge)
}

# One TMA model for both arms: pneumonia onset, then mortality. TMA charges
# the post-onset state, credits start-of-week survival and discounts at midpoint.

project_standard_and_intervention <- function(schedule, eligible_n, config) {
  economic <- config$economic
  n_cycles <- nrow(schedule)
  cycles <- seq_len(n_cycles) - 1L
  if (any(schedule$cost_no_pneumonia < 0 | schedule$cost_pneumonia < 0)) {
    stop("Mean weekly charges must be non-negative.")
  }
  states <- c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory")
  living <- states[1:2]
  death_states <- states[3:4]
  identity <- diag(length(states))
  dimnames(identity) <- list(states, states)
  period_labels <- paste0("source_cycle_", cycles)
  period_breaks <- c(cycles * economic$cycle_days, Inf)
  mortality <- lapply(seq_len(n_cycles), function(i) {
    matrix <- identity
    matrix["AliveNoPneum", "AliveNoPneum"] <- 1 - schedule$p_death_no_pneumonia[i]
    matrix["AliveNoPneum", "DeathNoHistory"] <- schedule$p_death_no_pneumonia[i]
    matrix["AlivePneum", "AlivePneum"] <- 1 - schedule$p_death_pneumonia[i]
    matrix["AlivePneum", "DeathHistory"] <- schedule$p_death_pneumonia[i]
    matrix
  })
  names(mortality) <- period_labels

  onset <- lapply(schedule$p_onset, function(probability) {
    matrix <- identity
    matrix["AliveNoPneum", "AliveNoPneum"] <- 1 - probability
    matrix["AliveNoPneum", "AlivePneum"] <- probability
    matrix
  })
  names(onset) <- period_labels
  model <- TrajectoryMarkovAnalysis::composeTransitionStages(
    stages = list(pneumonia_onset = onset, mortality = mortality),
    cycleLength = economic$cycle_days,
    absorbingStates = death_states,
    periodBreaks = period_breaks,
    periodLabels = period_labels
  )

  sum_by_cycle <- function(values, source_cycles) {
    result <- numeric(n_cycles)
    if (length(values)) {
      sums <- rowsum(values, source_cycles, reorder = FALSE)
      result[as.integer(rownames(sums)) + 1L] <- sums[, 1L]
    }
    result
  }
  project_arm <- function(reduction, startup_cost) {
    projection <- TrajectoryMarkovAnalysis::projectCohort(
      transitions = model, nCycles = n_cycles,
      initial = c(AliveNoPneum = 1), cohortSize = eligible_n,
      transitionEffects = if (economic$effect_cycles > 0 && reduction > 0) list(
        prevention = list(stage = "pneumonia_onset", from = "AliveNoPneum",
          to = "AlivePneum", multiplier = 1 - reduction, redistributeTo = "AliveNoPneum",
          endDay = economic$effect_cycles * economic$cycle_days)) else NULL
    )

    # Stage flows retain onsets followed by death within the same cycle.
    flows <- projection$stageFlows
    onset_flows <- flows[flows$stage == "pneumonia_onset", , drop = FALSE]
    new_onset <- onset_flows$from_state == "AliveNoPneum" &
      onset_flows$to_state == "AlivePneum"
    new_death <- flows$stage == "mortality" & flows$from_state %in% living &
      flows$to_state %in% death_states
    occupancy <- projection$occupancy
    alive <- occupancy$state %in% living
    alive_counts <- rowsum(occupancy$expected_count[alive], occupancy$cycle[alive],
                          reorder = FALSE)[as.character(0:n_cycles), 1L]
    rewards <- define_rewards(schedule, config, startup_cost / eligible_n)
    by_cycle <- data.frame(
      cycle = cycles,
      alive_end = unname(alive_counts[seq_len(n_cycles) + 1L]),
      pneumonia_onsets = sum_by_cycle(onset_flows$expected_count[new_onset], onset_flows$source_cycle[new_onset]),
      deaths = sum_by_cycle(flows$expected_count[new_death], flows$source_cycle[new_death])
    )
    cycles_run <- n_cycles
    if (!is.null(economic$stop_alive_threshold)) {
      reached <- which(by_cycle$alive_end <= economic$stop_alive_threshold)
      if (length(reached)) cycles_run <- reached[1L]
    }
    by_cycle <- by_cycle[seq_len(cycles_run), , drop = FALSE]
    list(projection = projection, by_cycle = by_cycle, rewards = rewards)
  }
  arms <- list(
    standard = project_arm(0, 0),
    intervention = project_arm(economic$reduction, eligible_n * economic$intervention_cost)
  )
  list(eligible_n = eligible_n, n_cycles = n_cycles, arms = arms)
}

# Explicit weekly rewards. The amounts are observed weekly means, including
# shortened final rows, rather than daily means multiplied by seven.
define_rewards <- function(schedule, config, upfront = 0) {
  states <- c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory")
  periods <- paste0("source_cycle_", schedule$cycle)
  rewards <- data.frame(period = rep(periods, each = 4L), state = rep(states, nrow(schedule)),
    cost = as.vector(t(cbind(schedule$cost_no_pneumonia, schedule$cost_pneumonia, 0, 0))))
  TrajectoryMarkovAnalysis::createRewardModel(
    stateCosts = rewards, upfrontCosts = c(program = upfront),
    aliveStates = states[1:2], utilities = stats::setNames(c(rep(config$economic$utility, 2), 0, 0), states),
    costTiming = "pneumonia_onset", healthTiming = "start", discountTiming = "midpoint",
    costDiscountRate = config$economic$discount_rate, healthDiscountRate = config$economic$discount_rate,
    daysPerYear = config$economic$days_per_year, currency = "EUR", priceYear = config$cost$price_year,
    periodBreaks = c(schedule$cycle * config$economic$cycle_days, Inf), periodLabels = periods)
}


estimate_model_inputs <- function(c2t, patients, config, cost_data) {
  states <- c("AliveNoPneum", "AlivePneum", "DeathNoHistory", "DeathHistory")
  eligible <- patients[patients$eligible, , drop = FALSE]
  if (!nrow(eligible)) stop("No patients met the bundled dysphagia visit-interval definition.")
  counts <- data.frame(
    all_patients = nrow(patients), eligible_patients = nrow(eligible),
    eligible_deaths = sum(!is.na(eligible$death_date)),
    eligible_pneumonia = sum(!is.na(eligible$pneumonia_date)),
    eligible_same_day_pneumonia_deaths = sum(eligible$pneumonia_date == eligible$death_date, na.rm = TRUE)
  )
  selected <- c2t$trajectories[c2t$trajectories$subject_id %in% eligible$subject_id, , drop = FALSE]
  histories <- TrajectoryMarkovAnalysis::prepareTrajectories(
    selected, cycleLength = 7L, states = states, absorbingStates = states[3:4]
  )
  if (histories$diagnostics$gap_count || histories$diagnostics$unresolved_intervals ||
      length(histories$diagnostics$subjects_without_start) || length(histories$diagnostics$subjects_without_exit)) {
    stop("Resolve missing boundaries, gaps or unresolved intervals before the primary analysis.")
  }
  cycles <- build_weekly_analysis_rows(eligible, histories, config)
  onset <- estimate_pneumonia_onset(cycles, config)
  mortality <- estimate_mortality(cycles, eligible, config$mortality$life_table_file, config)
  costs <- estimate_weekly_costs(histories, config, cost_data)
  schedule <- dplyr::left_join(onset$schedule, mortality$schedule, by = "cycle")
  rewards <- costs$rewards
  for (state in states[1:2]) {
    state_rewards <- rewards[rewards$state == state, ]
    values <- state_rewards$mean_cycle_charge[match(schedule$period, state_rewards$period)]
    if (anyNA(values)) stop("Missing weekly cost reward for ", state, ".")
    column <- if (state == "AliveNoPneum") "cost_no_pneumonia" else "cost_pneumonia"
    schedule[[column]] <- values
  }
  list(population = counts, onset = onset$fit, mortality = mortality,
       costs = costs, schedule = schedule)
}
