# Fixed thesis definitions and assumptions. Partners normally leave these unchanged.
study_settings <- function(study_dir) {
  list(
    output_dir = file.path(study_dir, "results"),
    validation_replicates = 10000L,
    economic_replicates = 10000L,
    seed = 120917L,
    broad_marker_file = file.path(study_dir, "inputs", "broad_marker_specification.csv"),
    cohort_dir = file.path(study_dir, "cohorts"),
    cohort_table_prefix = "tma_study",
    period_breaks = c(0, 7, 35, 91, 189, 371, 735, 1827, 3647, Inf),
    period_labels = c("Week 0", "Weeks 1-4", "Weeks 5-12", "Weeks 13-26",
      "Weeks 27-52", "Year 2", "Years 3-5", "Years 6-10", "Beyond observed follow-up"),
    economic = list(
      reduction = 0.25, effect_cycles = 4L, intervention_cost = 1529,
      utility = 0.68, discount_rate = 0.035, days_per_year = 365.25,
      cycle_days = 7L, horizon_years = 40, stop_alive_threshold = 0.5
    ),
    mortality = list(
      life_table_file = file.path(study_dir, "inputs", "estonia_life_table_2025_raw.csv"),
      taper_years = 5, max_age = 120L
    ),
    cost = list(
      price_index_file = file.path(study_dir, "inputs", "frozen_cpi.csv"),
      price_year = 2026L, currency_concept_id = "44818568",
      source = paste(
        "Frozen thesis index attributed to Statistics Estonia IA002, total CPI, 1997=100.",
        "1998-2025 annual means. 2026 January-August mean, not a full-year value."
      )
    )
  )
}
