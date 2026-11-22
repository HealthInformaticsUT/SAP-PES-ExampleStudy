# SAP and PES example study

This study uses Cohort2Trajectory and TrajectoryMarkovAnalysis to examine stroke-associated pneumonia (SAP) and an economic scenario for pharyngeal electrical stimulation (PES).

Use R 4.1 or later, Java 8 or later for CirceR and SqlRender, and the C++ build tools for your R installation. Install the tested package revisions before running the study:

```r
install.packages(c("remotes", "CDMConnector", "CirceR", "SqlRender",
                   "dplyr", "dbplyr", "jsonlite", "survival", "zip"))

remotes::install_github(
  "HealthInformaticsUT/Cohort2Trajectory@04e4a573c71fe1af919e20dc497dafb953851346",
  upgrade = "never"
)
remotes::install_github(
  "HealthInformaticsUT/TrajectoryMarkovAnalysis@21ca3daf104639ccfb901c145afb39867246255a",
  upgrade = "never"
)
```

These commits provide Cohort2Trajectory 1.3.0.9000 and TrajectoryMarkovAnalysis 1.1.0.9000, including `prepareStageObservations()` for the aligned trajectory views used here. Pinning both packages keeps later repository changes from altering the code used for this study. Restart R after updating packages that were already loaded.

Connect to your OMOP database in R so that `cdm` is available with a writable results schema. Keep database credentials in your own connection setup, then run:

```r
source("/path/to/SAP-PES-ExampleStudy/code2Run.R")

# Review the protected results
study_result$tables$economic_primary
study_result$archive
```

The script also uses CDMConnector, CirceR, SqlRender, dplyr, dbplyr, jsonlite, survival and zip. It checks dependencies before starting and does not install or update packages.

The study uses the bundled cohort definitions, life table, price index and fixed thesis assumptions. It estimates results from your CDM. Counts and estimates need not match the original study. The external mortality and price inputs are Estonian. The fixed utility is 0.68, intervention cost is EUR 1,529 and annual discount rate is 3.5%. Cost analysis requires `cost.total_charge` in EUR linked to `Visit Detail` records with service start and end dates. Recorded service durations are retained and charges are adjusted to 2026 prices.

NIHSS is optional. Before running, you may supply a data frame named `nihss_table` with columns `person_id`, `measurement_date` and `nihss_score`, or place `nihss.csv` with those columns beside `code2Run.R`. Keep person IDs as character values and use dates in `YYYY-MM-DD` format. A supplied `nihss_table` takes precedence over the CSV. Set `nihss_table <- NULL` to omit NIHSS even when the CSV is present. NIHSS is not extracted automatically. The analysis uses one distinct valid score on the stroke index date, excluding conflicting scores. Without this input, analyses that do not require NIHSS still run, and NIHSS-dependent results are reported as unavailable.

The run covers the primary analysis, strict pneumonia-before-death chronology, the broad swallowing-marker subgroup and the Day 8 scenario. It includes observed states, weekly charges, Cox models, economic sensitivities, survival validation and numeric inputs for Figures 5.3, 5.4, 5.6, C.1 and C.2. Day 8 analyses retain the fixed dysphagia group, restricted to patients alive, observed and without pneumonia through Day 7.

The 10,000 validation replicates and 10,000 paired economic simulations are scientific analyses and remain included. Developer tests are not included in this study repository.

Each run creates a new folder under `results/` with a `shareable/` directory and `shareable_results.zip`. Share only that ZIP. It contains aggregate result tables, figure data, `thesis_results.rds` and a short reading guide, with counts of 1–4 and linked values protected. Observed figure cells follow the agreed direct masking rule described in that guide. No patient rows or fitted model objects enter the shareable archive. Cohort2Trajectory may create temporary patient-level files locally while running. Keep the working environment on your secure server.

To follow the analysis, start with `R/run.R`. It reads in the same order as the study: define patients, create trajectories, estimate inputs, fit survival models, calculate economic outcomes, then export.

| File | What it does |
| --- | --- |
| `R/data.R` | Generates cohorts, defines states and reads patient characteristics. |
| `R/trajectories.R` | Estimates weekly transitions, mortality and costs, then projects trajectories. |
| `R/observed.R` | Summarises recorded care and compares observed with modelled survival. |
| `R/survival.R` | Fits the Cox models and optional NIHSS analyses. |
| `R/economics.R` | Calculates costs, QALYs, ICERs and economic sensitivities. |
| `R/export.R` | Applies disclosure rules and writes tables and figure data. |

C2T creates two aligned views of the same cohort facts in `R/data.R`, using the seven-day cycle length from `assumptions.R`. The final view retains the four clinical states and absorbing death states used for costs and follow-up summaries. The history view uses the two living states with pneumonia priority, exposing pneumonia even when death has the final label in that week. Its target ends at the final view's last clinical interval, preserving the same weeks without adding post-death follow-up. Both views carry their definitions into TMA. `prepareStageObservations()` supplies onset followed by mortality within each cycle. A source cycle enters transition estimation only when the next cycle is observed or the current cycle ends in death; a final censored week remains excluded even if it spans seven days. Exact event dates remain available for Cox analyses and chronology sensitivities.

The primary projection defines weekly charges and utilities beside the transition model, using the same states and periods. `summarise_economic_outcomes()` then sets the programme cost and discount rates for each scenario through `updateRewardModel()`. `outcomeTotals()` reads discounted costs, life-years and QALYs without accessing the result object's internal tables.

`assumptions.R` holds the fixed model settings. `cohorts/` and `inputs/` hold the definitions and external inputs. There are no recovery scripts or developer tests. A failed run is restarted from the beginning.

CohortContrast and previously published multi-database results are outside this local study.

## Licence and input sources

The study code is licensed under [Apache-2.0](LICENSE).

The bundled Statistics Estonia data retain their [CC BY-SA 4.0 licence](https://creativecommons.org/licenses/by-sa/4.0/):

- `inputs/estonia_life_table_2025_raw.csv` contains the 2025 male and female probability-of-dying values at ages 0 to 100 from [Statistics Estonia table RV046](https://andmed.stat.ee/en/stat/rahvastik__rahvastikunaitajad-ja-koosseis__demograafilised-pehinaitajad/RV046). Values beyond the published ages are model assumptions implemented in the analysis.
- `inputs/frozen_cpi.csv` contains averages derived from the monthly total consumer price index in [Statistics Estonia table IA002](https://andmed.stat.ee/en/stat/majandus__hinnad/IA002). The bundled 2026 value uses eight months, as recorded in `months_observed`. These inputs are kept fixed to reproduce the study.

`inputs/broad_marker_specification.csv` defines the study's recorded swallowing and care markers. The cohort JSON files contain definitions and concept metadata, not patient records.
