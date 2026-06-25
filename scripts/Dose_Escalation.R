###############################################################################
# CRM operating characteristics simulation
#
# Purpose:
#   This script evaluates the operating characteristics of a Bayesian CRM
#   dose-escalation design using crmPack.
#
#   The script:
#     1. Defines a discrete dose grid and several true toxicity scenarios.
#     2. Specifies the CRM model, dose-selection rule, escalation restrictions,
#        and stopping rules.
#     3. Runs simulations under each toxicity scenario.
#     4. Summarises operating characteristics, including:
#          - total sample size,
#          - escalation sample size,
#          - backfill sample size,
#          - true MTD,
#          - detected MTD distribution,
#          - probability of selecting the true MTD,
#          - probability of overdosing.
#     5. Produces figures and Excel outputs.
#
# Notes:
#   - The dose grid is assumed to be fixed and prespecified.
#   - Overdosing is defined as treatment at a dose whose true DLT probability
#     exceeds the target toxicity rate.
#   - The detected MTD is extracted from the crmPack simulation object when
#     available. If no explicit selected/recommended MTD slot is found, the
#     script falls back to the last non-backfill escalation dose.
#
# Required packages:
#   crmPack, dplyr, tidyr, tibble, ggplot2, purrr, scales, parallelly
#
# Optional package:
#   writexl, for exporting Excel files.
#
###############################################################################


# =============================================================================
# 0. Package loading
# =============================================================================

required_packages <- c(
  "crmPack",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "purrr",
  "scales",
  "parallelly"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are missing: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(crmPack)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(parallelly)
})


# =============================================================================
# 1. User-defined settings
# =============================================================================

# ---- General simulation settings ----

nsim <- 1000
seed <- 8052026

# Use all available cores by default.
# If running on a shared machine or cluster, consider replacing this with a
# smaller fixed value, e.g. ncores <- 4.
ncores <- parallelly::availableCores()

# Number of cohorts required near the recommended dose for some stopping rules.
num_cohorts_near_dose <- 3


# ---- Dose grid and target toxicity ----

doses <- c(50, 100, 200, 300, 400)

dose_grid <- doses
min_dose <- min(dose_grid)
max_dose <- max(dose_grid)

# Target DLT probability.
tox_thr <- 0.25


# ---- Scenario order for tables and plots ----

scenario_order <- c(
  "Toxic",
  "First",
  "Second",
  "Third",
  "Fourth",
  "Last",
  "Safe"
)


# ---- True DLT probability scenarios ----
#
# Each vector corresponds to the true DLT probabilities at the dose levels:
#   50, 100, 200, 300, 400 mg.
#
# Scenario interpretation:
#   Toxic : all doses are above target.
#   First : first dose is closest to target.
#   Second: second dose is closest to target.
#   Third : third dose is closest to target.
#   Fourth: fourth dose is closest to target.
#   Last  : last dose is closest to target.
#   Safe  : all doses are at or below target, with last dose close to target.

true_scenarios <- list(
  Toxic  = c(0.35, 0.45, 0.55, 0.65, 0.75),
  First  = c(0.25, 0.35, 0.45, 0.55, 0.65),
  Second = c(0.15, 0.25, 0.35, 0.45, 0.55),
  Third  = c(0.05, 0.15, 0.25, 0.35, 0.45),
  Fourth = c(0.04, 0.10, 0.18, 0.25, 0.35),
  Last   = c(0.02, 0.06, 0.12, 0.20, 0.25),
  Safe   = c(0.01, 0.04, 0.10, 0.15, 0.20)
)


# ---- Efficacy probabilities for optional truthResponse ----
#
# These values are only used if efficacy response simulation is enabled through
# truthResponse. They do not drive dose escalation decisions in this script.

eff_probs <- c(
  "50"  = 0.10,
  "100" = 0.35,
  "200" = 0.50,
  "300" = 0.60,
  "400" = 0.80
)


# ---- Corporate-style plotting palette ----

bayer_pal <- c(
  "Toxic"  = "#D30F4B",
  "First"  = "#10384F",
  "Second" = "#007CBF",
  "Third"  = "#00607E",
  "Fourth" = "#4F6D7A",
  "Last"   = "#286436",
  "Safe"   = "#624963"
)


# ---- Output folders ----

output_dir <- "outputs"
table_dir  <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE)
dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 2. Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# Safe slot extraction
# -----------------------------------------------------------------------------

slot_or_null <- function(object, slot_name) {
  if (slot_name %in% slotNames(object)) {
    slot(object, slot_name)
  } else {
    NULL
  }
}


# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------

fmt_pct <- function(x) {
  sprintf("%.1f%%", 100 * as.numeric(x))
}

fmt_num <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f"), as.numeric(x))
}


# -----------------------------------------------------------------------------
# Summary statistics helper
# -----------------------------------------------------------------------------

summarize_OC <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd   = stats::sd(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    p10  = as.numeric(stats::quantile(x, 0.10, na.rm = TRUE)),
    p90  = as.numeric(stats::quantile(x, 0.90, na.rm = TRUE))
  )
}


# -----------------------------------------------------------------------------
# crmPack truth-function builder
# -----------------------------------------------------------------------------

make_truth_function <- function(dose_grid, true_probs) {
  stopifnot(length(dose_grid) == length(true_probs))
  
  truth_matrix <- cbind(dose_grid, true_probs)
  
  function(dose) {
    truth_matrix[match(dose, truth_matrix[, 1]), 2]
  }
}


# -----------------------------------------------------------------------------
# Sample-size extraction helpers
# -----------------------------------------------------------------------------

extract_trial_N_total <- function(data_list) {
  vapply(
    data_list,
    function(obj) length(obj@y),
    integer(1)
  )
}

extract_trial_N_backfill <- function(data_list) {
  vapply(
    data_list,
    function(obj) {
      backfilled <- slot_or_null(obj, "backfilled")
      
      if (is.null(backfilled)) {
        return(0L)
      }
      
      as.integer(sum(backfilled, na.rm = TRUE))
    },
    integer(1)
  )
}

extract_trial_N_escalation <- function(data_list) {
  vapply(
    data_list,
    function(obj) {
      backfilled <- slot_or_null(obj, "backfilled")
      
      if (is.null(backfilled)) {
        return(as.integer(length(obj@y)))
      }
      
      as.integer(sum(!backfilled, na.rm = TRUE))
    },
    integer(1)
  )
}


# -----------------------------------------------------------------------------
# Approximate duration helper
# -----------------------------------------------------------------------------
#
# Approximation:
#   The escalation duration is approximated by the maximum observed escalation
#   cohort number times the cycle length.
#
#   If backfill is active, a fixed number of backfill cycles is included and the
#   total duration is approximated as the maximum of escalation and backfill
#   duration.
#
# This is a pragmatic approximation and should be adapted if trial-specific
# recruitment assumptions are available.

approx_duration_months_combined <- function(
    data_list,
    cycle_days = 28,
    bf_cycles = 5
) {
  
  vapply(
    data_list,
    function(obj) {
      
      backfilled <- slot_or_null(obj, "backfilled")
      cohort     <- slot_or_null(obj, "cohort")
      
      if (is.null(backfilled)) {
        backfilled <- rep(FALSE, length(obj@y))
      }
      
      if (is.null(cohort)) {
        n_cycles_escalation <- 0
      } else {
        escalation_cohorts <- cohort[!backfilled]
        
        n_cycles_escalation <- if (length(escalation_cohorts) == 0) {
          0
        } else {
          max(escalation_cohorts, na.rm = TRUE)
        }
      }
      
      has_backfill <- any(backfilled, na.rm = TRUE)
      n_cycles_backfill <- if (has_backfill) bf_cycles else 0
      
      total_cycles <- max(n_cycles_escalation, n_cycles_backfill)
      
      (total_cycles * cycle_days) / 30.4375
    },
    numeric(1)
  )
}


# -----------------------------------------------------------------------------
# Backfill dose distribution helper
# -----------------------------------------------------------------------------

extract_backfill_dose_distribution <- function(data_list, dose_grid) {
  
  lapply(
    data_list,
    function(obj) {
      
      backfilled <- slot_or_null(obj, "backfilled")
      
      if (is.null(backfilled)) {
        out <- rep(0, length(dose_grid))
        names(out) <- as.character(dose_grid)
        return(out)
      }
      
      backfill_doses <- obj@x[backfilled]
      
      if (length(backfill_doses) == 0) {
        out <- rep(0, length(dose_grid))
        names(out) <- as.character(dose_grid)
        return(out)
      }
      
      tab <- table(factor(backfill_doses, levels = dose_grid))
      out <- as.numeric(tab) / sum(tab)
      names(out) <- as.character(dose_grid)
      
      out
    }
  )
}


# -----------------------------------------------------------------------------
# True MTD helpers
# -----------------------------------------------------------------------------

get_true_mtd <- function(true_probs, dose_grid, tox_thr) {
  dose_grid[which.min(abs(true_probs - tox_thr))]
}

get_true_dlt_at_dose <- function(dose, true_probs, dose_grid) {
  true_probs[match(dose, dose_grid)]
}


# -----------------------------------------------------------------------------
# Detected MTD extraction helper
# -----------------------------------------------------------------------------
#
# The exact storage of selected/recommended doses may differ depending on the
# crmPack object version and simulation settings.
#
# This function first tries several plausible simulation-object slots. If none
# are found, the fallback is the last non-backfill escalation dose.

extract_detected_mtd_vector <- function(sims, dose_grid) {
  
  data_list <- slot_or_null(sims, "data")
  
  if (is.null(data_list)) {
    stop("Could not find @data slot in simulation object.")
  }
  
  nsim_expected <- length(data_list)
  
  candidate_slots <- c(
    "dosesSelected",
    "doseSelected",
    "selectedDoses",
    "selectedDose",
    "recommendedDoses",
    "recommendedDose",
    "mtd",
    "mtds"
  )
  
  available_slots <- slotNames(sims)
  
  for (sl in intersect(candidate_slots, available_slots)) {
    
    obj <- slot(sims, sl)
    
    if (is.numeric(obj) && length(obj) >= nsim_expected) {
      selected <- as.numeric(obj[seq_len(nsim_expected)])
      
      if (all(is.na(selected) | selected %in% dose_grid)) {
        message("Detected MTD extracted from slot: @", sl)
        return(selected)
      }
    }
    
    if (is.list(obj) && length(obj) >= nsim_expected) {
      selected <- purrr::map_dbl(
        obj[seq_len(nsim_expected)],
        function(x) {
          if (length(x) == 0) {
            return(NA_real_)
          }
          as.numeric(x[1])
        }
      )
      
      if (all(is.na(selected) | selected %in% dose_grid)) {
        message("Detected MTD extracted from list slot: @", sl)
        return(selected)
      }
    }
  }
  
  message(
    "No explicit detected/recommended MTD slot found. ",
    "Fallback used: last non-backfill escalation dose per simulation."
  )
  
  vapply(
    data_list,
    function(obj) {
      
      x <- obj@x
      backfilled <- slot_or_null(obj, "backfilled")
      
      if (is.null(backfilled)) {
        backfilled <- rep(FALSE, length(x))
      }
      
      x_escalation <- x[!backfilled]
      
      if (length(x_escalation) == 0) {
        return(as.numeric(tail(x, 1)))
      }
      
      as.numeric(tail(x_escalation, 1))
    },
    numeric(1)
  )
}


# -----------------------------------------------------------------------------
# Mode dose helper
# -----------------------------------------------------------------------------

get_mode_dose <- function(x, dose_grid) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  tab <- table(factor(x, levels = dose_grid))
  as.numeric(names(tab)[which.max(tab)])
}


# -----------------------------------------------------------------------------
# OC table row builder
# -----------------------------------------------------------------------------

make_OC_rows <- function(
    OC_raw,
    metric,
    criterion_label,
    integer_metric = TRUE
) {
  
  if (integer_metric) {
    
    fmt_mean_sd <- function(x) {
      sprintf("%.1f (%.1f)", x["mean"], x["sd"])
    }
    
    fmt_minmax <- function(x) {
      sprintf(
        "[%d–%d]",
        as.integer(x["min"]),
        as.integer(x["max"])
      )
    }
    
    fmt_p1090 <- function(x) {
      sprintf(
        "[%d–%d]",
        as.integer(round(x["p10"])),
        as.integer(round(x["p90"]))
      )
    }
    
  } else {
    
    fmt_mean_sd <- function(x) {
      sprintf("%.1f (%.1f)", x["mean"], x["sd"])
    }
    
    fmt_minmax <- function(x) {
      sprintf("[%.1f–%.1f]", x["min"], x["max"])
    }
    
    fmt_p1090 <- function(x) {
      sprintf("[%.1f–%.1f]", x["p10"], x["p90"])
    }
  }
  
  rbind(
    c(
      Criterion = criterion_label,
      Category = "mean (sd)",
      sapply(OC_raw, function(x) fmt_mean_sd(x[[metric]]))
    ),
    c(
      Criterion = criterion_label,
      Category = "min–max",
      sapply(OC_raw, function(x) fmt_minmax(x[[metric]]))
    ),
    c(
      Criterion = criterion_label,
      Category = "p10–p90",
      sapply(OC_raw, function(x) fmt_p1090(x[[metric]]))
    )
  )
}


# -----------------------------------------------------------------------------
# MTD distribution table helper
# -----------------------------------------------------------------------------

build_mtd_distribution_df <- function(
    scenario_results,
    dose_grid,
    scenario_order
) {
  
  purrr::imap_dfr(
    scenario_results,
    function(sim_obj, scenario_name) {
      
      selected_mtd <- extract_detected_mtd_vector(
        sims = sim_obj,
        dose_grid = dose_grid
      )
      
      tibble(
        scenario = scenario_name,
        dose = selected_mtd
      ) %>%
        count(scenario, dose, name = "n") %>%
        complete(scenario, dose = dose_grid, fill = list(n = 0)) %>%
        group_by(scenario) %>%
        mutate(percent = 100 * n / sum(n)) %>%
        ungroup()
    }
  ) %>%
    mutate(
      scenario = factor(scenario, levels = scenario_order)
    )
}


# -----------------------------------------------------------------------------
# Probability of overdosing based on selected MTD
# -----------------------------------------------------------------------------
#
# Definition used here:
#   PoD = P(selected MTD > true MTD)
#
# This is different from:
#   P(selected MTD has true DLT probability > target)
#
# Both are useful, but they answer slightly different questions.

compute_pod_selected_above_true_mtd <- function(mtd_df, true_mtd_tbl) {
  
  mtd_df %>%
    left_join(true_mtd_tbl, by = "scenario") %>%
    mutate(overdose = dose > true_mtd_dose) %>%
    group_by(scenario, true_mtd_dose, true_mtd_true_dlt) %>%
    summarise(
      PoD = sum(percent[overdose], na.rm = TRUE),
      .groups = "drop"
    )
}


# =============================================================================
# 3. Scenario preparation
# =============================================================================

truth_fun <- purrr::imap(
  true_scenarios,
  ~ make_truth_function(
    dose_grid = dose_grid,
    true_probs = .x
  )
)

truthResponse_fun <- function(dose) {
  eff_probs[as.character(dose)]
}


# =============================================================================
# 4. Plot true toxicity curves
# =============================================================================

tox_df <- purrr::imap_dfr(
  true_scenarios,
  ~ tibble(
    scenario = .y,
    dose = dose_grid,
    rate = .x
  )
) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    dose = as.numeric(dose)
  ) %>%
  arrange(scenario, dose)

p_true_toxicity_curves <- ggplot(
  tox_df,
  aes(
    x = dose,
    y = rate,
    color = scenario,
    group = scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = tox_thr,
    linetype = "dashed",
    color = "#D30F4B",
    linewidth = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_continuous(
    breaks = dose_grid
  ) +
  scale_color_manual(
    values = bayer_pal,
    breaks = scenario_order,
    drop = FALSE
  ) +
  labs(
    title = "True toxicity curves by scenario",
    x = "Dose (mg)",
    y = "True DLT probability",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(figure_dir, "true_toxicity_curves.png"),
  plot = p_true_toxicity_curves,
  width = 8,
  height = 5,
  dpi = 300
)


# =============================================================================
# 5. CRM design specification
# =============================================================================

empty_data <- Data(doseGrid = dose_grid)

# ---- CRM model ----
#
# LogisticLogNormalSub model used in crmPack.
# The prior parameters below should be aligned with the protocol/SAP and
# calibration work. If alternative prior calibrations are considered, document
# them explicitly.

initial_model <- LogisticLogNormalSub(
  mean = c(-1.33, -5.2),
  cov = matrix(
    c(
      0.7, 0,
      0,   0.7
    ),
    nrow = 2
  ),
  ref_dose = mean(dose_grid)
)


# ---- Next-best rule ----
#
# The next recommended dose is the dose with estimated toxicity probability
# closest to the target toxicity rate.

next_best <- NextBestMinDist(target = tox_thr)


# ---- Cohort size ----

cohort_size <- CohortSizeConst(size = 3)


# ---- Stopping rule 1: lowest dose too toxic ----
#
# Stop if:
#   - at the lowest dose, P(DLT probability in [0.40, 1]) > 0.70, and
#   - at least one cohort has been treated at that dose.

stopping_lowest_dose_nogo <- (
  StoppingSpecificDose(
    rule = StoppingTargetProb(target = c(0.40, 1), prob = 0.70),
    dose = min_dose
  ) &
    StoppingSpecificDose(
      rule = StoppingPatientsNearDose(nPatients = 3, percentage = 0),
      dose = min_dose
    )
)


# ---- Stopping rule 2: highest dose clearly safe ----
#
# Stop if:
#   - at the highest dose, P(DLT probability in [0, target]) > 0.80, and
#   - at least one cohort has been treated at that dose.

stopping_highest_dose_safe <- (
  StoppingSpecificDose(
    rule = StoppingTargetProb(target = c(0, tox_thr), prob = 0.80),
    dose = max_dose
  ) &
    StoppingSpecificDose(
      rule = StoppingPatientsNearDose(nPatients = 3, percentage = 0),
      dose = max_dose
    )
)


# ---- Stopping rule 3: sufficient information near dose ----

stopping_sufficient_information <- StoppingCohortsNearDose(
  nCohorts = num_cohorts_near_dose,
  percentage = 0
)


# ---- Stopping rule 4: precision-based stopping ----

stopping_precision <- (
  StoppingMTDCV(target = tox_thr, thresh_cv = 30) &
    StoppingCohortsNearDose(
      nCohorts = num_cohorts_near_dose,
      percentage = 0
    )
)


# ---- Combined stopping rule ----

stopping_trial <- (
  stopping_lowest_dose_nogo |
    stopping_highest_dose_safe |
    stopping_sufficient_information |
    stopping_precision |
    StoppingMissingDose()
)


# ---- Dose-increment restrictions ----
#
# No dose skipping:
#   The next dose cannot be more than one dose level above the last given dose.
#
# Hard safety restriction:
#   Escalation is limited based on posterior overdose risk.

increment_no_dose_skipping <- IncrementsDoseLevels(
  levels = 1,
  basis_level = "last"
)

increment_hard_safety <- IncrementsHSRBeta(
  target = tox_thr,
  prob = 0.95
)

increment_rule <- IncrementsMin(
  increments_list = list(
    increment_no_dose_skipping,
    increment_hard_safety
  )
)


# ---- Main design object ----

design <- Design(
  model = initial_model,
  nextBest = next_best,
  stopping = stopping_trial,
  increments = increment_rule,
  cohort_size = cohort_size,
  data = empty_data,
  startingDose = min_dose
)


# ---- Optional backfill specification ----
#
# Backfill is defined but not activated by default.
# To activate backfill, uncomment:
#   design@backfill <- backfill

backfill <- Backfill(
  cohort_size = CohortSizeConst(3),
  max_size = 2 * 15,
  opening = OpeningMinDose(min_dose = 200),
  recruitment = RecruitmentUnlimited(),
  priority = "lowest"
)

# design@backfill <- backfill


# ---- MCMC options ----
#
# rng_kind = "Mersenne-Twister" is a standard pseudo-random number generator in R.
# rng_seed fixes the MCMC seed to support reproducibility.

mcmc_options <- McmcOptions(
  burnin = 1000,
  step = 2,
  samples = 5000,
  rng_kind = "Mersenne-Twister",
  rng_seed = 12345
)


# =============================================================================
# 6. Simulation execution
# =============================================================================

run_one_scenario <- function(
    design,
    truth_fun,
    truthResponse_fun,
    nsim,
    seed,
    mcmc_options,
    ncores
) {
  
  simulate(
    design,
    truth = truth_fun,
    truthResponse = truthResponse_fun,
    nsim = nsim,
    seed = seed,
    mcmcOptions = mcmc_options,
    parallel = TRUE,
    nCores = ncores
  )
}


scenario_results <- purrr::imap(
  truth_fun,
  function(current_truth_fun, scenario_name) {
    
    message(">>> Running scenario: ", scenario_name)
    
    run_one_scenario(
      design = design,
      truth_fun = current_truth_fun,
      truthResponse_fun = truthResponse_fun,
      nsim = nsim,
      seed = seed,
      mcmc_options = mcmc_options,
      ncores = ncores
    )
  }
)

names(scenario_results) <- names(truth_fun)


# =============================================================================
# 7. Operating-characteristic extraction
# =============================================================================

OC_raw <- lapply(
  names(scenario_results),
  function(scenario_name) {
    
    sims <- scenario_results[[scenario_name]]
    data_list <- slot_or_null(sims, "data")
    
    if (is.null(data_list)) {
      stop("Could not find @data slot for scenario: ", scenario_name)
    }
    
    true_probs <- true_scenarios[[scenario_name]]
    
    # ---- True MTD ----
    
    true_mtd <- get_true_mtd(
      true_probs = true_probs,
      dose_grid = dose_grid,
      tox_thr = tox_thr
    )
    
    true_mtd_dlt <- get_true_dlt_at_dose(
      dose = true_mtd,
      true_probs = true_probs,
      dose_grid = dose_grid
    )
    
    
    # ---- Sample size ----
    
    N_total <- extract_trial_N_total(data_list)
    N_escal <- extract_trial_N_escalation(data_list)
    N_bf    <- extract_trial_N_backfill(data_list)
    
    duration_months <- approx_duration_months_combined(
      data_list,
      cycle_days = 28,
      bf_cycles = 5
    )
    
    
    # ---- Backfill dose distribution ----
    
    bf_dist_list <- extract_backfill_dose_distribution(
      data_list,
      dose_grid = dose_grid
    )
    
    bf_dist_mat <- do.call(rbind, bf_dist_list)
    bf_dist_mean <- colMeans(bf_dist_mat)
    
    
    # ---- Global overdosing ----
    #
    # Overdosing is defined as treatment at any dose whose true DLT probability
    # exceeds the target toxicity rate.
    
    any_overdose_total <- vapply(
      data_list,
      function(obj) {
        treated_doses <- obj@x
        true_dlt_at_treated_doses <- true_probs[match(treated_doses, dose_grid)]
        
        any(true_dlt_at_treated_doses > tox_thr, na.rm = TRUE)
      },
      logical(1)
    )
    
    prop_patients_overdose_total <- vapply(
      data_list,
      function(obj) {
        treated_doses <- obj@x
        true_dlt_at_treated_doses <- true_probs[match(treated_doses, dose_grid)]
        
        mean(true_dlt_at_treated_doses > tox_thr, na.rm = TRUE)
      },
      numeric(1)
    )
    
    any_overdose_escal <- vapply(
      data_list,
      function(obj) {
        
        treated_doses <- obj@x
        backfilled <- slot_or_null(obj, "backfilled")
        
        if (is.null(backfilled)) {
          backfilled <- rep(FALSE, length(treated_doses))
        }
        
        escalation_doses <- treated_doses[!backfilled]
        
        if (length(escalation_doses) == 0) {
          return(FALSE)
        }
        
        true_dlt_at_escalation_doses <- true_probs[
          match(escalation_doses, dose_grid)
        ]
        
        any(true_dlt_at_escalation_doses > tox_thr, na.rm = TRUE)
      },
      logical(1)
    )
    
    prop_patients_overdose_escal <- vapply(
      data_list,
      function(obj) {
        
        treated_doses <- obj@x
        backfilled <- slot_or_null(obj, "backfilled")
        
        if (is.null(backfilled)) {
          backfilled <- rep(FALSE, length(treated_doses))
        }
        
        escalation_doses <- treated_doses[!backfilled]
        
        if (length(escalation_doses) == 0) {
          return(NA_real_)
        }
        
        true_dlt_at_escalation_doses <- true_probs[
          match(escalation_doses, dose_grid)
        ]
        
        mean(true_dlt_at_escalation_doses > tox_thr, na.rm = TRUE)
      },
      numeric(1)
    )
    
    
    # ---- Detected MTD ----
    
    detected_mtd_vec <- extract_detected_mtd_vector(
      sims = sims,
      dose_grid = dose_grid
    )
    
    detected_mtd_true_dlt <- true_probs[
      match(detected_mtd_vec, dose_grid)
    ]
    
    detected_mtd_mode <- get_mode_dose(
      x = detected_mtd_vec,
      dose_grid = dose_grid
    )
    
    prob_detect_true_mtd <- mean(
      detected_mtd_vec == true_mtd,
      na.rm = TRUE
    )
    
    prob_detect_overdose <- mean(
      detected_mtd_true_dlt > tox_thr,
      na.rm = TRUE
    )
    
    detected_mtd_dist <- table(
      factor(detected_mtd_vec, levels = dose_grid)
    ) / length(detected_mtd_vec)
    
    
    # ---- Scenario-level OC object ----
    
    list(
      scenario = scenario_name,
      
      N_total = summarize_OC(N_total),
      N_escal = summarize_OC(N_escal),
      N_bf    = summarize_OC(N_bf),
      Dur     = summarize_OC(duration_months),
      
      bf_dist_mean = bf_dist_mean,
      
      true_mtd = true_mtd,
      true_mtd_dlt = true_mtd_dlt,
      
      prob_any_overdose_total =
        mean(any_overdose_total, na.rm = TRUE),
      
      mean_prop_patients_overdose_total =
        mean(prop_patients_overdose_total, na.rm = TRUE),
      
      prob_any_overdose_escal =
        mean(any_overdose_escal, na.rm = TRUE),
      
      mean_prop_patients_overdose_escal =
        mean(prop_patients_overdose_escal, na.rm = TRUE),
      
      detected_mtd_mode = detected_mtd_mode,
      prob_detect_true_mtd = prob_detect_true_mtd,
      prob_detect_overdose = prob_detect_overdose,
      detected_mtd_dist = detected_mtd_dist
    )
  }
)

names(OC_raw) <- names(scenario_results)


# =============================================================================
# 8. Operating-characteristic tables
# =============================================================================

# ---- Sample-size table ----

OC_N_total <- data.frame(
  make_OC_rows(
    OC_raw = OC_raw,
    metric = "N_total",
    criterion_label = "Total N (escalation + backfill)",
    integer_metric = TRUE
  ),
  check.names = FALSE
)

OC_N_escal <- data.frame(
  make_OC_rows(
    OC_raw = OC_raw,
    metric = "N_escal",
    criterion_label = "N in escalation (non-backfill)",
    integer_metric = TRUE
  ),
  check.names = FALSE
)

OC_N_bf <- data.frame(
  make_OC_rows(
    OC_raw = OC_raw,
    metric = "N_bf",
    criterion_label = "N in backfill",
    integer_metric = TRUE
  ),
  check.names = FALSE
)

OC_Dur <- data.frame(
  make_OC_rows(
    OC_raw = OC_raw,
    metric = "Dur",
    criterion_label = "Approx. duration (months)",
    integer_metric = FALSE
  ),
  check.names = FALSE
)

OC_table_sample_size <- dplyr::bind_rows(
  OC_N_total,
  OC_N_escal,
  OC_N_bf
  # Uncomment if duration should be included in the final table:
  # , OC_Dur
)


# ---- MTD detection and overdosing table ----

OC_mtd_overdose_summary <- data.frame(
  rbind(
    
    c(
      Criterion = "True scenario",
      Category  = "True MTD (mg)",
      sapply(
        OC_raw,
        function(x) as.character(x$true_mtd)
      )
    ),
    
    c(
      Criterion = "True scenario",
      Category  = "True DLT at true MTD",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$true_mtd_dlt)
      )
    ),
    
    c(
      Criterion = "MTD detection",
      Category  = "Most frequent detected MTD (mg)",
      sapply(
        OC_raw,
        function(x) as.character(x$detected_mtd_mode)
      )
    ),
    
    c(
      Criterion = "MTD detection",
      Category  = "P(detect true MTD)",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$prob_detect_true_mtd)
      )
    ),
    
    c(
      Criterion = "MTD detection",
      Category  = "P(detected MTD is overdosing)",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$prob_detect_overdose)
      )
    ),
    
    c(
      Criterion = "Global overdosing",
      Category  = "P(any patient treated at overdosing dose), total",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$prob_any_overdose_total)
      )
    ),
    
    c(
      Criterion = "Global overdosing",
      Category  = "Mean % patients treated at overdosing dose, total",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$mean_prop_patients_overdose_total)
      )
    ),
    
    c(
      Criterion = "Global overdosing",
      Category  = "P(any patient treated at overdosing dose), escalation only",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$prob_any_overdose_escal)
      )
    ),
    
    c(
      Criterion = "Global overdosing",
      Category  = "Mean % patients treated at overdosing dose, escalation only",
      sapply(
        OC_raw,
        function(x) fmt_pct(x$mean_prop_patients_overdose_escal)
      )
    )
  ),
  check.names = FALSE
)

OC_mtd_detected_distribution <- do.call(
  rbind,
  lapply(
    dose_grid,
    function(dose_value) {
      
      c(
        Criterion = "Detected MTD distribution",
        Category  = paste0("P(detected MTD = ", dose_value, " mg)"),
        sapply(
          OC_raw,
          function(x) {
            
            value <- x$detected_mtd_dist[as.character(dose_value)]
            
            if (is.na(value)) {
              value <- 0
            }
            
            fmt_pct(as.numeric(value))
          }
        )
      )
    }
  )
) %>%
  data.frame(check.names = FALSE)

OC_table_mtd_overdose <- dplyr::bind_rows(
  OC_mtd_overdose_summary,
  OC_mtd_detected_distribution
)


# ---- Combined table ----

OC_table_combined <- dplyr::bind_rows(
  OC_table_sample_size,
  OC_table_mtd_overdose
)

# Backward-compatible object name.
OC_table <- OC_table_combined


# =============================================================================
# 9. Export tables
# =============================================================================

csv_output_path <- file.path(
  table_dir,
  "CRM_operating_characteristics_combined.csv"
)

write.csv(
  OC_table_combined,
  file = csv_output_path,
  row.names = FALSE
)

message("CSV file exported: ", csv_output_path)


if (requireNamespace("writexl", quietly = TRUE)) {
  
  excel_output_path <- file.path(
    table_dir,
    "CRM_operating_characteristics_final.xlsx"
  )
  
  writexl::write_xlsx(
    list(
      "OC_sample_size" = OC_table_sample_size,
      "OC_MTD_overdose" = OC_table_mtd_overdose,
      "OC_combined" = OC_table_combined
    ),
    path = excel_output_path
  )
  
  message("Excel file exported: ", excel_output_path)
  
} else {
  
  message(
    "Package 'writexl' is not installed. ",
    "Install it with install.packages('writexl') if Excel export is required."
  )
}


# =============================================================================
# 10. Figure: distribution of total sample size
# =============================================================================

plot_df_total_N <- purrr::imap_dfr(
  scenario_results,
  function(sim_obj, scenario_name) {
    
    tibble(
      scenario = scenario_name,
      sim = seq_along(sim_obj@data),
      n_patients = purrr::map_int(sim_obj@data, ~ length(.x@y))
    )
  }
) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order)
  )

summary_df_total_N <- plot_df_total_N %>%
  group_by(scenario) %>%
  summarise(
    mean_n   = mean(n_patients),
    median_n = median(n_patients),
    p10      = as.numeric(quantile(n_patients, 0.10)),
    p90      = as.numeric(quantile(n_patients, 0.90)),
    .groups  = "drop"
  )

binwidth_use <- 3

scenario_pos <- tibble(
  scenario = factor(scenario_order, levels = scenario_order),
  scen_id  = seq_along(scenario_order)
)

hist_df <- plot_df_total_N %>%
  group_by(scenario) %>%
  group_modify(
    ~ {
      h <- hist(
        .x$n_patients,
        breaks = seq(
          floor(min(plot_df_total_N$n_patients)) - 1.5,
          ceiling(max(plot_df_total_N$n_patients)) + 1.5,
          by = binwidth_use
        ),
        plot = FALSE
      )
      
      tibble(
        ymin  = h$breaks[-length(h$breaks)],
        ymax  = h$breaks[-1],
        count = h$counts
      )
    }
  ) %>%
  ungroup()

max_count <- max(hist_df$count)

hist_df <- hist_df %>%
  left_join(scenario_pos, by = "scenario") %>%
  mutate(
    width_scaled = count / max_count * 0.38,
    x_left  = scen_id,
    x_right = scen_id + width_scaled
  )

summary_df_total_N <- summary_df_total_N %>%
  left_join(scenario_pos, by = "scenario")

p_total_N_distribution <- ggplot() +
  
  geom_segment(
    data = scenario_pos,
    aes(
      x = scen_id,
      xend = scen_id,
      y = min(plot_df_total_N$n_patients) - 1,
      yend = max(plot_df_total_N$n_patients) + 1
    ),
    colour = "grey70",
    linewidth = 0.5
  ) +
  
  geom_rect(
    data = hist_df,
    aes(
      xmin = x_left,
      xmax = x_right,
      ymin = ymin,
      ymax = ymax,
      fill = scenario
    ),
    colour = "white",
    linewidth = 0.25,
    alpha = 0.85
  ) +
  
  geom_segment(
    data = summary_df_total_N,
    aes(
      x = scen_id + 0.42,
      xend = scen_id + 0.42,
      y = p10,
      yend = p90
    ),
    colour = "#10384F",
    linewidth = 1.2
  ) +
  
  geom_point(
    data = summary_df_total_N,
    aes(x = scen_id + 0.42, y = p10),
    colour = "#10384F",
    size = 2.3
  ) +
  
  geom_point(
    data = summary_df_total_N,
    aes(x = scen_id + 0.42, y = p90),
    colour = "#10384F",
    size = 2.3
  ) +
  
  geom_point(
    data = summary_df_total_N,
    aes(x = scen_id + 0.42, y = mean_n),
    colour = "red",
    size = 3
  ) +
  
  geom_point(
    data = summary_df_total_N,
    aes(x = scen_id + 0.42, y = median_n),
    colour = "black",
    size = 3,
    shape = 18
  ) +
  
  geom_text(
    data = summary_df_total_N,
    aes(
      x = scen_id + 0.50,
      y = mean_n,
      label = round(mean_n, 1)
    ),
    hjust = 0,
    size = 3.8
  ) +
  
  scale_x_continuous(
    breaks = scenario_pos$scen_id,
    labels = scenario_pos$scenario,
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  
  scale_fill_manual(
    values = bayer_pal,
    breaks = scenario_order,
    drop = FALSE
  ) +
  
  labs(
    title = "Distribution of total number of patients by scenario",
    x = "Scenario",
    y = "Total number of patients"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12)
  )

ggsave(
  filename = file.path(figure_dir, "total_N_distribution_by_scenario.png"),
  plot = p_total_N_distribution,
  width = 9,
  height = 5.5,
  dpi = 300
)


# =============================================================================
# 11. Figure: detected MTD distribution
# =============================================================================

mtd_df <- build_mtd_distribution_df(
  scenario_results = scenario_results,
  dose_grid = dose_grid,
  scenario_order = scenario_order
)

mtd_df_plot <- mtd_df %>%
  group_by(scenario) %>%
  mutate(is_mtd_mode = percent == max(percent)) %>%
  ungroup()

p_mtd_distribution <- ggplot(
  mtd_df_plot,
  aes(
    x = factor(dose),
    y = percent,
    fill = is_mtd_mode
  )
) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    vjust = -0.3,
    size = 3.8,
    color = "black"
  ) +
  facet_wrap(~ scenario, ncol = 4) +
  scale_fill_manual(
    values = c("FALSE" = "grey40", "TRUE" = "#007CBF"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "MTD selection distribution by true toxicity scenario",
    subtitle = "Blue bar indicates the most frequently selected MTD in each scenario",
    x = "Selected MTD dose (mg)",
    y = "Percent (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = file.path(figure_dir, "mtd_selection_distribution.png"),
  plot = p_mtd_distribution,
  width = 10,
  height = 6,
  dpi = 300
)


# =============================================================================
# 12. Figure: probability of selecting a dose above the true MTD
# =============================================================================

true_dlt_tbl <- purrr::imap_dfr(
  true_scenarios,
  ~ tibble(
    scenario = .y,
    dose = dose_grid,
    true_dlt = .x
  )
)

true_mtd_tbl <- true_dlt_tbl %>%
  group_by(scenario) %>%
  slice_min(
    order_by = abs(true_dlt - tox_thr),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    scenario,
    true_mtd_dose = dose,
    true_mtd_true_dlt = true_dlt
  )

pod_tbl <- compute_pod_selected_above_true_mtd(
  mtd_df = mtd_df,
  true_mtd_tbl = true_mtd_tbl
) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = scenario_order,
      ordered = TRUE
    )
  )

p_pod <- ggplot(
  pod_tbl,
  aes(
    x = scenario,
    y = PoD
  )
) +
  geom_col(fill = "#4C78A8", width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f", PoD)),
    vjust = -0.35,
    size = 4.2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Scenario",
    y = "PoD (%)",
    title = "Pr(selected MTD > true MTD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(
  filename = file.path(figure_dir, "probability_selected_mtd_above_true_mtd.png"),
  plot = p_pod,
  width = 8,
  height = 5,
  dpi = 300
)


# =============================================================================
# 13. Optional individual trial plot
# =============================================================================
#
# This section is intentionally disabled by default because it is exploratory and
# depends on a manually selected scenario and simulation index.
#
# To inspect a single trial, uncomment and adapt the lines below.
#
# selected_scenario <- "Third"
# selected_trial <- 17
#
# trial_data <- scenario_results[[selected_scenario]]@data[[selected_trial]]
# plot(trial_data, mark_backfill = FALSE)


# =============================================================================
# 14. Session information
# =============================================================================

session_info_path <- file.path(output_dir, "sessionInfo.txt")

capture.output(
  sessionInfo(),
  file = session_info_path
)

message("Session information exported: ", session_info_path)

message("Script completed successfully.")