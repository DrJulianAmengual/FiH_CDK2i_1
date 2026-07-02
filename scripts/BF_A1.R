###############################################################
# Posterior predictive simulations across several sample sizes
# PROC backfill setting: ORR + chronic toxicity + durability
###############################################################
#
# Purpose:
#   Evaluate the operating characteristics of a GO / MODIFY / NO-GO
#   decision framework using Monte Carlo simulations.
#
# Endpoints:
#   1. ORR: binary efficacy endpoint
#   2. Chronic toxicity: binary safety endpoint
#   3. Durability success: binary durability endpoint
#
# Important:
#   This is not a time-to-event analysis.
#   There is no PFS6, OS12, Kaplan-Meier, censoring, or Weibull modelling.
#   Each endpoint is treated as a binomial proportion.
#
# Main output:
#   - A slide-ready table summarising decision probabilities by N
#   - A figure showing how GO / MODIFY / NO-GO probabilities evolve with N
#
###############################################################


# ============================================================
# 0) Load required packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
})


# ============================================================
# 1) Simulation settings
# ============================================================

set.seed(123)

# Number of simulated trials per true scenario and per N.
# A larger value gives more stable Monte Carlo estimates.
nsim <- 10000

# Sample sizes to evaluate.
# These are the possible backfilling cohort sizes.
N_grid <- c(10, 15, 20, 25, 30,35,40,45,50)


# ============================================================
# 2) Bayesian prior
# ============================================================
#
# We use Jeffreys prior for each binomial endpoint:
#
#   pi ~ Beta(0.5, 0.5)
#
# If r events are observed among N patients:
#
#   pi | data ~ Beta(0.5 + r, 0.5 + N - r)
#
# The same prior is used for ORR, toxicity, and durability.
# ============================================================

a0 <- 0.5
b0 <- 0.5


# ============================================================
# 3) Decision-rule thresholds
# ============================================================

# -----------------------------
# ORR thresholds
# -----------------------------
#
# pi_go:
#   Clinically relevant ORR threshold for GO.
#
# pi_mod:
#   Lower ORR threshold used for the NO-GO / MODIFY boundary.
#
# p_go:
#   Required posterior probability for GO:
#     Pr(ORR >= pi_go | data) > p_go
#
# p_nogo:
#   NO-GO rule for low activity:
#     Pr(ORR >= pi_mod | data) < p_nogo
#
# p_modL and p_modU:
#   Intermediate range interpreted as MODIFY.
# -----------------------------

pi_go  <- 0.30
pi_mod <- 0.25

p_go   <- 0.80
p_nogo <- 0.20

p_modL <- 0.50
p_modU <- 0.80


# -----------------------------
# Chronic toxicity thresholds
# -----------------------------
#
# tox_thr:
#   Chronic toxicity threshold considered clinically unacceptable.
#
# p_tox_bad:
#   If Pr(toxicity >= tox_thr | data) > p_tox_bad,
#   the decision is NO-GO due to safety.
#
# p_tox_ok:
#   For a GO decision, we require:
#     Pr(toxicity < tox_thr | data) > p_tox_ok
# -----------------------------

tox_thr <- 0.35

p_tox_bad <- 0.70
p_tox_ok  <- 0.85


# -----------------------------
# Durability thresholds
# -----------------------------
#
# dur_thr:
#   Clinically meaningful durability success rate.
#
# p_dur_bad:
#   If Pr(durability >= dur_thr | data) < p_dur_bad,
#   the decision is NO-GO due to poor durability.
#
# p_dur_ok:
#   For a GO decision, we require:
#     Pr(durability >= dur_thr | data) > p_dur_ok
# -----------------------------

dur_thr <- 0.70

p_dur_bad <- 0.20
p_dur_ok  <- 0.80


# ============================================================
# 4) True simulation scenarios
# ============================================================
#
# These values define the true data-generating probabilities.
#
# Each simulated scenario is one combination of:
#   - true ORR
#   - true chronic toxicity rate
#   - true durability success rate
#
# You can modify these values depending on the scenarios you want
# to show in the PROC/backfill discussion.
# ============================================================

true_orr_scenarios <- c(0.20, 0.25, 0.30, 0.40)
true_tox_scenarios <- c(0.10, 0.20, 0.30, 0.40)
true_dur_scenarios <- c(0.60, 0.80, 0.90)


# ============================================================
# 5) Posterior probability helper functions
# ============================================================

# Pr(pi >= threshold | data)
post_prob_ge <- function(thr, r, n) {
  pbeta(
    q = thr,
    shape1 = a0 + r,
    shape2 = b0 + n - r,
    lower.tail = FALSE
  )
}

# Pr(pi < threshold | data)
post_prob_lt <- function(thr, r, n) {
  pbeta(
    q = thr,
    shape1 = a0 + r,
    shape2 = b0 + n - r,
    lower.tail = TRUE
  )
}


# ============================================================
# 6) ORR-only decision rule
# ============================================================
#
# This function applies the efficacy rule only.
#
# GO:
#   Pr(ORR >= 30%) > 80%
#
# NO-GO:
#   Pr(ORR >= 25%) < 20%
#
# MODIFY:
#   Intermediate / borderline cases.
#
# Remaining cases are assigned to MODIFY by design, so we do not
# create a separate "indeterminate" category.
# ============================================================

decide_orr_only <- function(pr_ge_30, pr_ge_25) {
  
  case_when(
    pr_ge_30 > p_go ~ "GO",
    pr_ge_25 < p_nogo ~ "NO-GO",
    pr_ge_25 >= p_modL & pr_ge_25 <= p_modU ~ "MODIFY",
    TRUE ~ "MODIFY"
  )
}


# ============================================================
# 7) Joint decision rule
# ============================================================
#
# Final decision uses ORR, chronic toxicity, and durability.
#
# Conservative logic:
#
#   1. NO-GO if efficacy is clearly insufficient.
#   2. NO-GO if chronic toxicity is unacceptable.
#   3. NO-GO if durability is clearly poor.
#   4. GO only if efficacy, safety, and durability are all supportive.
#   5. MODIFY otherwise.
#
# This rule avoids declaring GO if one key dimension is not supportive.
# ============================================================

decide_joint <- function(pr_orr_ge_30,
                         pr_orr_ge_25,
                         pr_tox_ok,
                         pr_tox_bad,
                         pr_dur_ok) {
  
  case_when(
    
    # Lack of sufficient efficacy
    pr_orr_ge_25 < p_nogo ~ "NO-GO",
    
    # Unacceptable chronic toxicity
    pr_tox_bad > p_tox_bad ~ "NO-GO",
    
    # Poor durability
    pr_dur_ok < p_dur_bad ~ "NO-GO",
    
    # GO only if all components are supportive
    pr_orr_ge_30 > p_go &
      pr_tox_ok  > p_tox_ok &
      pr_dur_ok  > p_dur_ok ~ "GO",
    
    # Everything else is MODIFY
    TRUE ~ "MODIFY"
  )
}

# ============================================================
# 8) Simulation function for one N
# ============================================================
#
# For a fixed sample size N, this function:
#
#   1. Builds the full grid of true scenarios.
#   2. Simulates nsim trials for each scenario.
#   3. Generates binomial counts:
#        r_orr = number of responders
#        r_tox = number with chronic toxicity
#        r_dur = number with durability success
#   4. Computes posterior probabilities.
#   5. Applies the joint GO / MODIFY / NO-GO decision rule.
#   6. Returns decision probabilities for each true scenario.
#
# Important:
#   ORR, toxicity, and durability are simulated independently.
#   This is a marginal operating-characteristics assessment.
# ============================================================

run_simulation_for_N <- function(N, nsim = 10000) {
  
  # Create all combinations of true ORR, true toxicity, and true durability
  scenario_grid <- expand_grid(
    true_orr = true_orr_scenarios,
    true_tox = true_tox_scenarios,
    true_dur = true_dur_scenarios
  )
  
  # Simulate each scenario
  sim_counts <- pmap_dfr(
    scenario_grid,
    function(true_orr, true_tox, true_dur) {
      
      # ------------------------------------------------------
      # Simulate binomial data for nsim independent trials
      # ------------------------------------------------------
      
      r_orr <- rbinom(nsim, size = N, prob = true_orr)
      r_tox <- rbinom(nsim, size = N, prob = true_tox)
      r_dur <- rbinom(nsim, size = N, prob = true_dur)
      
      
      # ------------------------------------------------------
      # Compute posterior probabilities
      # ------------------------------------------------------
      
      # Efficacy posterior probabilities
      pr_orr_ge_30 <- post_prob_ge(pi_go,  r_orr, N)
      pr_orr_ge_25 <- post_prob_ge(pi_mod, r_orr, N)
      
      # Safety posterior probabilities
      pr_tox_ok_sim  <- post_prob_lt(tox_thr, r_tox, N)
      pr_tox_bad_sim <- post_prob_ge(tox_thr, r_tox, N)
      
      # Durability posterior probability
      pr_dur_ok_sim <- post_prob_ge(dur_thr, r_dur, N)
      
      
      # ------------------------------------------------------
      # Apply final joint decision rule
      # ------------------------------------------------------
      
      decision <- decide_joint(
        pr_orr_ge_30 = pr_orr_ge_30,
        pr_orr_ge_25 = pr_orr_ge_25,
        pr_tox_ok    = pr_tox_ok_sim,
        pr_tox_bad   = pr_tox_bad_sim,
        pr_dur_ok    = pr_dur_ok_sim
      )
      
      
      # ------------------------------------------------------
      # Count decisions across simulated trials
      # ------------------------------------------------------
      
      tibble(
        N        = N,
        true_orr = true_orr,
        true_tox = true_tox,
        true_dur = true_dur,
        decision = decision
      ) %>%
        count(
          N,
          true_orr,
          true_tox,
          true_dur,
          decision,
          name = "n"
        )
    }
  )
  
  
  # ----------------------------------------------------------
  # Convert counts into estimated decision probabilities
  # ----------------------------------------------------------
  
  sim_counts %>%
    group_by(N, true_orr, true_tox, true_dur) %>%
    summarise(
      GO = sum(n[decision == "GO"], na.rm = TRUE) / sum(n),
      MODIFY = sum(n[decision == "MODIFY"], na.rm = TRUE) / sum(n),
      NO_GO = sum(n[decision == "NO-GO"], na.rm = TRUE) / sum(n),
      .groups = "drop"
    )
}

# ============================================================
# 9) Run simulations for all selected sample sizes
# ============================================================

results_wide <- map_dfr(
  N_grid,
  ~ run_simulation_for_N(N = .x, nsim = nsim)
)

results_wide

# ============================================================
# 10) Define scenario groups
# ============================================================
#
# Favorable:
#   Scenarios where the true profile is clinically attractive:
#     - ORR >= 30%
#     - chronic toxicity <= 30%
#     - durability success >= 80%
#
# Unfavorable:
#   Scenarios where at least one key dimension is clearly poor:
#     - ORR <= 25%
#     - or chronic toxicity > 30%
#     - or durability success <= 60%
#
# Other:
#   Intermediate scenarios. These are excluded from the main slide summary
#   to keep the interpretation focused.
# ============================================================

results_grouped <- results_wide %>%
  mutate(
    scenario_group = case_when(
      true_orr >= 0.30 &
        true_tox <= 0.30 &
        true_dur >= 0.80 ~ "Favorable",
      
      true_orr <= 0.25 |
        true_tox > 0.30 |
        true_dur <= 0.60 ~ "Unfavorable",
      
      TRUE ~ "Other"
    )
  ) %>%
  filter(
    scenario_group %in% c("Favorable", "Unfavorable")
  )

results_grouped

# ============================================================
# 11) Table: decision probabilities by N and scenario group
# ============================================================
#
# For each N and each scenario group, we summarise the decision
# probabilities across all scenario combinations in that group.
#
# Median is used because each group contains multiple scenario combinations.
# Ranges are included to show variability across scenarios.
# ============================================================

summary_by_group <- results_grouped %>%
  group_by(scenario_group, N) %>%
  summarise(
    n_scenarios = n(),
    
    GO_median = median(GO, na.rm = TRUE),
    MODIFY_median = median(MODIFY, na.rm = TRUE),
    NO_GO_median = median(NO_GO, na.rm = TRUE),
    
    GO_min = min(GO, na.rm = TRUE),
    GO_max = max(GO, na.rm = TRUE),
    
    MODIFY_min = min(MODIFY, na.rm = TRUE),
    MODIFY_max = max(MODIFY, na.rm = TRUE),
    
    NO_GO_min = min(NO_GO, na.rm = TRUE),
    NO_GO_max = max(NO_GO, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    scenario_group = factor(
      scenario_group,
      levels = c("Favorable", "Unfavorable")
    ),
    
    `GO median (%)` = round(100 * GO_median, 1),
    `MODIFY median (%)` = round(100 * MODIFY_median, 1),
    `NO-GO median (%)` = round(100 * NO_GO_median, 1),
    
    `GO range (%)` = paste0(
      round(100 * GO_min, 1),
      "–",
      round(100 * GO_max, 1)
    ),
    
    `MODIFY range (%)` = paste0(
      round(100 * MODIFY_min, 1),
      "–",
      round(100 * MODIFY_max, 1)
    ),
    
    `NO-GO range (%)` = paste0(
      round(100 * NO_GO_min, 1),
      "–",
      round(100 * NO_GO_max, 1)
    )
  ) %>%
  arrange(scenario_group, N)

summary_table_slide <- summary_by_group %>%
  select(
    `Scenario group` = scenario_group,
    N,
    `Number of scenarios` = n_scenarios,
    `GO median (%)`,
    `MODIFY median (%)`,
    `NO-GO median (%)`,
    `GO range (%)`,
    `MODIFY range (%)`,
    `NO-GO range (%)`
  )

summary_table_slide

# ============================================================
# 12) Prepare data for figure
# ============================================================

plot_data <- summary_by_group %>%
  select(
    scenario_group,
    N,
    GO_median,
    MODIFY_median,
    NO_GO_median
  ) %>%
  pivot_longer(
    cols = c(GO_median, MODIFY_median, NO_GO_median),
    names_to = "decision",
    values_to = "prob"
  ) %>%
  mutate(
    decision = recode(
      decision,
      "GO_median" = "GO",
      "MODIFY_median" = "MODIFY",
      "NO_GO_median" = "NO-GO"
    ),
    decision = factor(decision, levels = c("GO", "MODIFY", "NO-GO")),
    prob_pct = 100 * prob
  )


# ============================================================
# 13) Figure: decision probability across N
# ============================================================

decision_palette <- c(
  "GO" = "#009E73",
  "MODIFY" = "#E69F00",
  "NO-GO" = "#D30F4B"
)

p_decision_by_N <- ggplot(
  plot_data,
  aes(
    x = N,
    y = prob_pct,
    color = decision,
    group = decision
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  facet_wrap(~ scenario_group, nrow = 1) +
  scale_color_manual(values = decision_palette) +
  scale_x_continuous(breaks = N_grid) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 100)
  ) +
  labs(
    title = "Decision probabilities across backfilling sample sizes",
    subtitle = "Median across simulated true scenarios within each scenario group",
    x = "Backfilling sample size (N)",
    y = "Decision probability (%)",
    color = "Decision"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_decision_by_N
