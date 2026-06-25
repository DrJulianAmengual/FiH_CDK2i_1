###############################################################
# Posterior predictive simulations across several sample sizes
# A2/B2 setting: ORR + CBR24 + chronic toxicity + durability
###############################################################
#
# Purpose:
#   Evaluate the operating characteristics of a GO / MODIFY / NO-GO
#   decision framework using Monte Carlo simulations.
#
# Endpoints:
#   1. ORR:
#        Binary response endpoint.
#
#   2. CBR24:
#        Binary clinical benefit endpoint at 24 weeks.
#
#   3. Chronic toxicity:
#        Binary G3+ chronic toxicity endpoint.
#
#   4. Durability:
#        Binary durability success endpoint.
#
# Important:
#   This is not a time-to-event analysis.
#   There is no PFS6, OS12, Kaplan-Meier, censoring, or Weibull modelling.
#   All endpoints are treated as binary binomial outcomes.
#
# Main outputs:
#   - results_allN_A2B2:
#       Scenario-level decision probabilities for each N.
#
#   - summary_A2B2_table_pct:
#       Slide-ready table by scenario group and sample size.
#
#   - p_decision_by_N_A2B2:
#       Figure showing how GO / MODIFY / NO-GO probabilities change with N.
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

# Number of simulated trials per true scenario and per sample size.
# Increase this value for more stable Monte Carlo estimates.
nsim <- 10000

# Sample sizes to evaluate.
# These are the possible backfilling cohort sizes.
N_values <- c(10, 15, 20, 25, 30)


# ============================================================
# 2) Bayesian prior
# ============================================================
#
# Jeffreys prior for each binary endpoint:
#
#   pi ~ Beta(0.5, 0.5)
#
# If r events are observed among N patients:
#
#   pi | data ~ Beta(0.5 + r, 0.5 + N - r)
#
# The same prior is used for ORR, CBR24, toxicity, and durability.
# ============================================================

a0 <- 0.5
b0 <- 0.5


# ============================================================
# 3) Decision-rule thresholds
# ============================================================

# ------------------------------------------------------------
# Efficacy thresholds
# ------------------------------------------------------------
#
# ORR criterion:
#   Pr(ORR >= orr_thr | data)
#
# CBR24 criterion:
#   Pr(CBR24 >= cbr_thr | data)
#
# Efficacy is considered strong if either ORR or CBR24 is strong.
# ------------------------------------------------------------

orr_thr <- 0.25      # ORR threshold: ORR >= 25%
cbr_thr <- 0.40      # CBR24 threshold: CBR24 >= 40%

p_eff_go   <- 0.80   # Posterior probability required for efficacy GO
p_eff_modL <- 0.60   # Lower bound for borderline efficacy
p_eff_modU <- 0.80   # Upper bound for borderline efficacy
p_eff_nogo <- 0.35   # Futility cutoff for low efficacy


# ------------------------------------------------------------
# Chronic toxicity thresholds
# ------------------------------------------------------------
#
# Safety has two components:
#
#   1. Safety acceptability:
#        Pr(toxicity < 35% | data) > 80%
#
#   2. Toxicity concern:
#        Pr(toxicity >= 30% | data) > 70%
#
# The second criterion is a conservative safety override.
# ------------------------------------------------------------

tox_thr <- 0.35       # Acceptability threshold: tox < 35%
p_tox_ok <- 0.80      # Required posterior probability for acceptable tox

tox_thr_bad <- 0.30   # Toxicity concern threshold: tox >= 30%
p_tox_bad <- 0.70     # Posterior probability triggering safety NO-GO


# ------------------------------------------------------------
# Durability thresholds
# ------------------------------------------------------------
#
# Durability success criterion:
#   Pr(durability >= 70% | data)
#
# Durability is required for GO.
# Very poor durability can trigger NO-GO.
# ------------------------------------------------------------

dur_thr <- 0.70

p_dur_ok <- 0.80        # Required posterior probability for durability support
p_dur_futile <- 0.20    # Futility cutoff for poor durability


# ============================================================
# 4) True simulation scenarios
# ============================================================
#
# These are the true data-generating values used in the simulations.
#
# Each true scenario is one combination of:
#   - true ORR
#   - true CBR24
#   - true chronic toxicity
#   - true durability success
#
# You can edit these values depending on what you want to show
# in the PROC/backfilling discussion.
# ============================================================

true_orr_scenarios <- c(0.20, 0.25, 0.30, 0.40)
true_cbr_scenarios <- c(0.30, 0.40, 0.50)
true_tox_scenarios <- c(0.10, 0.20, 0.30, 0.40)
true_dur_scenarios <- c(0.60, 0.80, 0.90)


# Create full scenario grid.
scenario_grid <- expand_grid(
  true_orr = true_orr_scenarios,
  true_cbr = true_cbr_scenarios,
  true_tox = true_tox_scenarios,
  true_dur = true_dur_scenarios
)


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
# 6) Joint A2/B2 decision rule
# ============================================================
#
# Inputs:
#   pr_orr:
#     Posterior probability Pr(ORR >= orr_thr | data)
#
#   pr_cbr:
#     Posterior probability Pr(CBR24 >= cbr_thr | data)
#
#   pr_dur_ok:
#     Posterior probability Pr(durability >= dur_thr | data)
#
#   pr_tox_ok:
#     Posterior probability Pr(toxicity < tox_thr | data)
#
#   pr_tox_bad:
#     Posterior probability Pr(toxicity >= tox_thr_bad | data)
#
# Decision logic:
#
#   1. Safety-first:
#        NO-GO if Pr(toxicity >= 30%) > 70%.
#        NO-GO if Pr(toxicity < 35%) <= 80%.
#
#      Note:
#        The second rule is very conservative because failure to demonstrate
#        acceptable safety leads to NO-GO rather than MODIFY.
#
#   2. Efficacy futility:
#        NO-GO if both ORR and CBR24 posterior probabilities are low.
#
#   3. Durability futility:
#        NO-GO if durability support is very weak.
#
#   4. GO:
#        GO if efficacy is strong and durability is supportive.
#
#   5. MODIFY:
#        MODIFY for efficacy GO with insufficient durability,
#        or for borderline efficacy,
#        or for all remaining intermediate cases.
# ============================================================

decide_A2B2 <- function(pr_orr,
                        pr_cbr,
                        pr_dur_ok,
                        pr_tox_ok,
                        pr_tox_bad) {
  
  # Strong efficacy if either ORR or CBR24 is strong.
  eff_go <- (pr_cbr >= p_eff_go) | (pr_orr >= p_eff_go)
  
  # Borderline efficacy if either ORR or CBR24 is in the MODIFY range.
  eff_mod <- (pr_orr >= p_eff_modL & pr_orr <= p_eff_modU) |
    (pr_cbr >= p_eff_modL & pr_cbr <= p_eff_modU)
  
  # Efficacy futility:
  # both ORR and CBR24 must be low to declare low efficacy.
  eff_low <- (pr_orr < p_eff_nogo) & (pr_cbr < p_eff_nogo)
  
  case_when(
    
    # Safety override: high probability of unacceptable chronic toxicity
    pr_tox_bad > p_tox_bad ~ "NO-GO",
    
    # Conservative safety rule:
    # if we do not have enough evidence that toxicity is below 35%,
    # classify as NO-GO.
    pr_tox_ok <= p_tox_ok ~ "NO-GO",
    
    # Efficacy futility
    eff_low ~ "NO-GO",
    
    # Durability futility
    pr_dur_ok < p_dur_futile ~ "NO-GO",
    
    # Strong efficacy but durability not strong enough
    eff_go & pr_dur_ok < p_dur_ok ~ "MODIFY",
    
    # GO only if efficacy and durability are both supportive
    eff_go & pr_dur_ok >= p_dur_ok ~ "GO",
    
    # Borderline efficacy
    eff_mod ~ "MODIFY",
    
    # Remaining intermediate cases
    TRUE ~ "MODIFY"
  )
}


# ============================================================
# 7) Simulation function for one sample size N
# ============================================================
#
# For a fixed N:
#
#   1. Simulate nsim trials for every true scenario.
#   2. Generate binomial counts for:
#        - ORR
#        - CBR24
#        - chronic toxicity
#        - durability success
#   3. Compute posterior probabilities.
#   4. Apply the A2/B2 decision rule.
#   5. Return decision probabilities by true scenario.
#
# Important:
#   This simulation treats ORR, CBR24, toxicity, and durability as
#   marginally independent binary endpoints.
#
#   If a patient-level correlation structure is required, the data-generating
#   mechanism would need to be modified.
# ============================================================

run_simulation_for_N_A2B2 <- function(N, nsim = 10000) {
  
  sim_counts <- pmap_dfr(
    scenario_grid,
    function(true_orr, true_cbr, true_tox, true_dur) {
      
      # ------------------------------------------------------
      # Simulate endpoint counts for nsim independent trials
      # ------------------------------------------------------
      
      r_orr <- rbinom(nsim, size = N, prob = true_orr)
      r_cbr <- rbinom(nsim, size = N, prob = true_cbr)
      r_tox <- rbinom(nsim, size = N, prob = true_tox)
      r_dur <- rbinom(nsim, size = N, prob = true_dur)
      
      
      # ------------------------------------------------------
      # Posterior probabilities
      # ------------------------------------------------------
      
      pr_orr <- post_prob_ge(
        thr = orr_thr,
        r = r_orr,
        n = N
      )
      
      pr_cbr <- post_prob_ge(
        thr = cbr_thr,
        r = r_cbr,
        n = N
      )
      
      pr_tox_ok <- post_prob_lt(
        thr = tox_thr,
        r = r_tox,
        n = N
      )
      
      pr_tox_bad <- post_prob_ge(
        thr = tox_thr_bad,
        r = r_tox,
        n = N
      )
      
      pr_dur_ok <- post_prob_ge(
        thr = dur_thr,
        r = r_dur,
        n = N
      )
      
      
      # ------------------------------------------------------
      # Apply decision rule
      # ------------------------------------------------------
      
      decision <- decide_A2B2(
        pr_orr = pr_orr,
        pr_cbr = pr_cbr,
        pr_dur_ok = pr_dur_ok,
        pr_tox_ok = pr_tox_ok,
        pr_tox_bad = pr_tox_bad
      )
      
      
      # ------------------------------------------------------
      # Count decisions across simulated trials
      # ------------------------------------------------------
      
      tibble(
        true_orr = true_orr,
        true_cbr = true_cbr,
        true_tox = true_tox,
        true_dur = true_dur,
        decision = decision
      ) %>%
        count(
          true_orr,
          true_cbr,
          true_tox,
          true_dur,
          decision,
          name = "n"
        )
    }
  )
  
  
  # ----------------------------------------------------------
  # Convert counts into decision probabilities
  # ----------------------------------------------------------
  
  sim_counts %>%
    group_by(true_orr, true_cbr, true_tox, true_dur) %>%
    summarise(
      GO = sum(n[decision == "GO"], na.rm = TRUE) / sum(n),
      MODIFY = sum(n[decision == "MODIFY"], na.rm = TRUE) / sum(n),
      NO_GO = sum(n[decision == "NO-GO"], na.rm = TRUE) / sum(n),
      .groups = "drop"
    ) %>%
    mutate(N = N)
}


# ============================================================
# 8) Run simulations for all sample sizes
# ============================================================

results_allN_A2B2 <- map_dfr(
  N_values,
  ~ run_simulation_for_N_A2B2(N = .x, nsim = nsim)
)

results_allN_A2B2


# ============================================================
# 9) Define scenario groups
# ============================================================
#
# For slide-level interpretation, we collapse the full scenario grid
# into favourable and unfavourable groups.
#
# Favorable:
#   - true ORR >= 30%
#   - true CBR24 >= 40%
#   - true chronic toxicity <= 20%
#   - true durability success >= 80%
#
# Unfavorable:
#   - true ORR <= 25%, or
#   - true CBR24 <= 30%, or
#   - true chronic toxicity >= 30%, or
#   - true durability success <= 60%
#
# Other:
#   Intermediate scenarios, excluded from the main summary table and plot.
# ============================================================

results_grouped_A2B2 <- results_allN_A2B2 %>%
  mutate(
    scenario_group = case_when(
      true_orr >= 0.30 &
        true_cbr >= 0.40 &
        true_tox <= 0.20 &
        true_dur >= 0.80 ~ "Favorable",
      
      true_orr <= 0.25 |
        true_cbr <= 0.30 |
        true_tox >= 0.30 |
        true_dur <= 0.60 ~ "Unfavorable",
      
      TRUE ~ "Other"
    )
  ) %>%
  filter(
    scenario_group %in% c("Favorable", "Unfavorable")
  )


# ============================================================
# 10) Summary table by scenario group and N
# ============================================================
#
# Each group contains multiple true scenario combinations.
#
# We summarise decision probabilities using the median across
# scenario combinations in each group.
#
# Ranges are included to show heterogeneity within each group.
# ============================================================

summary_A2B2_by_group <- results_grouped_A2B2 %>%
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


# ------------------------------------------------------------
# Slide-ready percentage table
# ------------------------------------------------------------

summary_A2B2_table_pct <- summary_A2B2_by_group %>%
  select(
    `Scenario group` = scenario_group,
    `Sample size` = N,
    `Number of scenarios` = n_scenarios,
    `GO median (%)`,
    `MODIFY median (%)`,
    `NO-GO median (%)`,
    `GO range (%)`,
    `MODIFY range (%)`,
    `NO-GO range (%)`
  )

summary_A2B2_table_pct


# ------------------------------------------------------------
# Compact numeric table, if preferred
# ------------------------------------------------------------

summary_A2B2_table_numeric <- summary_A2B2_by_group %>%
  transmute(
    Scenario = scenario_group,
    `Sample size` = N,
    `P(GO)` = round(GO_median, 2),
    `P(MODIFY)` = round(MODIFY_median, 2),
    `P(NO-GO)` = round(NO_GO_median, 2)
  )

summary_A2B2_table_numeric


# ============================================================
# 11) Prepare data for plotting
# ============================================================

plot_data_A2B2 <- summary_A2B2_by_group %>%
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
# 12) Figure: decision probabilities across N
# ============================================================
#
# This figure shows how decision probabilities change across
# increasing backfilling sample size.
#
# Interpretation:
#
#   - In favourable scenarios:
#       desirable behaviour is higher GO probability.
#
#   - In unfavourable scenarios:
#       desirable behaviour is higher NO-GO probability.
#
#   - MODIFY captures intermediate or insufficiently conclusive cases.
# ============================================================

decision_palette <- c(
  "GO" = "#009E73",
  "MODIFY" = "#E69F00",
  "NO-GO" = "#D30F4B"
)

p_decision_by_N_A2B2 <- ggplot(
  plot_data_A2B2,
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
  scale_x_continuous(breaks = N_values) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 100)
  ) +
  labs(
    title = "Decision probabilities across backfilling sample sizes",
    subtitle = paste0(
      "A2/B2 framework: ORR, CBR24, chronic toxicity and durability; ",
      "median across scenario combinations"
    ),
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

p_decision_by_N_A2B2


# ============================================================
# 13) Alternative figure: one panel per decision
# ============================================================
#
# This version is often clearer for slides because it separates
# GO, MODIFY, and NO-GO into different panels.
# ============================================================

p_decision_by_N_A2B2_faceted <- ggplot(
  plot_data_A2B2,
  aes(
    x = N,
    y = prob_pct,
    color = scenario_group,
    group = scenario_group
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  facet_wrap(~ decision, nrow = 1) +
  scale_x_continuous(breaks = N_values) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 100)
  ) +
  labs(
    title = "Decision probability by sample size and scenario group",
    subtitle = "Median probability across simulated true scenarios",
    x = "Backfilling sample size (N)",
    y = "Decision probability (%)",
    color = "Scenario group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_decision_by_N_A2B2_faceted