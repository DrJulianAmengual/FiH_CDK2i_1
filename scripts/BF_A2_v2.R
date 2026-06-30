# ============================================================
# Bayesian operating characteristics simulation
# ABC Dose Escalation Back-Fills
#
# Features:
# - ORR simulated among measurable patients
# - CBR24 simulated coherently as ORR + prolonged SD
# - DOR8 simulated as binomial among responders only
# - Chronic G3+ toxicity simulated as binomial among treated patients
# - GO / MODIFY / NO-GO decision rule aligned with slide criteria
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

set.seed(123)

# ============================================================
# PARAMETERS
# ============================================================

N    <- 15
N_meas <- N
nsim <- 10000

# ------------------------------------------------------------
# Priors
# ------------------------------------------------------------
# Jeffreys prior: weakly informative, symmetric, commonly used
# for binomial endpoints when we do not want to strongly favour
# either success or failure.

a_orr <- 0.5
b_orr <- 0.5

a_cbr <- 0.5
b_cbr <- 0.5

#a_tox <- 0.5
#b_tox <- 0.5

a_dor <- 0.5
b_dor <- 0.5

# Optional more conservative toxicity prior:
# centred at 35% with ESS = 4
# a_tox <- 1.4
# b_tox <- 2.6


# ============================================================
# DECISION THRESHOLDS
# ============================================================

# ------------------------------------------------------------
# Efficacy
# ------------------------------------------------------------

cbr_thr <- 0.40      # CBR at 24 weeks >= 40%

# The slide says ORR >= 25-30%.
# For simulation, one value should be pre-specified.
# I use 0.275 as base case.
# You can run sensitivity with 0.25 and 0.30.
orr_thr <- 0.275

p_eff_go <- 0.70    # Pr(CBR >= 40%) > 80% OR Pr(ORR >= theta) > 80%

# Futility threshold for efficacy
p_eff_nogo <- 0.10


# ------------------------------------------------------------
# DOR
# ------------------------------------------------------------

dor8_thr <- 0.70     # DOR at 8 months > 70%

# Since no censoring is included, DOR8 is binary:
# among responders, whether DOR >= 8 months.
p_dor_go <- 0.70     # Pr(p_DOR8 >= 70%) > 70%

# Futility threshold for DOR
p_dor_nogo <- 0.20


# ------------------------------------------------------------
# Chronic toxicity
# ------------------------------------------------------------

tox_thr <- 0.35      # chronic Grade 3+ toxicity < 35%

p_tox_go <- 0.80     # Pr(p_tox < 35%) > 80%

# Safety NO-GO:
# If probability that tox is above/equal 35% is high, stop/no-go.
p_tox_bad <- 0.80    # Pr(p_tox >= 35%) > 80% => NO-GO


# ============================================================
# SCENARIO GRID
# ============================================================

true_orr_scenarios <- c(0.20, 0.25, 0.30, 0.40)
true_cbr_scenarios <- c(0.30, 0.40, 0.50, 0.60)
true_tox_scenarios <- c(0.10, 0.20, 0.30, 0.40)
true_dor8_scenarios <- c(0.60, 0.70, 0.80, 0.90)

grid <- expand_grid(
  true_orr  = true_orr_scenarios,
  true_cbr  = true_cbr_scenarios,
  true_tox  = true_tox_scenarios,
  true_dor8 = true_dor8_scenarios
) %>%
  # CBR includes ORR, so scenarios with CBR < ORR are incoherent.
  filter(true_cbr >= true_orr)


# ============================================================
# POSTERIOR HELPERS
# ============================================================

post_prob_ge_beta <- function(thr, r, n, a0, b0) {
  pbeta(
    q = thr,
    shape1 = a0 + r,
    shape2 = b0 + n - r,
    lower.tail = FALSE
  )
}

post_prob_lt_beta <- function(thr, r, n, a0, b0) {
  pbeta(
    q = thr,
    shape1 = a0 + r,
    shape2 = b0 + n - r,
    lower.tail = TRUE
  )
}




# ============================================================
# SIMULATE ORR AND CBR COHERENTLY
# ============================================================

simulate_orr_cbr <- function(nsim, N_meas, true_orr, true_cbr) {
  
  # ----------------------------------------------------------
  # Step 1: simulate objective responses
  # ----------------------------------------------------------
  # ORR patients are automatically part of CBR.
  
  r_orr <- rbinom(
    n = nsim,
    size = N_meas,
    prob = true_orr
  )
  
  # ----------------------------------------------------------
  # Step 2: simulate prolonged stable disease among
  # non-responders only.
  # ----------------------------------------------------------
  #
  # We want:
  #
  #   P(CBR) = P(ORR) + P(no ORR) * P(prolonged SD | no ORR)
  #
  # Therefore:
  #
  #   P(prolonged SD | no ORR)
  #   = (true_cbr - true_orr) / (1 - true_orr)
  #
  # This guarantees that CBR >= ORR at patient level.
  
  p_sd_prolonged_given_no_orr <- ifelse(
    true_orr < 1,
    (true_cbr - true_orr) / (1 - true_orr),
    0
  )
  
  p_sd_prolonged_given_no_orr <- pmin(
    pmax(p_sd_prolonged_given_no_orr, 0),
    1
  )
  
  n_nonresponders <- N_meas - r_orr
  
  r_sd_prolonged <- rbinom(
    n = nsim,
    size = n_nonresponders,
    prob = p_sd_prolonged_given_no_orr
  )
  
  r_cbr <- r_orr + r_sd_prolonged
  
  tibble(
    r_orr = r_orr,
    r_sd_prolonged = r_sd_prolonged,
    r_cbr = r_cbr
  )
}


# ============================================================
# SIMULATE DOR8 AMONG RESPONDERS
# ============================================================

simulate_dor8 <- function(r_orr, true_dor8) {
  
  # DOR8 is only defined among responders.
  # If r_orr = 0, then there are no responders contributing to DOR8.
  
  rbinom(
    n = length(r_orr),
    size = r_orr,
    prob = true_dor8
  )
}


# ============================================================
# DECISION FUNCTION
# ============================================================

decide_A2B2 <- function(pr_orr,
                        pr_cbr,
                        pr_dor8_ok,
                        pr_tox_ok,
                        pr_tox_bad) {
  
  # ----------------------------------------------------------
  # Efficacy gate
  # ----------------------------------------------------------
  # Slide-aligned:
  # Pr(CBR24 >= 40%) > 80%
  # OR
  # Pr(ORR >= theta) > 80%
  
  eff_go <- (pr_cbr > p_eff_go) | (pr_orr > p_eff_go)
  
  # Futility if both efficacy signals are clearly weak.
  eff_low <- (pr_cbr < p_eff_nogo) & (pr_orr < p_eff_nogo)
  
  
  # ----------------------------------------------------------
  # DOR gate
  # ----------------------------------------------------------
  # Pr(p_DOR8 >= 70%) > 70%
  
  dor_go <- pr_dor8_ok > p_dor_go
  
  
  # ----------------------------------------------------------
  # Toxicity gate
  # ----------------------------------------------------------
  # Pr(p_tox < 35%) > 80%
  
  tox_go <- pr_tox_ok > p_tox_go
  
  
  # ----------------------------------------------------------
  # Composite decision
  # ----------------------------------------------------------
  
  case_when(
    
    # Safety first
    pr_tox_bad > p_tox_bad ~ "NO-GO",
    
    # Efficacy futility
    eff_low ~ "NO-GO",
    
    # DOR futility
    pr_dor8_ok < p_dor_nogo ~ "NO-GO",
    
    # Full GO requires efficacy + durability + safety
    eff_go & dor_go & tox_go ~ "GO",
    
    # Otherwise grey zone
    TRUE ~ "MODIFY"
  )
}


# ============================================================
# RUN SIMULATION FOR ONE SAMPLE SIZE
# ============================================================

run_simulation_for_N <- function(N,
                                 N_meas = N,
                                 nsim = 10000) {
  
  pmap_dfr(
    grid,
    function(true_orr, true_cbr, true_tox, true_dor8) {
      
      # ------------------------------------------------------
      # Simulate ORR and CBR coherently
      # ------------------------------------------------------
      
      eff_dat <- simulate_orr_cbr(
        nsim = nsim,
        N_meas = N_meas,
        true_orr = true_orr,
        true_cbr = true_cbr
      )
      
      r_orr <- eff_dat$r_orr
      r_cbr <- eff_dat$r_cbr
      
      
      # ------------------------------------------------------
      # Simulate DOR8 among responders only
      # ------------------------------------------------------
      
      r_dor8 <- simulate_dor8(
        r_orr = r_orr,
        true_dor8 = true_dor8
      )
      
      
      # ------------------------------------------------------
      # Simulate chronic Grade 3+ toxicity
      # ------------------------------------------------------
      
      r_tox <- rbinom(
        n = nsim,
        size = N,
        prob = true_tox
      )
      
      
      # ------------------------------------------------------
      # Posterior probabilities
      # ------------------------------------------------------
      
      pr_orr <- post_prob_ge_beta(
        thr = orr_thr,
        r   = r_orr,
        n   = N_meas,
        a0  = a_orr,
        b0  = b_orr
      )
      
      pr_cbr <- post_prob_ge_beta(
        thr = cbr_thr,
        r   = r_cbr,
        n   = N_meas,
        a0  = a_cbr,
        b0  = b_cbr
      )
      
      # Important:
      # DOR8 denominator is number of responders, not total N.
      #
      # If there are zero responders, the posterior remains equal
      # to the prior Beta(a_dor, b_dor).
      #
      # That is acceptable mathematically, but clinically it means
      # DOR cannot rescue a scenario with no responses.
      
      pr_dor8_ok <- post_prob_ge_beta(
        thr = dor8_thr,
        r   = r_dor8,
        n   = r_orr,
        a0  = a_dor,
        b0  = b_dor
      )
      
      pr_tox_ok <- post_prob_lt_beta(
        thr = tox_thr,
        r   = r_tox,
        n   = N,
        a0  = a_tox,
        b0  = b_tox
      )
      
      pr_tox_bad <- post_prob_ge_beta(
        thr = tox_thr,
        r   = r_tox,
        n   = N,
        a0  = a_tox,
        b0  = b_tox
      )
      
      
      # ------------------------------------------------------
      # Decision
      # ------------------------------------------------------
      
      decision <- decide_A2B2(
        pr_orr     = pr_orr,
        pr_cbr     = pr_cbr,
        pr_dor8_ok = pr_dor8_ok,
        pr_tox_ok  = pr_tox_ok,
        pr_tox_bad = pr_tox_bad
      )
      
      
      # ------------------------------------------------------
      # Store decisions
      # ------------------------------------------------------
      
      tibble(
        true_orr  = true_orr,
        true_cbr  = true_cbr,
        true_tox  = true_tox,
        true_dor8 = true_dor8,
        decision  = decision
      ) %>%
        count(
          true_orr,
          true_cbr,
          true_tox,
          true_dor8,
          decision,
          name = "n"
        )
    }
  ) %>%
    group_by(true_orr, true_cbr, true_tox, true_dor8) %>%
    summarise(
      GO     = sum(n[decision == "GO"])     / sum(n),
      MODIFY = sum(n[decision == "MODIFY"]) / sum(n),
      NO_GO  = sum(n[decision == "NO-GO"])  / sum(n),
      .groups = "drop"
    ) %>%
    mutate(N = N)
}


# ============================================================
# MAIN RUN
# ============================================================

summary_tbl <- run_simulation_for_N(
  N = N,
  N_meas = N_meas,
  nsim = nsim
)

summary_tbl


# Different backfilling sample sizes to evaluate
N_values <- c(10, 15, 20, 25, 30, 35)

# Run simulations for all N values
results_allN <- map_dfr(
  N_values,
  ~ run_simulation_for_N(
    N = .x,
    N_meas = .x,
    nsim = nsim
  )
)

# Convert to long format
results_long <- results_allN %>%
  pivot_longer(
    cols = c(GO, MODIFY, NO_GO),
    names_to = "decision",
    values_to = "prob"
  ) %>%
  mutate(
    prob = 100 * prob
  )


# ============================================================
# DEFINE 4 SCENARIO GROUPS
# ============================================================


results_long_grouped <- results_long %>%
  mutate(
    
    scenario_group = case_when(
      
      true_orr  >= 0.40 &
        true_cbr  >= 0.50 &
        true_tox  <= 0.10 &
        true_dor8 >= 0.90
      ~ "Very favourable",
      
      true_orr  >= 0.30 &
        true_cbr  >= 0.40 &
        true_tox  <= 0.20 &
        true_dor8 >= 0.80
      ~ "Favourable",
      
      true_orr  <= 0.20 &
        true_cbr  <= 0.30 &
        true_tox  >= 0.40 &
        true_dor8 <= 0.60
      ~ "Very unfavourable",
      
      true_orr  <= 0.25 |
        true_cbr  <= 0.30 |
        true_tox  >= 0.30 |
        true_dor8 <= 0.60
      ~ "Unfavourable",
      
      TRUE ~ "Intermediate"
    )
  ) %>%
  filter(
    scenario_group != "Intermediate"
  ) %>%
  mutate(
    scenario_group = factor(
      scenario_group,
      levels = c(
        "Very unfavourable",
        "Unfavourable",
        "Favourable",
        "Very favourable"
      )
    )
  )


results_long_grouped
# ============================================================
# SUMMARY BY N AND SCENARIO GROUP
# ============================================================

summary_by_group_N <- results_long_grouped %>%
  group_by(
    N,
    scenario_group,
    decision
  ) %>%
  summarise(
    mean_prob   = mean(prob),
    median_prob = median(prob),
    min_prob    = min(prob),
    max_prob    = max(prob),
    .groups = "drop"
  )
summary_by_group_N

facet_labels <- c(
  
  "Very unfavourable" =
    "Very unfavourable\nORR ≤20%, CBR ≤30%\nTox ≥40%, DOR8 ≤60%",
  
  "Unfavourable" =
    "Unfavourable\nORR ≤25% or CBR ≤30%\nor Tox ≥30% or DOR8 ≤60%",
  
  "Favourable" =
    "Favourable\nORR ≥30%, CBR ≥40%\nTox ≤20%, DOR8 ≥80%",
  
  "Very favourable" =
    "Very favourable\nORR ≥40%, CBR ≥50%\nTox ≤10%, DOR8 ≥90%"
)

summary_by_group_N <- summary_by_group_N %>%
  mutate(
    scenario_group = factor(
      scenario_group,
      levels = c(
        "Very unfavourable",
        "Unfavourable",
        "Favourable",
        "Very favourable"
      )
    )
  )

decision_palette <- c(
  "GO"     = "#009E73",
  "MODIFY" = "#E69F00",
  "NO_GO"  = "#D22B2B"
)

p_decision_by_N <- ggplot(
  summary_by_group_N,
  aes(
    x = N,
    y = mean_prob,      # yo usaría mean_prob
    colour = decision,
    group = decision
  )
) +
  
  geom_line(linewidth = 1.3) +
  
  geom_point(size = 3) +
  
  facet_wrap(
    ~ scenario_group,
    ncol = 2,
    labeller = labeller(
      scenario_group = facet_labels
    )
  ) +
  
  scale_colour_manual(
    values = decision_palette,
    labels = c(
      "GO",
      "MODIFY",
      "NO-GO"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0,100)
  ) +
  
  labs(
    title = "Decision probabilities by sample size",
    subtitle = "Average probability across simulated scenarios",
    x = "Backfilling sample size (N)",
    y = "Decision probability (%)",
    colour = "Decision"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    
    panel.grid.minor = element_blank()
  )

p_decision_by_N
# ============================================================
# HEATMAPS FOR EACH N
# ============================================================

heatmaps_N <- map(
  N_values,
  function(N_current){
    
    summary_tbl <- run_simulation_for_N(
      N = N_current,
      nsim = nsim
    )
    
    # --------------------------------------------------------
    # Dominant decision
    # --------------------------------------------------------
    
    summary_tbl_dec <- summary_tbl %>%
      rowwise() %>%
      mutate(
        
        decision = {
          
          vals <- c(
            GO      = GO,
            MODIFY  = MODIFY,
            `NO-GO` = NO_GO
          )
          
          mx <- max(vals)
          
          winners <- names(vals)[vals == mx]
          
          if(length(winners) == 1){
            
            winners
            
          } else {
            
            c("NO-GO","MODIFY","GO")[
              c("NO-GO","MODIFY","GO") %in% winners
            ][1]
          }
        },
        
        decision_strength =
          100 * max(GO, MODIFY, NO_GO)
      ) %>%
      ungroup()
    
    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------
    
    summary_tbl_dec <- summary_tbl_dec %>%
      mutate(
        
        cbr_label =
          paste0(
            "True CBR24 = ",
            percent(true_cbr, accuracy = 1)
          ),
        
        dor_label =
          paste0(
            "True DOR8 = ",
            percent(true_dor8, accuracy = 1)
          ),
        
        true_orr =
          factor(
            true_orr,
            levels = sort(unique(true_orr)),
            labels = percent(
              sort(unique(true_orr)),
              accuracy = 1
            )
          ),
        
        true_tox =
          factor(
            true_tox,
            levels = sort(unique(true_tox)),
            labels = percent(
              sort(unique(true_tox)),
              accuracy = 1
            )
          )
      )
    
    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    
    ggplot(
      summary_tbl_dec,
      aes(
        x = true_tox,
        y = true_orr,
        fill = decision
      )
    ) +
      
      geom_tile(
        color = "white",
        linewidth = 0.5
      ) +
      
      geom_text(
        aes(
          label =
            round(
              decision_strength,
              0
            )
        ),
        size = 2.8
      ) +
      
      facet_grid(
        rows = vars(cbr_label),
        cols = vars(dor_label)
      ) +
      
      scale_fill_manual(
        name = "Decision",
        values = c(
          "GO" = "#009E73",
          "MODIFY" = "#F0E442",
          "NO-GO" = "#D22B2B"
        )
      ) +
      
      labs(
        title = paste0("N = ", N_current),
        x = "True chronic G3+ toxicity",
        y = "True ORR"
      ) +
      
      theme_minimal(base_size = 6) +
      
      theme(
        panel.grid = element_blank(),
        
        strip.text =
          element_text(
            face = "bold",
            size = 6
          ),
        
        plot.title =
          element_text(
            face = "bold",
            size = 6,
            hjust = 0.5
          ),
        
        legend.position = "none"
      )
  }
)

# ============================================================
# COMBINE ALL N PANELS
# ============================================================

combined_heatmaps <-
  wrap_plots(
    heatmaps_N,
    ncol = 2,
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )

combined_heatmaps
