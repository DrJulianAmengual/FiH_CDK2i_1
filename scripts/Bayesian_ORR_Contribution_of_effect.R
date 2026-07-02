## ============================================
## Bayesian ORR vs fixed p0 (Simulates external control arm RWE)
## ============================================

set.seed(123)

# -------------------------------
# 1. Inputs
# -------------------------------
N_vec   <- c(20, 40, 60,80,100)  ## Different samples size (per arm)
p0_vec  <- c(0.2, 0.3, 0.4,0.5)  ## True ORR (historical benchmark)
pt_vec  <- c(0.25, 0.3, 0.4, 0.5, 0.6,0.7) # True ORR (experimental arm)

nsim <- 5000 ## Number of repetitions

# Jeffreys prior
a_prior <- 0.5
b_prior <- 0.5

##P(ORR_exp > p0_vec + delta) > go_treshold

go_threshold <- 0.8


# -------------------------------
# 2. grid with all values for P0 and P_true
# -------------------------------
grid <- expand.grid(N = N_vec,
                    p0 = p0_vec,
                    p_true = pt_vec)

# We analyze cases such as pt > p0
grid <- subset(grid, p_true > p0)

##P(ORR_exp > p0_vec + delta) > go_treshold
delta <- 0.05
# -------------------------------
# 3. Posterior probability
# -------------------------------
simulate_many <- function(N, p_true, p0, nsim, a_prior, b_prior,go_threshold, delta ) {
  
  y <- rbinom(nsim, N, p_true) # Binomial N, P_true for the experimental arm
  
  a_post <- a_prior + y
  b_post <- b_prior + N - y
  
  prob_sup <- 1 - pbeta(p0, a_post, b_post)
  
  mean(prob_sup+delta > go_threshold ) 
}


# -------------------------------
# 4. Simulation
# -------------------------------
grid$prob_go <- NA

for (i in 1:nrow(grid)) {
  
  grid$prob_go[i] <- simulate_many(
    N        = grid$N[i],
    p_true   = grid$p_true[i],
    p0       = grid$p0[i],
    go_threshold = go_threshold,
    nsim     = nsim,
    a_prior  = a_prior,
    b_prior  = b_prior,
    delta = delta
  )
}

# -------------------------------
# 5.Final output
# -------------------------------
grid$delta <- grid$p_true - grid$p0

print(grid)

# -------------------------------
# 6. Heatmap: p0 vs p_true by N
# -------------------------------

library(ggplot2)


grid$p0_f     <- factor(grid$p0)
grid$p_true_f <- factor(grid$p_true)
grid$N_f      <- factor(grid$N)

p_heatmap <- ggplot(
  grid,
  aes(x = p0_f, y = p_true_f, fill = prob_go)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", prob_go)),
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~ N_f, nrow = 1) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    limits = c(0, 1),
    name = "P(GO)"
  ) +
  labs(
    title = "Bayesian ORR design: P(GO) by fixed p0 (historical RWE benchmark) and true experimental ORR",
    subtitle = paste0(
      "Only scenarios with p_true > p0 + 0.05; GO if P(p_t > p0 + 0.05 | data) > ",
      go_threshold
    ),
    x = "Fixed control ORR (historical RWE Benchmark, p0)",
    y = "True experimental ORR (p_true)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_heatmap)



## ============================================================
## Bayesian randomized add-on ORR design
## Control: backbone
## Experimental: backbone + new compound
## Decision: GO if P(p_exp > p_ctrl + delta_min | data) > go_threshold
## ============================================================

set.seed(123)

# -------------------------------
# 1. Inputs
# -------------------------------

# Sample size per arm
N_per_arm_vec <- c(20, 40, 60, 80, 100)

# True ORR in control arm: backbone
p_ctrl_vec <- c(0.20, 0.30, 0.40, 0.50)

# True ORR in experimental arm: backbone + new compound
p_exp_vec <- c(0.25, 0.30, 0.40, 0.50, 0.60, 0.70)

# Number of simulated trials per scenario
nsim <- 5000

# Prior for both arms
a_prior <- 0.5
b_prior <- 0.5

# Bayesian decision rule
go_threshold <- 0.90

# Clinically meaningful add-on margin
# delta_min = 0 means any improvement
# delta_min = 0.05 means at least +5% ORR improvement
delta_min <- 0.05

# Number of posterior samples used to estimate
# P(p_exp > p_ctrl + delta_min | data)
n_post <- 5000

# -------------------------------
# 2. Function: simulate many randomized trials
# -------------------------------

simulate_randomized_addon <- function(
    N_per_arm,
    p_ctrl_true,
    p_exp_true,
    nsim,
    a_prior,
    b_prior,
    go_threshold,
    delta_min,
    n_post
) {
  
  # Simulate observed responders in each arm
  y_ctrl <- rbinom(nsim, N_per_arm, p_ctrl_true)
  y_exp  <- rbinom(nsim, N_per_arm, p_exp_true)
  
  # Posterior parameters
  a_ctrl_post <- a_prior + y_ctrl
  b_ctrl_post <- b_prior + N_per_arm - y_ctrl
  
  a_exp_post <- a_prior + y_exp
  b_exp_post <- b_prior + N_per_arm - y_exp
  
  # Store posterior probability for each simulated trial
  post_prob_addon <- numeric(nsim)
  
  for (j in seq_len(nsim)) {
    
    p_ctrl_samp <- rbeta(n_post, a_ctrl_post[j], b_ctrl_post[j])
    p_exp_samp  <- rbeta(n_post, a_exp_post[j],  b_exp_post[j])
    
    post_prob_addon[j] <- mean(p_exp_samp > p_ctrl_samp + delta_min)
  }
  
  # Probability of GO across simulated trials
  prob_go <- mean(post_prob_addon > go_threshold)
  
  # Optional summaries
  mean_post_prob <- mean(post_prob_addon)
  median_post_prob <- median(post_prob_addon)
  
  return(
    data.frame(
      prob_go = prob_go,
      mean_post_prob = mean_post_prob,
      median_post_prob = median_post_prob
    )
  )
}

# -------------------------------
# 3. Scenario grid
# -------------------------------

grid <- expand.grid(
  N_per_arm = N_per_arm_vec,
  p_ctrl_true = p_ctrl_vec,
  p_exp_true = p_exp_vec
)

# Keep only add-on scenarios
grid <- subset(grid, p_exp_true > p_ctrl_true)

# Add true treatment effect
grid$delta_true <- grid$p_exp_true - grid$p_ctrl_true


# -------------------------------
# 4. Run simulations
# -------------------------------

grid$prob_go <- NA
grid$mean_post_prob <- NA
grid$median_post_prob <- NA

for (i in seq_len(nrow(grid))) {
  
  out_i <- simulate_randomized_addon(
    N_per_arm    = grid$N_per_arm[i],
    p_ctrl_true  = grid$p_ctrl_true[i],
    p_exp_true   = grid$p_exp_true[i],
    nsim         = nsim,
    a_prior      = a_prior,
    b_prior      = b_prior,
    go_threshold = go_threshold,
    delta_min    = delta_min,
    n_post       = n_post
  )
  
  grid$prob_go[i]          <- out_i$prob_go
  grid$mean_post_prob[i]   <- out_i$mean_post_prob
  grid$median_post_prob[i] <- out_i$median_post_prob
}

# Total sample size
grid$N_total <- 2 * grid$N_per_arm

print(grid)

# -------------------------------
# 5. Heatmap by sample size
# -------------------------------

library(ggplot2)

grid$p_ctrl_f <- factor(grid$p_ctrl_true)
grid$p_exp_f  <- factor(grid$p_exp_true)
grid$N_f      <- factor(grid$N_per_arm)

p_heatmap <- ggplot(
  grid,
  aes(x = p_ctrl_f, y = p_exp_f, fill = prob_go)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", prob_go)),
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~ N_f, nrow = 1) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    limits = c(0, 1),
    name = "P(GO)"
  ) +
  labs(
    title = "Bayesian randomized add-on ORR design",
    subtitle = paste0(
      "GO if P(p_exp > p_ctrl + ",
      delta_min,
      " | data) > ",
      go_threshold,
      "; N shown is per arm"
    ),
    x = "True ORR control arm: backbone",
    y = "True ORR experimental arm: backbone + new compound"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_heatmap)

ggsave(
  filename = "Bayesian_randomized_addon_ORR_heatmap_by_N.png",
  plot = p_heatmap,
  width = 14,
  height = 4.5,
  dpi = 300
)


# ============================================================
# Plot: P(GO) by relative ORR increase vs control
# ============================================================

library(ggplot2)
library(dplyr)

# -------------------------------
# 1. Create relative increase variable
# -------------------------------

grid_plot <- grid %>%
  mutate(
    rel_increase = (p_exp_true - p_ctrl_true) / p_ctrl_true,
    rel_increase_pct = 100 * rel_increase,
    abs_increase = p_exp_true - p_ctrl_true,
    N_label = paste0("N = ", N_per_arm, " per arm"),
    p_ctrl_label = paste0("Control ORR = ", round(100 * p_ctrl_true), "%")
  )

# Optional: keep only relative increases close to 10%, 20%, 30%, 40%, 50%
# This is useful if your grid contains exactly those scenarios.
grid_plot <- grid_plot %>%
  filter(round(rel_increase_pct, 0) %in% c(10, 20, 30, 40, 50))


# -------------------------------
# 2. Main plot
# -------------------------------

p_rel <- ggplot(
  grid_plot,
  aes(
    x = rel_increase_pct,
    y = prob_go,
    color = N_label,
    group = N_label
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  facet_wrap(~ p_ctrl_label, nrow = 1) +
  scale_x_continuous(
    breaks = c(10, 20, 30, 40, 50),
    limits = c(10, 50)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "Bayesian randomized add-on ORR design",
    subtitle = paste0(
      "GO if P(p_exp > p_ctrl + ",
      delta_min,
      " | data) > ",
      go_threshold,
      "; experimental ORR expressed as relative increase vs control"
    ),
    x = "Relative increase in experimental ORR vs control (%)",
    y = "Probability of GO",
    color = "Sample size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_rel)


# -------------------------------
# 3. Save plot
# -------------------------------

ggsave(
  filename = "Bayesian_randomized_addon_ORR_PGO_by_relative_increase.png",
  plot = p_rel,
  width = 14,
  height = 5,
  dpi = 300
)
