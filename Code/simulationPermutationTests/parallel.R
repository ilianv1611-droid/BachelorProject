################################################################################
# Title: Analysis of delayed-start design
# Author: Achille Demares
# Date: 24-02-2026
################################################################################

rm(list = ls())

library(MASS)
library(GGally)
library(nlme)
library(mmrm)
library(tidyverse)
library(SuperLearner)
library(gridExtra)
library(parallel)
library(pbapply)

simulateDSD <- function(n.arm = 10,
                        arms  = 4,
                        tau   = 7,
                        beta0 = 2,
                        beta_time = 0.3,
                        sd_b0 = 0,
                        rho   = 0.6,
                        sd_eps = 0.2) {
  
  times <- 0:tau
  Tn    <- length(times)
  
  df <- data.frame(
    id        = 1:(n.arm * arms),
    crossover = rep(1:arms, each = n.arm)
  )
  
  # Store y, trt_t, trt_t_plot  -> 3*Tn columns
  y_mat <- matrix(NA_real_, nrow = nrow(df), ncol = 3 * Tn)
  
  for (i in seq_len(nrow(df))) {
    
    crossover_i <- df$crossover[i]
    
    # Random intercept
    b0 <- rnorm(1, 0, sd_b0)
    
    # AR(1) residuals
    eps <- numeric(Tn)
    eps[1] <- rnorm(1, 0, sd_eps / sqrt(1 - rho^2))
    
    if (Tn > 1) {
      for (k in 2:Tn) {
        eps[k] <- rho * eps[k - 1] + rnorm(1, 0, sd_eps)
      }
    }
    
    # Treatment starts once time >= crossover_i
    beta_vec  <- as.integer(times >= crossover_i)
    
    times_vec <- c(rep(0, times = crossover_i + 1),
                   1:(tau - crossover_i))
    
    y <- beta0 + (beta_vec * times_vec * beta_time) + b0 + eps
    
    y_mat[i, ] <- c(
      y,
      c(0, beta_vec[1:(length(beta_vec) - 1)]),  # trt_t  (lagged)
      c(beta_vec)                                # trt_t_plot (current)
    )
  }
  
  colnames(y_mat) <- c(
    paste0("y_t", times),
    paste0("trt_t", times),
    paste0("trt_t_plot", times)
  )
  
  # Wide data (contains y_t*, trt_t*, trt_t_plot*)
  df_wide <- cbind(df, as.data.frame(y_mat)) %>%
    mutate(baseline = .data[["y_t0"]])
  
  # --- Long: y
  y_long <- df_wide %>%
    pivot_longer(
      cols = starts_with("y_t"),
      names_to = "time",
      values_to = "y"
    ) %>%
    mutate(time = as.integer(sub("^y_t", "", time)))
  
  # --- Long: trt (ONLY trt_t0, trt_t1, ... not trt_t_plot*)
  trt_long <- df_wide %>%
    pivot_longer(
      cols = matches("^trt_t[0-9]+$"),
      names_to = "time",
      values_to = "trt"
    ) %>%
    mutate(time = as.integer(sub("^trt_t", "", time)))
  
  # --- Long: trt_t_plot
  trt_plot_long <- df_wide %>%
    pivot_longer(
      cols = starts_with("trt_t_plot"),
      names_to = "time",
      values_to = "trt_t_plot"
    ) %>%
    mutate(time = as.integer(sub("^trt_t_plot", "", time)))
  
  # Join y + trt + trt_t_plot
  df_long <- y_long %>%
    left_join(trt_long,      by = c("id", "crossover", "time")) %>%
    left_join(trt_plot_long, by = c("id", "crossover", "time")) %>%
    arrange(id, time) %>%
    select(id, crossover, y, trt, trt_t_plot, time)
  
  # Add baseline (time 0) to long
  baseline_df <- df_long %>%
    filter(time == 0) %>%
    select(id, baseline = y)
  
  df_long <- df_long %>%
    left_join(baseline_df, by = "id")
  
  # baseline_pre = y at crossover time
  baseline_pre_df <- df_long %>%
    filter(time == crossover) %>%
    select(id, baseline_pre = y)
  
  # Add to long
  df_long <- df_long %>%
    left_join(baseline_pre_df, by = "id")
  
  # Add to wide
  df_wide <- df_wide %>%
    left_join(baseline_pre_df, by = "id")
  
  return(list(df_wide = df_wide,
              df_long = df_long))
}

################################################################################
# VISUALISATIONS
################################################################################

sim <- simulateDSD()

df_wide <- sim$df_wide
df_long <- sim$df_long

df_arm_mean <- df_long %>%
  group_by(crossover, time) %>%
  summarise(
    y_mean = mean(y),
    trt_t_plot = max(trt_t_plot),   # or first(trt_t_plot)
    .groups = "drop"
  )

ggplot(df_arm_mean, aes(x = time, y = y_mean, group = crossover, color = factor(trt_t_plot))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("0" = "black", "1" = "red")) +
  labs(x = "Time", y = "Mean Y", color = "Treatment") +
  theme_bw()

df_long %>%
  ggplot(aes(x = time, y = y, group = id, color = factor(trt_t_plot))) +
  geom_line(alpha = 0.6, linewidth = 0.7) +
  facet_wrap(~ crossover, labeller = label_both) +
  scale_color_manual(
    values = c("0" = "black", "1" = "red")   # explicitly map trt=0 → black, trt=1 → red
  ) +
  labs(
    x = "Time",
    y = "y"
  ) + theme_bw() + 
  theme(
    legend.position = "none"
  )

df_long %>%
  ggplot(aes(x = time, y = y - baseline, group = id, color = factor(trt_t_plot))) +
  geom_line(alpha = 0.6, linewidth = 0.7) +
  facet_wrap(~ crossover, labeller = label_both) +
  scale_color_manual(
    values = c("0" = "black", "1" = "red")   # explicitly map trt=0 → black, trt=1 → red
  ) +
  labs(
    x = "Time",
    y = "cfb (start follow-up)"
  ) + theme_bw() + 
  theme(
    legend.position = "none"
  )

df_long %>%
  ggplot(aes(x = time, y = y - baseline_pre, group = id, color = factor(trt_t_plot))) +
  geom_line(alpha = 0.6, linewidth = 0.7) +
  facet_wrap(~ crossover, labeller = label_both) +
  scale_color_manual(
    values = c("0" = "black", "1" = "red")   # explicitly map trt=0 → black, trt=1 → red
  ) +
  labs(
    x = "Time",
    y = "cfb (start treatment)"
  ) +  theme_bw() + 
  theme(
    legend.position = "none"
  )


################################################################################
# SIMULATIONS
################################################################################

# ----- Estimator -----
estimators <- function(data){
  
  fit_lme <- lme(
    y ~  time + trt + time:trt,
    random = ~ 1 | id,
    correlation = corAR1(form = ~ time | id),
    data = data$df_long,
    method = "REML"
  )
  
  est <- summary(fit_lme)$tTable["time:trt", c(1, 2, 5)]
  est
}

# Test
set.seed(1)
sim <- simulateDSD()
estimators(sim)
sim$df_long

# ----- Wrapper for parallel -----
powerDSD.seed <- function(param){
  set.seed(param$seed)
  data <- simulateDSD()
  estimators(data)
}

# ----- Parallel setup -----
n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(nlme)
  library(tidyverse)
})

clusterExport(
  cl,
  varlist = c("simulateDSD", "powerDSD.seed", "estimators")
)

# ----- Seeds -----
n.seed <- 500

params <- expand_grid(
  seed = 1:n.seed
)

param_list <- split(params, seq_len(nrow(params)))

# ----- Run simulation -----
result <- pblapply(param_list, powerDSD.seed, cl = cl)

stopCluster(cl)

# ----- Bind results -----
results.sim <- data.frame(do.call(rbind, result))
results.sim$seed <- 1:n.seed
rownames(results.sim) <- NULL
colnames(results.sim) <- c("slope", "se", "p.val", "seed")

head(results.sim)
