# Project 01 / Analysis 03
# Cauchy Location Model with Grid Approximation
#
# Purpose:
#   Approximate the posterior distribution for the unknown center of a Cauchy
#   distribution using a grid, then simulate posterior and posterior predictive draws.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/03-cauchy-location-grid-approximation"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(927001)

# -----------------------------
# Data and posterior grid
# -----------------------------

observed_y <- c(53, 54, 55, 56.5, 57.5)
n_grid <- 1000
theta_grid <- seq(0, 100, length.out = n_grid)
dtheta <- theta_grid[2] - theta_grid[1]

# Log unnormalized posterior:
# p(theta | y) proportional to prod_i 1 / (1 + (y_i - theta)^2)
# for theta in [0, 100].
log_weight <- sapply(theta_grid, function(theta) {
  -sum(log(1 + (observed_y - theta)^2))
})

# Stabilize before exponentiating.
log_weight <- log_weight - max(log_weight)
weight <- exp(log_weight)

posterior_density <- weight / sum(weight * dtheta)
posterior_mass <- weight / sum(weight)

posterior_grid <- data.frame(
  theta = theta_grid,
  posterior_density = posterior_density,
  posterior_mass = posterior_mass
)

write.csv(
  posterior_grid,
  file = file.path(results_dir, "cauchy_posterior_grid.csv"),
  row.names = FALSE
)

# -----------------------------
# Posterior and posterior predictive simulation
# -----------------------------

theta_draws <- sample(
  theta_grid,
  size = 1000,
  replace = TRUE,
  prob = posterior_mass
)

y6_draws <- rcauchy(
  n = length(theta_draws),
  location = theta_draws,
  scale = 1
)

posterior_summary <- data.frame(
  quantity = c("mean", "median", "sd", "q025", "q975", "grid_mode"),
  value = c(
    mean(theta_draws),
    median(theta_draws),
    sd(theta_draws),
    quantile(theta_draws, 0.025),
    quantile(theta_draws, 0.975),
    theta_grid[which.max(posterior_density)]
  )
)

predictive_summary <- data.frame(
  quantity = c("mean", "median", "sd", "q025", "q975"),
  value = c(
    mean(y6_draws),
    median(y6_draws),
    sd(y6_draws),
    quantile(y6_draws, 0.025),
    quantile(y6_draws, 0.975)
  )
)

write.csv(
  posterior_summary,
  file = file.path(results_dir, "theta_posterior_summary.csv"),
  row.names = FALSE
)

write.csv(
  predictive_summary,
  file = file.path(results_dir, "y6_predictive_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(theta = theta_draws),
  file = file.path(results_dir, "theta_posterior_draws.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(y6 = y6_draws),
  file = file.path(results_dir, "y6_predictive_draws.csv"),
  row.names = FALSE
)

# -----------------------------
# Plots
# -----------------------------

png(file.path(fig_dir, "01_cauchy_posterior_density.png"), width = 1200, height = 800, res = 150)
plot(
  theta_grid,
  posterior_density,
  type = "l",
  xlab = expression(theta),
  ylab = "Posterior density",
  main = "Cauchy location posterior: grid approximation"
)
rug(observed_y)
dev.off()

png(file.path(fig_dir, "02_theta_posterior_draws.png"), width = 1200, height = 800, res = 150)
hist(
  theta_draws,
  breaks = 30,
  xlab = expression(theta),
  main = expression("Posterior draws of " * theta)
)
dev.off()

png(file.path(fig_dir, "03_y6_posterior_predictive_draws.png"), width = 1200, height = 800, res = 150)
hist(
  y6_draws,
  breaks = 40,
  xlab = expression(y[6]),
  main = expression("Posterior predictive draws of " * y[6])
)
dev.off()

# A zoomed version is useful because Cauchy predictive draws can contain very large outliers.
y6_limits <- quantile(y6_draws, probs = c(0.01, 0.99))

png(file.path(fig_dir, "04_y6_posterior_predictive_draws_zoom.png"), width = 1200, height = 800, res = 150)
hist(
  y6_draws[y6_draws >= y6_limits[1] & y6_draws <= y6_limits[2]],
  breaks = 40,
  xlab = expression(y[6]),
  main = expression("Posterior predictive draws of " * y[6] * " (middle 98%)")
)
dev.off()

cat("Analysis 03 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")