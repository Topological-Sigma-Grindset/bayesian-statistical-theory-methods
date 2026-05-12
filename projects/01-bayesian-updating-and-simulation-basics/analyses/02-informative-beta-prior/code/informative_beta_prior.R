# Project 01 / Analysis 02
# Informative Beta Prior for a Survey Proportion
#
# Purpose:
#   Construct an informative Beta(a,b) prior from a prior mean and standard
#   deviation, update it with binomial survey data, and summarize the posterior.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/02-informative-beta-prior"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Prior specification
# -----------------------------

prior_mean_target <- 0.4
prior_sd_target <- 0.1
prior_var_target <- prior_sd_target^2

# For theta ~ Beta(a,b), let total_prior_size = a + b.
# Then Var(theta) = m(1-m) / (total_prior_size + 1).
total_prior_size <- prior_mean_target * (1 - prior_mean_target) / prior_var_target - 1
prior_a <- prior_mean_target * total_prior_size
prior_b <- (1 - prior_mean_target) * total_prior_size

# -----------------------------
# Binomial data and posterior update
# -----------------------------

n_survey <- 100
y_delivery <- 25

posterior_a <- prior_a + y_delivery
posterior_b <- prior_b + n_survey - y_delivery

posterior_mean <- posterior_a / (posterior_a + posterior_b)
posterior_var <- (posterior_a * posterior_b) / (
  (posterior_a + posterior_b)^2 * (posterior_a + posterior_b + 1)
)

summary_table <- data.frame(
  quantity = c(
    "target_prior_mean",
    "target_prior_sd",
    "prior_a",
    "prior_b",
    "survey_n",
    "survey_y",
    "posterior_a",
    "posterior_b",
    "posterior_mean",
    "posterior_variance",
    "posterior_sd"
  ),
  value = c(
    prior_mean_target,
    prior_sd_target,
    prior_a,
    prior_b,
    n_survey,
    y_delivery,
    posterior_a,
    posterior_b,
    posterior_mean,
    posterior_var,
    sqrt(posterior_var)
  )
)

write.csv(
  summary_table,
  file = file.path(results_dir, "beta_prior_posterior_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Plots
# -----------------------------

theta_grid <- seq(0, 1, length.out = 1000)
prior_density <- dbeta(theta_grid, prior_a, prior_b)
posterior_density <- dbeta(theta_grid, posterior_a, posterior_b)

png(file.path(fig_dir, "01_beta_prior_density.png"), width = 1200, height = 800, res = 150)
plot(
  theta_grid,
  prior_density,
  type = "l",
  xlab = expression(theta),
  ylab = "Density",
  main = "Informative prior: Beta(9.2, 13.8)"
)
dev.off()

png(file.path(fig_dir, "02_beta_posterior_density.png"), width = 1200, height = 800, res = 150)
plot(
  theta_grid,
  posterior_density,
  type = "l",
  xlab = expression(theta),
  ylab = "Density",
  main = "Posterior after 25 successes out of 100"
)
dev.off()

png(file.path(fig_dir, "03_prior_posterior_overlay.png"), width = 1200, height = 800, res = 150)
plot(
  theta_grid,
  prior_density,
  type = "l",
  ylim = range(c(prior_density, posterior_density)),
  xlab = expression(theta),
  ylab = "Density",
  main = "Prior and posterior densities"
)
lines(theta_grid, posterior_density, lty = 2)
abline(v = y_delivery / n_survey, lty = 3)
legend(
  "topright",
  legend = c("Prior", "Posterior", "Observed proportion"),
  lty = c(1, 2, 3),
  bty = "n"
)
dev.off()

cat("Analysis 02 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")