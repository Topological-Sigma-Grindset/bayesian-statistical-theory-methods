# Project 01 / Analysis 01
# Beta-Binomial Foundations
#
# Purpose:
#   Demonstrate basic beta-binomial calculations used in Bayesian updating:
#   prior predictive probabilities, conjugate posterior updating, posterior
#   means as weighted averages, and variance behavior before/after observing data.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/01-beta-binomial-foundations"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

beta_mean <- function(a, b) a / (a + b)

beta_var <- function(a, b) {
  (a * b) / ((a + b)^2 * (a + b + 1))
}

posterior_beta_params <- function(a, b, y, n) {
  c(a_post = a + y, b_post = b + n - y)
}

# -----------------------------
# 1. Prior predictive under Uniform(0,1)
# -----------------------------

# If theta ~ Uniform(0,1) = Beta(1,1) and Y | theta ~ Binomial(n, theta),
# then P(Y = k) = 1 / (n + 1), for k = 0, ..., n.
n_demo <- 10
k_values <- 0:n_demo
prior_predictive <- choose(n_demo, k_values) * beta(k_values + 1, n_demo - k_values + 1)

prior_predictive_results <- data.frame(
  n = n_demo,
  k = k_values,
  prior_predictive_probability = prior_predictive,
  expected_uniform_probability = 1 / (n_demo + 1)
)

write.csv(
  prior_predictive_results,
  file = file.path(results_dir, "prior_predictive_uniform_binomial.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "01_prior_predictive_uniform_binomial.png"), width = 1200, height = 800, res = 150)
barplot(
  height = prior_predictive,
  names.arg = k_values,
  xlab = "Number of successes k",
  ylab = "Prior predictive probability",
  main = "Prior predictive distribution under Uniform(0,1) prior"
)
dev.off()

# -----------------------------
# 2. Posterior mean as weighted average
# -----------------------------

# Example values used only to numerically illustrate the algebra.
a <- 2
b <- 5
n <- 20
y <- 8

post <- posterior_beta_params(a, b, y, n)
a_post <- post["a_post"]
b_post <- post["b_post"]

prior_mean <- beta_mean(a, b)
sample_proportion <- y / n
posterior_mean <- beta_mean(a_post, b_post)

weight_prior <- (a + b) / (a + b + n)
weight_data <- n / (a + b + n)
posterior_mean_from_weights <- weight_prior * prior_mean + weight_data * sample_proportion

posterior_mean_results <- data.frame(
  a = a,
  b = b,
  n = n,
  y = y,
  prior_mean = prior_mean,
  sample_proportion = sample_proportion,
  posterior_mean = posterior_mean,
  prior_weight = weight_prior,
  data_weight = weight_data,
  posterior_mean_from_weighted_average = posterior_mean_from_weights
)

write.csv(
  posterior_mean_results,
  file = file.path(results_dir, "posterior_mean_weighted_average.csv"),
  row.names = FALSE
)

# -----------------------------
# 3. Uniform prior: posterior variance decreases for n >= 1
# -----------------------------

variance_rows <- list()
row_id <- 1

for (n_i in 1:25) {
  for (y_i in 0:n_i) {
    a_prior <- 1
    b_prior <- 1
    params <- posterior_beta_params(a_prior, b_prior, y_i, n_i)
    variance_rows[[row_id]] <- data.frame(
      n = n_i,
      y = y_i,
      prior_variance = beta_var(a_prior, b_prior),
      posterior_variance = beta_var(params["a_post"], params["b_post"])
    )
    row_id <- row_id + 1
  }
}

uniform_variance_results <- do.call(rbind, variance_rows)
uniform_variance_results$posterior_less_than_prior <- (
  uniform_variance_results$posterior_variance < uniform_variance_results$prior_variance
)

write.csv(
  uniform_variance_results,
  file = file.path(results_dir, "uniform_prior_variance_check.csv"),
  row.names = FALSE
)

# Show, for each n, the maximum posterior variance over y.
max_var_by_n <- aggregate(
  posterior_variance ~ n,
  data = uniform_variance_results,
  FUN = max
)
max_var_by_n$prior_variance <- beta_var(1, 1)

png(file.path(fig_dir, "02_uniform_prior_max_posterior_variance.png"), width = 1200, height = 800, res = 150)
plot(
  max_var_by_n$n,
  max_var_by_n$posterior_variance,
  type = "b",
  xlab = "Number of trials n",
  ylab = "Variance",
  main = "Maximum posterior variance under Uniform(0,1) prior"
)
abline(h = beta_var(1, 1), lty = 2)
legend(
  "topright",
  legend = c("max posterior variance", "prior variance"),
  lty = c(1, 2),
  pch = c(1, NA),
  bty = "n"
)
dev.off()

# -----------------------------
# 4. Example where posterior variance increases
# -----------------------------

# With theta ~ Beta(1,3), observing y = 1 success out of n = 1 trial gives
# theta | y ~ Beta(2,3). The posterior variance is larger than the prior variance.
a_example <- 1
b_example <- 3
n_example <- 1
y_example <- 1
post_example <- posterior_beta_params(a_example, b_example, y_example, n_example)

variance_increase_example <- data.frame(
  a_prior = a_example,
  b_prior = b_example,
  n = n_example,
  y = y_example,
  a_post = post_example["a_post"],
  b_post = post_example["b_post"],
  prior_mean = beta_mean(a_example, b_example),
  posterior_mean = beta_mean(post_example["a_post"], post_example["b_post"]),
  prior_variance = beta_var(a_example, b_example),
  posterior_variance = beta_var(post_example["a_post"], post_example["b_post"])
)

variance_increase_example$posterior_variance_larger <- (
  variance_increase_example$posterior_variance > variance_increase_example$prior_variance
)

write.csv(
  variance_increase_example,
  file = file.path(results_dir, "variance_increase_example.csv"),
  row.names = FALSE
)

x_grid <- seq(0, 1, length.out = 1000)

png(file.path(fig_dir, "03_variance_increase_example.png"), width = 1200, height = 800, res = 150)
plot(
  x_grid,
  dbeta(x_grid, a_example, b_example),
  type = "l",
  xlab = expression(theta),
  ylab = "Density",
  main = "Prior and posterior when posterior variance increases"
)
lines(x_grid, dbeta(x_grid, post_example["a_post"], post_example["b_post"]), lty = 2)
legend(
  "topright",
  legend = c("Prior: Beta(1,3)", "Posterior: Beta(2,3)"),
  lty = c(1, 2),
  bty = "n"
)
dev.off()

cat("Analysis 01 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")