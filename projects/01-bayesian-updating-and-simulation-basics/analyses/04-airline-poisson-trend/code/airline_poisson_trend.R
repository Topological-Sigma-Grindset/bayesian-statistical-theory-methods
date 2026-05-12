# Project 01 / Analysis 04
# Poisson Trend Model for Airline Fatal Accident Counts
#
# Purpose:
#   Use grid sampling for a two-parameter Poisson trend model and simulate the
#   posterior predictive distribution for the next year.
#
# Model:
#   y_t | alpha0, beta ~ Poisson(lambda_t)
#   lambda_t = alpha0 + beta * (year_t - 1980)
#   with a flat prior over parameter values that keep all lambda_t positive.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/04-airline-poisson-trend"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")
data_dir <- file.path(out_dir, "data")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(927001)

# -----------------------------
# Data
# -----------------------------

load_airline_data <- function() {
  data_path <- file.path(data_dir, "planes.txt")

  if (file.exists(data_path)) {
    dat <- read.table(data_path, header = TRUE, sep = "\t")

    if (!all(c("year", "fatal") %in% names(dat))) {
      stop("planes.txt must contain columns named 'year' and 'fatal'.")
    }

    dat$source <- "data/planes.txt"
    return(dat)
  }

  # Fallback data used to make the analysis reproducible without storing the
  # original course data file. These are fatal accident counts for 1976--1985.
  data.frame(
    year = 1976:1985,
    fatal = c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22),
    source = "embedded_fatal_accident_counts"
  )
}

dat <- load_airline_data()
y <- dat$fatal
year <- dat$year
x <- year - 1980

write.csv(
  dat,
  file = file.path(results_dir, "airline_data_used.csv"),
  row.names = FALSE
)

# -----------------------------
# Posterior grid
# -----------------------------

log_posterior <- function(alpha0, beta) {
  lambda <- alpha0 + beta * x
  if (any(lambda <= 0)) return(-Inf)

  # Dropping Poisson normalizing constants that do not depend on parameters.
  sum(y * log(lambda) - lambda)
}

alpha0_grid <- seq(10, 50, length.out = 250)
beta_grid <- seq(-2, 0.6, length.out = 250)

log_post_grid <- outer(
  alpha0_grid,
  beta_grid,
  Vectorize(function(alpha0, beta) log_posterior(alpha0, beta))
)

log_post_grid <- log_post_grid - max(log_post_grid[is.finite(log_post_grid)])
posterior_weight <- exp(log_post_grid)
posterior_weight[!is.finite(posterior_weight)] <- 0

if (sum(posterior_weight) == 0) {
  stop("All grid points have zero posterior weight. Widen or shift the grid.")
}

posterior_prob <- posterior_weight / sum(posterior_weight)

posterior_grid <- data.frame(
  alpha0 = rep(alpha0_grid, times = length(beta_grid)),
  beta = rep(beta_grid, each = length(alpha0_grid)),
  posterior_mass = as.vector(posterior_prob)
)

write.csv(
  posterior_grid,
  file = file.path(results_dir, "airline_posterior_grid.csv"),
  row.names = FALSE
)

# -----------------------------
# Posterior sampling
# -----------------------------

n_draws <- 1000

idx <- sample.int(
  length(posterior_prob),
  size = n_draws,
  replace = TRUE,
  prob = as.vector(posterior_prob)
)

alpha_index <- ((idx - 1) %% length(alpha0_grid)) + 1
beta_index <- ((idx - 1) %/% length(alpha0_grid)) + 1

alpha0_draws <- alpha0_grid[alpha_index]
beta_draws <- beta_grid[beta_index]

posterior_draws <- data.frame(
  alpha0 = alpha0_draws,
  beta = beta_draws
)

write.csv(
  posterior_draws,
  file = file.path(results_dir, "airline_alpha0_beta_draws.csv"),
  row.names = FALSE
)

# -----------------------------
# Posterior predictive simulation for 1986
# -----------------------------

year_star <- 1986
x_star <- year_star - 1980
lambda_star <- alpha0_draws + beta_draws * x_star

if (any(lambda_star <= 0)) {
  stop("Some posterior predictive Poisson rates are non-positive. Check the grid.")
}

y_star <- rpois(n = n_draws, lambda = lambda_star)

predictive_draws <- data.frame(
  year = year_star,
  lambda = lambda_star,
  y_star = y_star
)

write.csv(
  predictive_draws,
  file = file.path(results_dir, "airline_predictive_draws_1986.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  quantity = c(
    "alpha0_mean", "alpha0_median", "alpha0_q025", "alpha0_q975",
    "beta_mean", "beta_median", "beta_q025", "beta_q975",
    "lambda_1986_mean", "lambda_1986_median", "lambda_1986_q025", "lambda_1986_q975",
    "y_1986_mean", "y_1986_median", "y_1986_q025", "y_1986_q975"
  ),
  value = c(
    mean(alpha0_draws), median(alpha0_draws), quantile(alpha0_draws, 0.025), quantile(alpha0_draws, 0.975),
    mean(beta_draws), median(beta_draws), quantile(beta_draws, 0.025), quantile(beta_draws, 0.975),
    mean(lambda_star), median(lambda_star), quantile(lambda_star, 0.025), quantile(lambda_star, 0.975),
    mean(y_star), median(y_star), quantile(y_star, 0.025), quantile(y_star, 0.975)
  )
)

write.csv(
  summary_table,
  file = file.path(results_dir, "airline_poisson_trend_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Plots
# -----------------------------

png(file.path(fig_dir, "01_joint_posterior_contour.png"), width = 1200, height = 900, res = 150)
contour(
  x = beta_grid,
  y = alpha0_grid,
  z = t(posterior_prob),
  nlevels = 10,
  xlab = expression(beta),
  ylab = expression(alpha[0]),
  main = expression("Joint posterior grid for " * alpha[0] * " and " * beta)
)
dev.off()

png(file.path(fig_dir, "02_posterior_draws_alpha0_beta.png"), width = 1200, height = 900, res = 150)
plot(
  beta_draws,
  alpha0_draws,
  pch = 16,
  cex = 0.5,
  xlab = expression(beta),
  ylab = expression(alpha[0]),
  main = expression("Posterior draws of " * alpha[0] * " and " * beta)
)
dev.off()

png(file.path(fig_dir, "03_predictive_distribution_1986.png"), width = 1200, height = 800, res = 150)
hist(
  y_star,
  breaks = 20,
  xlab = expression(y^"*"),
  main = "Posterior predictive distribution for fatal accidents in 1986"
)
dev.off()

png(file.path(fig_dir, "04_fitted_trend_draws.png"), width = 1200, height = 800, res = 150)
plot(
  year,
  y,
  pch = 16,
  xlab = "Year",
  ylab = "Fatal accidents",
  main = "Observed counts and posterior trend draws"
)
for (s in sample(seq_len(n_draws), size = min(100, n_draws))) {
  lines(year, alpha0_draws[s] + beta_draws[s] * x, lwd = 0.5)
}
points(year, y, pch = 16)
dev.off()

cat("Analysis 04 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")  