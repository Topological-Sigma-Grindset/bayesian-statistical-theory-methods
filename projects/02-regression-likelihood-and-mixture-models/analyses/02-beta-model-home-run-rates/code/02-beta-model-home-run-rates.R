# Project 02 - Analysis 02
# Beta model for Lahman home-run rates
#
# Required package: Lahman

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/02-beta-model-home-run-rates"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("Lahman", quietly = TRUE)) {
  stop("Analysis 02 requires the Lahman package. Install it with install.packages(\"Lahman\").")
}

data("Batting", package = "Lahman")

valid <- rep(TRUE, nrow(Batting))
valid[is.na(Batting$yearID) | Batting$yearID < 1970] <- FALSE
valid[is.na(Batting$HR) | Batting$HR < 1] <- FALSE
valid[is.na(Batting$AB) | Batting$AB < 200] <- FALSE

Batting.restricted <- Batting[valid, ]
Batting.restricted$HRrate <- Batting.restricted$HR / Batting.restricted$AB
Batting.restricted <- Batting.restricted[
  is.finite(Batting.restricted$HRrate) & Batting.restricted$HRrate > 0 & Batting.restricted$HRrate < 1,
]

HRrate <- Batting.restricted$HRrate
n <- length(HRrate)
S1 <- sum(log(HRrate))
S2 <- sum(log1p(-HRrate))

loglik_beta <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) return(-Inf)
  -n * lbeta(alpha, beta) + (alpha - 1) * S1 + (beta - 1) * S2
}

# A quick optim fit is used only to center the contour grid.
negloglik_on_logscale <- function(par) {
  alpha <- exp(par[1])
  beta <- exp(par[2])
  -loglik_beta(alpha, beta)
}
fit <- optim(par = log(c(2, 50)), fn = negloglik_on_logscale, method = "BFGS")
alpha_center <- exp(fit$par[1])
beta_center <- exp(fit$par[2])

alpha_grid <- seq(max(0.2, alpha_center / 4), alpha_center * 3, length.out = 160)
beta_grid <- seq(max(0.2, beta_center / 4), beta_center * 3, length.out = 160)
LL <- outer(alpha_grid, beta_grid, Vectorize(loglik_beta))

summary_table <- data.frame(
  n_player_seasons = n,
  min_HRrate = min(HRrate),
  q25_HRrate = as.numeric(quantile(HRrate, 0.25)),
  median_HRrate = median(HRrate),
  mean_HRrate = mean(HRrate),
  q75_HRrate = as.numeric(quantile(HRrate, 0.75)),
  max_HRrate = max(HRrate),
  sample_variance = var(HRrate),
  contour_center_alpha = alpha_center,
  contour_center_beta = beta_center,
  contour_center_logLik = loglik_beta(alpha_center, beta_center),
  row.names = NULL
)

keep_cols <- intersect(c("playerID", "yearID", "stint", "teamID", "lgID", "AB", "HR", "HRrate"),
                       names(Batting.restricted))
home_run_rates <- Batting.restricted[, keep_cols]

contour_grid <- data.frame(
  alpha = rep(alpha_grid, times = length(beta_grid)),
  beta = rep(beta_grid, each = length(alpha_grid)),
  log_likelihood = as.vector(LL)
)

write.csv(summary_table, file.path(results_dir, "home_run_rate_summary.csv"), row.names = FALSE)
write.csv(home_run_rates, file.path(results_dir, "restricted_home_run_rates.csv"), row.names = FALSE)
write.csv(contour_grid, file.path(results_dir, "beta_loglik_contour_grid.csv"), row.names = FALSE)

png(file.path(fig_dir, "home_run_rate_histogram.png"), width = 900, height = 700)
hist(HRrate, breaks = 40, main = "Histogram of Home Run Rate (HR/AB)",
     xlab = "HR/AB", ylab = "Count")
dev.off()

png(file.path(fig_dir, "beta_loglik_contour.png"), width = 900, height = 700)
contour(alpha_grid, beta_grid, LL,
        nlevels = 25,
        xlab = expression(alpha),
        ylab = expression(beta),
        main = expression("Beta log-likelihood " * ell(alpha, beta)))
points(alpha_center, beta_center, pch = 19)
text(alpha_center, beta_center, labels = "optim center", pos = 4)
dev.off()

cat("Analysis 02 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")
