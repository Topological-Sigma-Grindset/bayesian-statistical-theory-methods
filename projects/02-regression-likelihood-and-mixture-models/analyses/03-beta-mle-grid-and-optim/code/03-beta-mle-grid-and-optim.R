# Project 02 - Analysis 03
# Beta MLE by grid search and optim
#
# Required package: Lahman

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/03-beta-mle-grid-and-optim"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("Lahman", quietly = TRUE)) {
  stop("Analysis 03 requires the Lahman package. Install it with install.packages(\"Lahman\").")
}

data("Batting", package = "Lahman")

valid <- rep(TRUE, nrow(Batting))
valid[is.na(Batting$yearID) | Batting$yearID < 1970] <- FALSE
valid[is.na(Batting$HR) | Batting$HR < 1] <- FALSE
valid[is.na(Batting$AB) | Batting$AB < 200] <- FALSE

Batting.restricted <- Batting[valid, ]
Batting.restricted$HRrate <- Batting.restricted$HR / Batting.restricted$AB
HRrate <- Batting.restricted$HRrate
HRrate <- HRrate[is.finite(HRrate) & HRrate > 0 & HRrate < 1]

n <- length(HRrate)
S1 <- sum(log(HRrate))
S2 <- sum(log1p(-HRrate))

loglik_beta <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) return(-Inf)
  -n * lbeta(alpha, beta) + (alpha - 1) * S1 + (beta - 1) * S2
}

# Two-stage grid search.
alpha_grid1 <- seq(0.5, 10, length.out = 200)
beta_grid1 <- seq(5, 200, length.out = 200)
LL1 <- outer(alpha_grid1, beta_grid1, Vectorize(loglik_beta))
idx1 <- arrayInd(which.max(LL1), dim(LL1))
alpha_hat_stage1 <- alpha_grid1[idx1[1]]
beta_hat_stage1 <- beta_grid1[idx1[2]]

alpha_grid2 <- seq(max(0.1, alpha_hat_stage1 * 0.5), alpha_hat_stage1 * 1.5, length.out = 200)
beta_grid2 <- seq(max(0.1, beta_hat_stage1 * 0.5), beta_hat_stage1 * 1.5, length.out = 200)
LL2 <- outer(alpha_grid2, beta_grid2, Vectorize(loglik_beta))
idx2 <- arrayInd(which.max(LL2), dim(LL2))
alpha_hat_grid <- alpha_grid2[idx2[1]]
beta_hat_grid <- beta_grid2[idx2[2]]
logLik_grid <- loglik_beta(alpha_hat_grid, beta_hat_grid)

# Numerical MLE using log alpha and log beta to enforce positivity.
negloglik_uv <- function(par) {
  alpha <- exp(par[1])
  beta <- exp(par[2])
  -loglik_beta(alpha, beta)
}
fit <- optim(par = log(c(2, 50)), fn = negloglik_uv, method = "BFGS")
alpha_hat_optim <- exp(fit$par[1])
beta_hat_optim <- exp(fit$par[2])
logLik_optim <- loglik_beta(alpha_hat_optim, beta_hat_optim)

beta_variance <- function(alpha, beta) {
  (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
}

empirical_variance <- var(HRrate)

comparison <- data.frame(
  method = c("grid_search", "optim_logscale"),
  alpha_hat = c(alpha_hat_grid, alpha_hat_optim),
  beta_hat = c(beta_hat_grid, beta_hat_optim),
  log_likelihood = c(logLik_grid, logLik_optim),
  model_implied_variance = c(beta_variance(alpha_hat_grid, beta_hat_grid),
                             beta_variance(alpha_hat_optim, beta_hat_optim)),
  empirical_variance = empirical_variance,
  delta_alpha_vs_grid = c(0, alpha_hat_optim - alpha_hat_grid),
  delta_beta_vs_grid = c(0, beta_hat_optim - beta_hat_grid),
  row.names = NULL
)

grid_details <- data.frame(
  stage = c("coarse", "refined"),
  alpha_hat = c(alpha_hat_stage1, alpha_hat_grid),
  beta_hat = c(beta_hat_stage1, beta_hat_grid),
  log_likelihood = c(loglik_beta(alpha_hat_stage1, beta_hat_stage1), logLik_grid),
  alpha_min = c(min(alpha_grid1), min(alpha_grid2)),
  alpha_max = c(max(alpha_grid1), max(alpha_grid2)),
  beta_min = c(min(beta_grid1), min(beta_grid2)),
  beta_max = c(max(beta_grid1), max(beta_grid2)),
  grid_points = c(length(alpha_grid1) * length(beta_grid1), length(alpha_grid2) * length(beta_grid2)),
  row.names = NULL
)

write.csv(comparison, file.path(results_dir, "beta_mle_comparison.csv"), row.names = FALSE)
write.csv(grid_details, file.path(results_dir, "beta_grid_search_details.csv"), row.names = FALSE)

png(file.path(fig_dir, "beta_mle_refined_loglik_contour.png"), width = 900, height = 700)
contour(alpha_grid2, beta_grid2, LL2,
        nlevels = 25,
        xlab = expression(alpha),
        ylab = expression(beta),
        main = expression("Refined beta log-likelihood " * ell(alpha, beta)))
points(alpha_hat_grid, beta_hat_grid, pch = 19)
text(alpha_hat_grid, beta_hat_grid, labels = "grid", pos = 4)
points(alpha_hat_optim, beta_hat_optim, pch = 4, lwd = 2)
text(alpha_hat_optim, beta_hat_optim, labels = "optim", pos = 2)
dev.off()

cat("Analysis 03 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")
