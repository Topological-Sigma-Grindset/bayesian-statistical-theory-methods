# Project 02 - Analysis 04
# EM for normal / half-normal mixture of Lahman home-run rates
#
# Required package: Lahman

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/04-normal-halfnormal-mixture-em"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("Lahman", quietly = TRUE)) {
  stop("Analysis 04 requires the Lahman package. Install it with install.packages(\"Lahman\").")
}

set.seed(927004)
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

x <- Batting.restricted$HRrate
n <- length(x)

logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

d_halfnorm <- function(x, sigma) {
  out <- rep(0, length(x))
  ok <- x >= 0
  out[ok] <- 2 * dnorm(x[ok], mean = 0, sd = sigma)
  out
}

log_d_halfnorm <- function(x, sigma) {
  out <- rep(-Inf, length(x))
  ok <- x >= 0
  out[ok] <- log(2) + dnorm(x[ok], mean = 0, sd = sigma, log = TRUE)
  out
}

loglik_halfnormal_mix <- function(pi, mu, sigma, x) {
  if (pi <= 0 || pi >= 1 || sigma <= 0) return(-Inf)
  lf1 <- dnorm(x, mean = mu, sd = sigma, log = TRUE)
  lf0 <- log_d_halfnorm(x, sigma)
  sum(logsumexp2(log(pi) + lf1, log1p(-pi) + lf0))
}

elite_prob_halfnormal <- function(x, pi, mu, sigma) {
  lf1 <- dnorm(x, mean = mu, sd = sigma, log = TRUE)
  lf0 <- log_d_halfnorm(x, sigma)
  lognum1 <- log(pi) + lf1
  lognum0 <- log1p(-pi) + lf0
  delta <- pmin(700, lognum0 - lognum1)
  1 / (1 + exp(delta))
}

run_em_halfnormal <- function(x, pi, mu, sigma,
                              maxit = 2000, tol = 1e-8,
                              sigma_min = 0.008, pi_min = 1e-4) {
  n <- length(x)
  ll_old <- loglik_halfnormal_mix(pi, mu, sigma, x)

  for (it in seq_len(maxit)) {
    r <- elite_prob_halfnormal(x, pi, mu, sigma)

    pi_new <- mean(r)
    pi_new <- min(max(pi_new, pi_min), 1 - pi_min)

    sr <- sum(r)
    if (sr < .Machine$double.eps) sr <- .Machine$double.eps
    mu_new <- sum(r * x) / sr

    sigma2_new <- (1 / n) * sum(r * (x - mu_new)^2 + (1 - r) * x^2)
    sigma2_new <- max(sigma2_new, sigma_min^2)
    sigma_new <- sqrt(sigma2_new)

    ll_new <- loglik_halfnormal_mix(pi_new, mu_new, sigma_new, x)

    if (abs(ll_new - ll_old) < tol) {
      return(list(pi = pi_new, mu = mu_new, sigma = sigma_new,
                  ll = ll_new, it = it, conv = TRUE))
    }

    pi <- pi_new
    mu <- mu_new
    sigma <- sigma_new
    ll_old <- ll_new
  }

  list(pi = pi, mu = mu, sigma = sigma, ll = ll_old, it = maxit, conv = FALSE)
}

n_starts <- 30
fits <- vector("list", n_starts)
for (k in seq_len(n_starts)) {
  pi0 <- runif(1, 0.02, 0.20)
  mu0 <- as.numeric(quantile(x, runif(1, 0.85, 0.995)))
  sigma0 <- sd(x)
  fits[[k]] <- run_em_halfnormal(x, pi0, mu0, sigma0)
}

lls <- sapply(fits, function(f) f$ll)
best_idx <- which.max(lls)
best <- fits[[best_idx]]

start_results <- data.frame(
  start = seq_len(n_starts),
  pi_hat = sapply(fits, function(f) f$pi),
  mu_hat = sapply(fits, function(f) f$mu),
  sigma_hat = sapply(fits, function(f) f$sigma),
  sigma2_hat = sapply(fits, function(f) f$sigma^2),
  log_likelihood = lls,
  iterations = sapply(fits, function(f) f$it),
  converged = sapply(fits, function(f) f$conv),
  row.names = NULL
)
start_results <- start_results[order(-start_results$log_likelihood), ]

best_parameters <- data.frame(
  model = "normal_halfnormal_mixture",
  n_observations = n,
  n_starts = n_starts,
  best_start = best_idx,
  pi_hat = best$pi,
  mu_hat = best$mu,
  sigma_hat = best$sigma,
  sigma2_hat = best$sigma^2,
  log_likelihood = best$ll,
  iterations = best$it,
  converged = best$conv,
  row.names = NULL
)

Batting.restricted$p_elite <- elite_prob_halfnormal(x, best$pi, best$mu, best$sigma)
keep_cols <- intersect(c("playerID", "yearID", "stint", "teamID", "lgID", "AB", "HR", "HRrate", "p_elite"),
                       names(Batting.restricted))
responsibilities <- Batting.restricted[, keep_cols]

write.csv(start_results, file.path(results_dir, "em_start_results.csv"), row.names = FALSE)
write.csv(best_parameters, file.path(results_dir, "best_fit_parameters.csv"), row.names = FALSE)
write.csv(responsibilities, file.path(results_dir, "home_run_rate_elite_probabilities.csv"), row.names = FALSE)

xg <- seq(0, max(x), length.out = 600)
f_elite <- dnorm(xg, mean = best$mu, sd = best$sigma)
f_non <- d_halfnorm(xg, best$sigma)
f_mix <- best$pi * f_elite + (1 - best$pi) * f_non

png(file.path(fig_dir, "normal_halfnormal_mixture_density.png"), width = 900, height = 700)
hist(x, breaks = 50, freq = FALSE,
     main = "HR/AB with fitted Normal/Half-Normal mixture densities",
     xlab = "HR/AB", ylab = "Density")
lines(xg, f_mix, lwd = 2)
lines(xg, best$pi * f_elite, lwd = 2, lty = 2)
lines(xg, (1 - best$pi) * f_non, lwd = 2, lty = 3)
legend("topright",
       legend = c("Mixture", "Elite (weighted Normal)", "Non-elite (weighted Half-Normal)"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")
mtext(sprintf("pi=%.4f, mu=%.4f, sigma=%.4f, logLik=%.2f",
              best$pi, best$mu, best$sigma, best$ll),
      side = 3, line = 0.2, cex = 0.85)
dev.off()

cat("Analysis 04 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")
