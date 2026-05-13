# Project 02 - Analysis 05
# David Ortiz elite probabilities under the normal / half-normal mixture
#
# Required package: Lahman

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/05-ortiz-elite-probabilities-halfnormal-model"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("Lahman", quietly = TRUE)) {
  stop("Analysis 05 requires the Lahman package. Install it with install.packages(\"Lahman\").")
}

set.seed(927005)
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

fit_used <- data.frame(
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

ortiz <- Batting.restricted[Batting.restricted$playerID == "ortizda01", ]
ortiz <- ortiz[order(ortiz$yearID), ]
ortiz$p_elite <- elite_prob_halfnormal(ortiz$HRrate, best$pi, best$mu, best$sigma)

ortiz_out <- ortiz[, intersect(c("playerID", "yearID", "stint", "teamID", "lgID", "AB", "HR", "HRrate", "p_elite"), names(ortiz))]

write.csv(start_results, file.path(results_dir, "em_start_results_halfnormal.csv"), row.names = FALSE)
write.csv(fit_used, file.path(results_dir, "halfnormal_fit_parameters_used.csv"), row.names = FALSE)
write.csv(ortiz_out, file.path(results_dir, "ortiz_elite_probabilities_halfnormal.csv"), row.names = FALSE)

png(file.path(fig_dir, "ortiz_elite_probabilities_halfnormal.png"), width = 900, height = 700)
plot(ortiz$yearID, ortiz$p_elite, type = "b", pch = 19,
     ylim = c(0, 1), xlab = "Season", ylab = "P(elite | HR/AB)",
     main = "David Ortiz elite probabilities: Normal/Half-Normal mixture")
abline(h = c(0.5, 0.9), lty = c(2, 3))
dev.off()

cat("Analysis 05 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")
