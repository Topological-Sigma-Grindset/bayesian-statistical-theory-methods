# Project 02 - Analysis 06
# EM for normal / beta mixture of Lahman home-run rates
#
# Required package: Lahman
#
# Testing-speed version:
# - n_starts set to to 5
# - maxit set to 200

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/06-normal-beta-mixture-em"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("Lahman", quietly = TRUE)) {
  stop("Analysis 06 requires the Lahman package. Install it with install.packages(\"Lahman\").")
}

set.seed(927006)
data("Batting", package = "Lahman")

valid <- rep(TRUE, nrow(Batting))
valid[is.na(Batting$yearID) | Batting$yearID < 1970] <- FALSE
valid[is.na(Batting$HR) | Batting$HR < 1] <- FALSE
valid[is.na(Batting$AB) | Batting$AB < 200] <- FALSE

Batting.restricted <- Batting[valid, ]
Batting.restricted$HRrate <- Batting.restricted$HR / Batting.restricted$AB
Batting.restricted <- Batting.restricted[
  is.finite(Batting.restricted$HRrate) &
    Batting.restricted$HRrate > 0 &
    Batting.restricted$HRrate < 1,
]

x <- Batting.restricted$HRrate
n <- length(x)

logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

loglik_normal_beta_mix <- function(pi, mu, sigma, alpha, beta, x) {
  if (pi <= 0 || pi >= 1 || sigma <= 0 || alpha <= 0 || beta <= 0) {
    return(-Inf)
  }

  lf1 <- dnorm(x, mean = mu, sd = sigma, log = TRUE)
  lf0 <- dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE)

  sum(logsumexp2(log(pi) + lf1, log1p(-pi) + lf0))
}

elite_prob_normal_beta <- function(x, pi, mu, sigma, alpha, beta) {
  lf1 <- dnorm(x, mean = mu, sd = sigma, log = TRUE)
  lf0 <- dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE)

  lognum1 <- log(pi) + lf1
  lognum0 <- log1p(-pi) + lf0

  delta <- pmin(700, lognum0 - lognum1)

  1 / (1 + exp(delta))
}

run_em_normal_beta <- function(x, pi, mu, sigma, alpha, beta,
                               maxit = 200,
                               tol = 1e-6,
                               pi_min = 1e-4,
                               sigma_min = 1e-6) {
  ll_old <- loglik_normal_beta_mix(pi, mu, sigma, alpha, beta, x)

  for (it in seq_len(maxit)) {
    r <- elite_prob_normal_beta(x, pi, mu, sigma, alpha, beta)

    pi_new <- mean(r)
    pi_new <- min(max(pi_new, pi_min), 1 - pi_min)

    R <- sum(r)

    if (R < .Machine$double.eps) {
      R <- .Machine$double.eps
    }

    mu_new <- sum(r * x) / R

    sigma2_new <- sum(r * (x - mu_new)^2) / R
    sigma2_new <- max(sigma2_new, sigma_min^2)
    sigma_new <- sqrt(sigma2_new)

    w <- 1 - r
    W <- sum(w)

    if (W < 1e-8) {
      alpha_new <- alpha
      beta_new <- beta
    } else {
      S1 <- sum(w * log(x))
      S2 <- sum(w * log1p(-x))

      negQ_ab <- function(par) {
        a <- exp(par[1])
        b <- exp(par[2])

        Q <- (a - 1) * S1 + (b - 1) * S2 - W * lbeta(a, b)

        -Q
      }

      opt <- tryCatch(
        optim(
          par = log(c(alpha, beta)),
          fn = negQ_ab,
          method = "L-BFGS-B",
          lower = log(c(1e-4, 1e-4)),
          upper = log(c(1e6, 1e6))
        ),
        error = function(e) NULL
      )

      if (!is.null(opt) && is.finite(opt$value)) {
        alpha_new <- exp(opt$par[1])
        beta_new <- exp(opt$par[2])
      } else {
        alpha_new <- alpha
        beta_new <- beta
      }
    }

    ll_new <- loglik_normal_beta_mix(
      pi_new,
      mu_new,
      sigma_new,
      alpha_new,
      beta_new,
      x
    )

    if (abs(ll_new - ll_old) < tol) {
      return(
        list(
          pi = pi_new,
          mu = mu_new,
          sigma = sigma_new,
          alpha = alpha_new,
          beta = beta_new,
          ll = ll_new,
          it = it,
          conv = TRUE
        )
      )
    }

    pi <- pi_new
    mu <- mu_new
    sigma <- sigma_new
    alpha <- alpha_new
    beta <- beta_new
    ll_old <- ll_new
  }

  list(
    pi = pi,
    mu = mu,
    sigma = sigma,
    alpha = alpha,
    beta = beta,
    ll = ll_old,
    it = maxit,
    conv = FALSE
  )
}

m <- mean(x)
v <- var(x)
tmp <- m * (1 - m) / v - 1

if (!is.finite(tmp) || tmp <= 0) {
  tmp <- 100
}

alpha_mom <- max(0.5, m * tmp)
beta_mom <- max(0.5, (1 - m) * tmp)

# Reduced for faster local testing.
# For the final full version, you may increase this back to 30.
n_starts <- 5
fits <- vector("list", n_starts)

cat("Analysis 06: fitting Normal/Beta mixture by multi-start EM.\n")
cat("Number of observations: ", n, "\n", sep = "")
cat("Number of starts: ", n_starts, "\n", sep = "")

for (k in seq_len(n_starts)) {
  cat("Running start ", k, " of ", n_starts, "...\n", sep = "")

  pi0 <- runif(1, 0.02, 0.20)
  mu0 <- as.numeric(quantile(x, runif(1, 0.85, 0.995)))
  sigma0 <- sd(x)

  alpha0 <- max(1e-3, alpha_mom * exp(rnorm(1, 0, 0.4)))
  beta0 <- max(1e-3, beta_mom * exp(rnorm(1, 0, 0.4)))

  fits[[k]] <- run_em_normal_beta(
    x = x,
    pi = pi0,
    mu = mu0,
    sigma = sigma0,
    alpha = alpha0,
    beta = beta0,
    maxit = 200,
    tol = 1e-6
  )

  cat(
    "  start ", k,
    " complete: logLik = ", round(fits[[k]]$ll, 3),
    ", iterations = ", fits[[k]]$it,
    ", converged = ", fits[[k]]$conv,
    "\n",
    sep = ""
  )
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
  alpha_hat = sapply(fits, function(f) f$alpha),
  beta_hat = sapply(fits, function(f) f$beta),
  log_likelihood = lls,
  iterations = sapply(fits, function(f) f$it),
  converged = sapply(fits, function(f) f$conv),
  row.names = NULL
)

start_results <- start_results[order(-start_results$log_likelihood), ]

best_parameters <- data.frame(
  model = "normal_beta_mixture",
  n_observations = n,
  n_starts = n_starts,
  best_start = best_idx,
  pi_hat = best$pi,
  mu_hat = best$mu,
  sigma_hat = best$sigma,
  sigma2_hat = best$sigma^2,
  alpha_hat = best$alpha,
  beta_hat = best$beta,
  log_likelihood = best$ll,
  iterations = best$it,
  converged = best$conv,
  row.names = NULL
)

Batting.restricted$p_elite <- elite_prob_normal_beta(
  x,
  best$pi,
  best$mu,
  best$sigma,
  best$alpha,
  best$beta
)

keep_cols <- intersect(
  c(
    "playerID",
    "yearID",
    "stint",
    "teamID",
    "lgID",
    "AB",
    "HR",
    "HRrate",
    "p_elite"
  ),
  names(Batting.restricted)
)

responsibilities <- Batting.restricted[, keep_cols]

write.csv(
  start_results,
  file.path(results_dir, "em_start_results.csv"),
  row.names = FALSE
)

write.csv(
  best_parameters,
  file.path(results_dir, "best_fit_parameters.csv"),
  row.names = FALSE
)

write.csv(
  responsibilities,
  file.path(results_dir, "home_run_rate_elite_probabilities.csv"),
  row.names = FALSE
)

xg <- seq(0, max(x), length.out = 600)

f_elite <- dnorm(xg, mean = best$mu, sd = best$sigma)
f_non <- dbeta(xg, shape1 = best$alpha, shape2 = best$beta)
f_mix <- best$pi * f_elite + (1 - best$pi) * f_non

png(
  file.path(fig_dir, "normal_beta_mixture_density.png"),
  width = 900,
  height = 700
)

hist(
  x,
  breaks = 50,
  freq = FALSE,
  main = "HR/AB with fitted Normal/Beta mixture densities",
  xlab = "HR/AB",
  ylab = "Density"
)

lines(xg, f_mix, lwd = 2)
lines(xg, best$pi * f_elite, lwd = 2, lty = 2)
lines(xg, (1 - best$pi) * f_non, lwd = 2, lty = 3)

legend(
  "topright",
  legend = c(
    "Mixture",
    "Elite (weighted Normal)",
    "Non-elite (weighted Beta)"
  ),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)

mtext(
  sprintf(
    "pi=%.4f, mu=%.4f, sigma=%.4f, alpha=%.2f, beta=%.2f, logLik=%.2f",
    best$pi,
    best$mu,
    best$sigma,
    best$alpha,
    best$beta,
    best$ll
  ),
  side = 3,
  line = 0.2,
  cex = 0.80
)

dev.off()

cat("Analysis 06 complete.\n")
cat("Best start: ", best_idx, "\n", sep = "")
cat("Best log-likelihood: ", round(best$ll, 3), "\n", sep = "")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")