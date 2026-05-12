# Project 01 / Analysis 06
# Baseball Inning Simulation
#
# Purpose:
#   Simulate the distribution of runs scored in a simplified baseball inning.
#   Compare the standard three-out inning with a hypothetical four-out inning.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/06-baseball-inning-simulation"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(927001)

# -----------------------------
# Inning simulator
# -----------------------------

simulate_inning <- function(
  n_outs = 3,
  outcome_prob = c(out = 0.6, single = 0.2, double = 0.1, hr = 0.1)
) {
  outs <- 0
  runs <- 0

  # Base states: TRUE means a runner occupies that base.
  b1 <- FALSE
  b2 <- FALSE
  b3 <- FALSE

  while (outs < n_outs) {
    outcome <- sample(names(outcome_prob), size = 1, prob = outcome_prob)

    if (outcome == "out") {
      outs <- outs + 1

    } else if (outcome == "single") {
      if (b3) runs <- runs + 1

      b3_new <- b2
      b2_new <- b1
      b1_new <- TRUE

      b1 <- b1_new
      b2 <- b2_new
      b3 <- b3_new

    } else if (outcome == "double") {
      if (b2) runs <- runs + 1
      if (b3) runs <- runs + 1

      b3_new <- b1
      b2_new <- TRUE
      b1_new <- FALSE

      b1 <- b1_new
      b2 <- b2_new
      b3 <- b3_new

    } else if (outcome == "hr") {
      runs <- runs + 1 + as.integer(b1) + as.integer(b2) + as.integer(b3)

      b1 <- FALSE
      b2 <- FALSE
      b3 <- FALSE
    }
  }

  runs
}

summarize_runs <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    max = max(x)
  )
}

# -----------------------------
# Simulation
# -----------------------------

n_sim <- 1000
runs_3_outs <- replicate(n_sim, simulate_inning(n_outs = 3))
runs_4_outs <- replicate(n_sim, simulate_inning(n_outs = 4))

simulation_draws <- data.frame(
  simulation = seq_len(n_sim),
  runs_3_outs = runs_3_outs,
  runs_4_outs = runs_4_outs
)

write.csv(
  simulation_draws,
  file = file.path(results_dir, "baseball_simulation_draws.csv"),
  row.names = FALSE
)

summary_table <- rbind(
  data.frame(inning_type = "3_outs", t(summarize_runs(runs_3_outs))),
  data.frame(inning_type = "4_outs", t(summarize_runs(runs_4_outs)))
)

write.csv(
  summary_table,
  file = file.path(results_dir, "baseball_simulation_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Plots
# -----------------------------

png(file.path(fig_dir, "01_runs_3_out_inning.png"), width = 1200, height = 800, res = 150)
hist(
  runs_3_outs,
  breaks = seq(-0.5, max(runs_3_outs) + 0.5, by = 1),
  xlab = "Runs",
  ylab = "Frequency",
  main = "Runs scored in a 3-out inning"
)
dev.off()

png(file.path(fig_dir, "02_runs_4_out_inning.png"), width = 1200, height = 800, res = 150)
hist(
  runs_4_outs,
  breaks = seq(-0.5, max(runs_4_outs) + 0.5, by = 1),
  xlab = "Runs",
  ylab = "Frequency",
  main = "Runs scored in a 4-out inning"
)
dev.off()

png(file.path(fig_dir, "03_runs_comparison_boxplot.png"), width = 1200, height = 800, res = 150)
boxplot(
  runs_3_outs,
  runs_4_outs,
  names = c("3 outs", "4 outs"),
  ylab = "Runs",
  main = "Comparison of simulated runs by inning length"
)
dev.off()

cat("Analysis 06 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")