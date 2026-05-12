# Project 01 / Analysis 05
# Train Station Queue Simulation
#
# Purpose:
#   Simulate a four-hour ticket-counter queue with three agents, exponential
#   customer interarrival times, and uniform service times.
#
# How to run:
#   Open this file in VS Code and press Source.
#   The script assumes VS Code is opened at the repository root:
#   bayesian-statistical-theory-methods/

# -----------------------------
# Project setup
# -----------------------------

out_dir <- "projects/01-bayesian-updating-and-simulation-basics/analyses/05-train-station-queue-simulation"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(927001)

# -----------------------------
# Queue simulator
# -----------------------------

simulate_shift <- function(
  max_time = 240,
  mean_interarrival = 10,
  service_min = 8,
  service_max = 20,
  n_agents = 3
) {
  current_time <- 0
  arrivals <- numeric(0)

  while (TRUE) {
    current_time <- current_time + rexp(1, rate = 1 / mean_interarrival)
    if (current_time > max_time) break
    arrivals <- c(arrivals, current_time)
  }

  n_customers <- length(arrivals)

  if (n_customers == 0) {
    return(list(
      n_customers = 0L,
      n_wait = 0L,
      avg_wait_if_wait = 0,
      wait_times = numeric(0),
      arrivals = numeric(0)
    ))
  }

  free_time <- rep(0, n_agents)
  wait_times <- numeric(n_customers)

  for (i in seq_len(n_customers)) {
    arrival_time <- arrivals[i]
    agent <- which.min(free_time)
    service_start <- max(arrival_time, free_time[agent])
    wait_times[i] <- service_start - arrival_time

    service_time <- runif(1, min = service_min, max = service_max)
    free_time[agent] <- service_start + service_time
  }

  n_wait <- sum(wait_times > 0)
  avg_wait_if_wait <- if (n_wait == 0) 0 else mean(wait_times[wait_times > 0])

  list(
    n_customers = n_customers,
    n_wait = n_wait,
    avg_wait_if_wait = avg_wait_if_wait,
    wait_times = wait_times,
    arrivals = arrivals
  )
}

summarize_draws <- function(x) {
  c(
    median = median(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975))
  )
}

# -----------------------------
# One simulated shift
# -----------------------------

single_shift <- simulate_shift()

single_shift_summary <- data.frame(
  quantity = c("n_customers", "n_wait", "avg_wait_if_wait"),
  value = c(
    single_shift$n_customers,
    single_shift$n_wait,
    single_shift$avg_wait_if_wait
  )
)

write.csv(
  single_shift_summary,
  file = file.path(results_dir, "single_shift_summary.csv"),
  row.names = FALSE
)

if (single_shift$n_customers > 0) {
  single_shift_customer_level <- data.frame(
    customer = seq_len(single_shift$n_customers),
    arrival_time = single_shift$arrivals,
    wait_time = single_shift$wait_times
  )

  write.csv(
    single_shift_customer_level,
    file = file.path(results_dir, "single_shift_customer_level.csv"),
    row.names = FALSE
  )
}

# -----------------------------
# Monte Carlo simulation over 100 shifts
# -----------------------------

n_sim <- 100
n_customers <- numeric(n_sim)
n_wait <- numeric(n_sim)
avg_wait <- numeric(n_sim)

for (sim in seq_len(n_sim)) {
  shift <- simulate_shift()
  n_customers[sim] <- shift$n_customers
  n_wait[sim] <- shift$n_wait
  avg_wait[sim] <- shift$avg_wait_if_wait
}

simulation_draws <- data.frame(
  simulation = seq_len(n_sim),
  n_customers = n_customers,
  n_wait = n_wait,
  avg_wait_if_wait = avg_wait
)

write.csv(
  simulation_draws,
  file = file.path(results_dir, "queue_simulation_draws.csv"),
  row.names = FALSE
)

summary_table <- rbind(
  data.frame(quantity = "n_customers", t(summarize_draws(n_customers))),
  data.frame(quantity = "n_wait", t(summarize_draws(n_wait))),
  data.frame(quantity = "avg_wait_if_wait", t(summarize_draws(avg_wait)))
)

write.csv(
  summary_table,
  file = file.path(results_dir, "queue_simulation_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Plots
# -----------------------------

png(file.path(fig_dir, "01_customers_per_shift.png"), width = 1200, height = 800, res = 150)
hist(
  n_customers,
  breaks = 15,
  xlab = "Number of customers",
  main = "Customers arriving during a four-hour shift"
)
dev.off()

png(file.path(fig_dir, "02_customers_who_waited.png"), width = 1200, height = 800, res = 150)
hist(
  n_wait,
  breaks = 15,
  xlab = "Number of customers who waited",
  main = "Customers who waited before service"
)
dev.off()

png(file.path(fig_dir, "03_average_wait_if_waited.png"), width = 1200, height = 800, res = 150)
hist(
  avg_wait,
  breaks = 15,
  xlab = "Average waiting time among customers who waited",
  main = "Average waiting time across simulated shifts"
)
dev.off()

cat("Analysis 05 complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Results saved to: ", results_dir, "\n", sep = "")