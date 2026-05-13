# Project 02 - Analysis 01
# Radon regression and posterior prediction
#
# Required packages: base R only
# Required local data: radon CSV in this analysis data/ folder

out_dir <- "projects/02-regression-likelihood-and-mixture-models/analyses/01-radon-regression-posterior-prediction"
fig_dir <- file.path(out_dir, "figures")
results_dir <- file.path(out_dir, "results")
data_dir <- file.path(out_dir, "data")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(927001)

candidate_files <- c(
  file.path(data_dir, "radon.table.7.3.csv"),
  file.path(data_dir, "radon_table_7_3.csv"),
  file.path(data_dir, "table_7_3_radon.csv"),
  file.path(data_dir, "radon.csv")
)

existing_candidates <- candidate_files[file.exists(candidate_files)]
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(existing_candidates) > 0) {
  radon_file <- existing_candidates[1]
} else if (length(csv_files) == 1) {
  radon_file <- csv_files[1]
} else {
  radon_file <- NA_character_
}

if (is.na(radon_file)) {
  requirement_text <- c(
    "Analysis 01 requires a local radon CSV that is not stored in this public repository.",
    "Place the local radon file in this analysis data directory:",
    data_dir,
    "Accepted filenames:",
    "  - radon.table.7.3.csv",
    "  - radon_table_7_3.csv",
    "  - table_7_3_radon.csv",
    "  - radon.csv",
    "If exactly one CSV is present in the data folder, the script will use that file.",
    "Expected variables are radon measurement, county, and floor."
  )

  writeLines(requirement_text, file.path(results_dir, "radon_data_requirement.txt"))

  cat("Analysis 01 complete: local radon CSV not found, so calculation was skipped.\n")
  cat("Requirement note saved to: ", file.path(results_dir, "radon_data_requirement.txt"), "\n", sep = "")
  cat("Figures saved to: ", fig_dir, "\n", sep = "")
  cat("Results saved to: ", results_dir, "\n", sep = "")

} else {

  raw_dat <- read.csv(radon_file, stringsAsFactors = FALSE)

  normalize_name <- function(x) {
    tolower(gsub("[^a-z0-9]", "", x))
  }

  find_col <- function(possible_names, dat, label) {
    actual_names <- names(dat)
    actual_normalized <- normalize_name(actual_names)
    possible_normalized <- normalize_name(possible_names)

    match_index <- match(possible_normalized, actual_normalized)
    match_index <- match_index[!is.na(match_index)]

    if (length(match_index) == 0) {
      stop(
        "Could not find ", label, " column in the radon CSV.\n",
        "Found columns: ", paste(actual_names, collapse = ", "),
        call. = FALSE
      )
    }

    actual_names[match_index[1]]
  }

  radon_col <- find_col(
    c(
      "radon",
      "Radon",
      "radon measurement",
      "radonmeasurement",
      "radon measure",
      "radonmeasure",
      "activity"
    ),
    raw_dat,
    "a radon measurement"
  )

  county_col <- find_col(
    c(
      "county",
      "County",
      "county name",
      "countyname"
    ),
    raw_dat,
    "a county"
  )

  floor_col <- find_col(
    c(
      "floor",
      "Floor",
      "measurement floor",
      "measurementfloor",
      "story"
    ),
    raw_dat,
    "a floor"
  )

  dat <- data.frame(
    Radon = suppressWarnings(as.numeric(raw_dat[[radon_col]])),
    County = raw_dat[[county_col]],
    Floor = raw_dat[[floor_col]],
    stringsAsFactors = FALSE
  )

  dat <- dat[is.finite(dat$Radon) & !is.na(dat$County) & !is.na(dat$Floor), ]

  if (nrow(dat) == 0) {
    stop("No usable radon rows remain after removing missing values.")
  }

  if (any(dat$Radon <= 0)) {
    stop("All radon measurements must be positive because the model uses log(Radon).")
  }

  # County indicators. Relevel Blue Earth / BlueEarth as baseline if present.
  dat$County <- factor(trimws(as.character(dat$County)))
  county_levels <- levels(dat$County)

  blue_earth_level <- county_levels[
    grepl("^blue\\s*earth$|^blueearth$", county_levels, ignore.case = TRUE)
  ]

  if (length(blue_earth_level) > 0) {
    dat$County <- relevel(dat$County, ref = blue_earth_level[1])
  }

  # Floor indicator. Accept common string labels and 0/1 coding.
  if (is.numeric(dat$Floor) || is.integer(dat$Floor)) {
    floor_label <- ifelse(
      dat$Floor == 0,
      "Basement",
      ifelse(dat$Floor == 1, "First", as.character(dat$Floor))
    )
  } else {
    floor_text <- trimws(as.character(dat$Floor))
    floor_lower <- tolower(floor_text)

    floor_label <- ifelse(
      grepl("base|basement|^0$", floor_lower),
      "Basement",
      ifelse(grepl("first|1st|^1$", floor_lower), "First", floor_text)
    )
  }

  dat$Floor <- factor(floor_label)
  floor_levels <- levels(dat$Floor)

  basement_level <- floor_levels[
    grepl("base", floor_levels, ignore.case = TRUE)
  ]

  if (length(basement_level) > 0) {
    dat$Floor <- relevel(dat$Floor, ref = basement_level[1])
  }

  if (length(levels(dat$Floor)) < 2) {
    stop("The floor variable must contain at least basement and first-floor measurements.")
  }

  first_floor_level <- setdiff(levels(dat$Floor), levels(dat$Floor)[1])[1]

  # Bayesian linear regression with p(beta, sigma^2) proportional to 1 / sigma^2.
  y <- log(dat$Radon)
  X <- model.matrix(~ County + Floor, data = dat)

  n <- nrow(X)
  p <- ncol(X)
  df <- n - p

  if (df <= 2) {
    stop("The regression requires more observations than model parameters, with df > 2 for sigma summaries.")
  }

  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% y
  resid <- as.vector(y - X %*% beta_hat)
  SSE <- sum(resid^2)
  s2 <- SSE / df

  # Marginal posterior intervals for beta_j under the Student-t form.
  tcrit <- qt(0.975, df = df)
  se_beta <- sqrt(diag(XtX_inv) * s2)
  beta_lower <- as.vector(beta_hat) - tcrit * se_beta
  beta_upper <- as.vector(beta_hat) + tcrit * se_beta

  beta_table <- data.frame(
    parameter = colnames(X),
    posterior_mean = as.vector(beta_hat),
    posterior_sd_approx = se_beta,
    ci_2.5 = beta_lower,
    ci_97.5 = beta_upper,
    row.names = NULL
  )

  sigma2_mean <- (df * s2) / (df - 2)

  sigma2_ci <- c(
    (df * s2) / qchisq(0.975, df = df),
    (df * s2) / qchisq(0.025, df = df)
  )

  sigma_ci <- sqrt(sigma2_ci)

  sigma_table <- data.frame(
    parameter = c("sigma2", "sigma"),
    posterior_summary = c(sigma2_mean, sqrt(sigma2_mean)),
    ci_2.5 = c(sigma2_ci[1], sigma_ci[1]),
    ci_97.5 = c(sigma2_ci[2], sigma_ci[2]),
    row.names = NULL
  )

  # Posterior predictive simulation for a new Clay County house.
  S <- 5000
  L <- chol(XtX_inv)

  county_levels <- levels(dat$County)

  clay_level <- county_levels[
    grepl("^clay$", county_levels, ignore.case = TRUE)
  ]

  if (length(clay_level) == 0) {
    stop("Could not find Clay County in the county factor levels: ", paste(county_levels, collapse = ", "))
  }

  clay_level <- clay_level[1]

  new_base <- data.frame(
    County = factor(clay_level, levels = levels(dat$County)),
    Floor = factor(levels(dat$Floor)[1], levels = levels(dat$Floor))
  )

  new_first <- data.frame(
    County = factor(clay_level, levels = levels(dat$County)),
    Floor = factor(first_floor_level, levels = levels(dat$Floor))
  )

  make_design_row <- function(new_data, reference_columns) {
    mm <- model.matrix(~ County + Floor, data = new_data)

    missing_cols <- setdiff(reference_columns, colnames(mm))

    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        mm <- cbind(mm, 0)
        colnames(mm)[ncol(mm)] <- col
      }
    }

    mm[, reference_columns, drop = FALSE]
  }

  x_base <- make_design_row(new_base, colnames(X))
  x_first <- make_design_row(new_first, colnames(X))

  radon_base <- numeric(S)
  radon_first <- numeric(S)

  for (s in seq_len(S)) {
    sigma2_draw <- (df * s2) / rchisq(1, df = df)
    sigma_draw <- sqrt(sigma2_draw)

    beta_draw <- as.vector(beta_hat) + sigma_draw * as.vector(t(L) %*% rnorm(p))

    y_base <- rnorm(
      1,
      mean = as.vector(x_base %*% beta_draw),
      sd = sigma_draw
    )

    y_first <- rnorm(
      1,
      mean = as.vector(x_first %*% beta_draw),
      sd = sigma_draw
    )

    radon_base[s] <- exp(y_base)
    radon_first[s] <- exp(y_first)
  }

  predictive_table <- data.frame(
    scenario = c("Clay County basement", "Clay County first floor"),
    posterior_predictive_mean = c(mean(radon_base), mean(radon_first)),
    posterior_predictive_median = c(median(radon_base), median(radon_first)),
    ci_2.5 = c(
      as.numeric(quantile(radon_base, 0.025)),
      as.numeric(quantile(radon_first, 0.025))
    ),
    ci_97.5 = c(
      as.numeric(quantile(radon_base, 0.975)),
      as.numeric(quantile(radon_first, 0.975))
    ),
    row.names = NULL
  )

  predictive_draws <- data.frame(
    draw = seq_len(S),
    clay_basement = radon_base,
    clay_first_floor = radon_first
  )

  data_summary <- data.frame(
    radon_file = radon_file,
    n_observations = nrow(dat),
    n_parameters = p,
    residual_df = df,
    county_levels = paste(levels(dat$County), collapse = "; "),
    floor_levels = paste(levels(dat$Floor), collapse = "; "),
    row.names = NULL
  )

  write.csv(
    data_summary,
    file.path(results_dir, "radon_data_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    beta_table,
    file.path(results_dir, "regression_coefficient_posterior_summaries.csv"),
    row.names = FALSE
  )

  write.csv(
    sigma_table,
    file.path(results_dir, "sigma_posterior_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    predictive_table,
    file.path(results_dir, "clay_county_posterior_predictive_intervals.csv"),
    row.names = FALSE
  )

  write.csv(
    predictive_draws,
    file.path(results_dir, "clay_county_posterior_predictive_draws.csv"),
    row.names = FALSE
  )

  png(
    file.path(fig_dir, "clay_county_posterior_predictive_radon.png"),
    width = 1200,
    height = 600
  )

  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  hist(
    radon_base,
    breaks = 40,
    main = "Clay County basement",
    xlab = "Predicted radon",
    ylab = "Count"
  )

  hist(
    radon_first,
    breaks = 40,
    main = "Clay County first floor",
    xlab = "Predicted radon",
    ylab = "Count"
  )

  par(op)
  dev.off()

  cat("Analysis 01 complete.\n")
  cat("Radon file used: ", radon_file, "\n", sep = "")
  cat("Figures saved to: ", fig_dir, "\n", sep = "")
  cat("Results saved to: ", results_dir, "\n", sep = "")
}