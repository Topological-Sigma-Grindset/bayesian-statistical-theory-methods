# Project 02 runner
# Source this file from the repository root:
# source("projects/02-regression-likelihood-and-mixture-models/run_all.R")

project_dir <- "projects/02-regression-likelihood-and-mixture-models"

scripts <- c(
  file.path(project_dir, "analyses/01-radon-regression-posterior-prediction/code/01-radon-regression-posterior-prediction.R"),
  file.path(project_dir, "analyses/02-beta-model-home-run-rates/code/02-beta-model-home-run-rates.R"),
  file.path(project_dir, "analyses/03-beta-mle-grid-and-optim/code/03-beta-mle-grid-and-optim.R"),
  file.path(project_dir, "analyses/04-normal-halfnormal-mixture-em/code/04-normal-halfnormal-mixture-em.R"),
  file.path(project_dir, "analyses/05-ortiz-elite-probabilities-halfnormal-model/code/05-ortiz-elite-probabilities-halfnormal-model.R"),
  file.path(project_dir, "analyses/06-normal-beta-mixture-em/code/06-normal-beta-mixture-em.R"),
  file.path(project_dir, "analyses/07-ortiz-elite-probabilities-beta-mixture/code/07-ortiz-elite-probabilities-beta-mixture.R")
)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing analysis scripts:\n", paste(missing_scripts, collapse = "\n"))
}

for (script in scripts) {
  cat("\n--- Running ", script, " ---\n", sep = "")
  source(script)
}

cat("\nProject 02 complete.\n")
