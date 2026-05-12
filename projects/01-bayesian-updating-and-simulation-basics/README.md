# Project 01: Bayesian Updating and Simulation Basics

This project contains a set of self-contained Bayesian analyses and simulations. The analyses move from closed-form Bayesian updating to grid-based posterior approximation, posterior predictive simulation, and Monte Carlo simulation of stochastic systems.

Each analysis is written as an independent project with its own code, figures, results, and short explanation.

## Analyses

| Analysis | Topic | Main methods |
|---|---|---|
| 01 | Beta-binomial foundations | Prior predictive distribution, conjugacy, posterior variance |
| 02 | Informative beta prior for a survey proportion | Prior elicitation, beta-binomial updating, posterior summaries |
| 03 | Cauchy location model | Grid approximation, posterior sampling, posterior predictive simulation |
| 04 | Airline fatal-accident counts | Poisson trend model, two-dimensional grid posterior, posterior prediction |
| 05 | Train-station queue | Exponential arrivals, uniform service times, discrete-event simulation |
| 06 | Baseball inning scoring | State-based simulation, run distribution comparison |

## Folder structure

```text
analyses/
  01-beta-binomial-foundations/
    README.md
    code/
    figures/
    results/
    notes/
  02-informative-beta-prior/
  03-cauchy-location-grid-approximation/
  04-airline-poisson-trend/
  05-train-station-queue-simulation/
  06-baseball-inning-simulation/

run_all.R
```

## How to run

Run one analysis from its folder:

```bash
cd analyses/03-cauchy-location-grid-approximation
Rscript code/cauchy_location_grid_approximation.R
```

Or run all analyses from the project folder:

```bash
Rscript run_all.R
```

Each analysis writes generated plots to its own `figures/` folder and numerical summaries to its own `results/` folder.

## Reproducibility

The code uses base R only. Simulation scripts set random seeds so that results are reproducible. If simulation counts or seeds are changed, numerical summaries may change slightly.

The airline analysis looks for a local file at `analyses/04-airline-poisson-trend/data/planes.txt` with columns `year` and `fatal`. If that file is not present, the script uses an inline table of fatal-accident counts for 1976--1985 so the analysis remains reproducible without redistributing restricted files.

## Disclaimer

This folder contains original code, original explanations, generated figures, and data that may be redistributed. Restricted prompts, solution PDFs, and private data files should not be committed.
