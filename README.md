# Bayesian Statistical Theory and Methods

This repository is a chronological portfolio of Bayesian statistical modeling projects.

The projects are in probability modeling, posterior computation, uncertainty quantification, simulation, hierarchical modeling, model checking, and reproducible analysis in R.

## Portfolio goals

This repository is designed to show:

- Clear statistical reasoning from model assumptions to posterior inference
- Reproducible Bayesian computation in R
- Clean project organization and documentation
- Interpretation of posterior estimates, intervals, and predictive quantities
- Model checking and sensitivity analysis where appropriate

## Bayesian workflow

Each project is organized around the core Bayesian workflow:

-> Specify a probability model for observed and unobserved quantities 
-> Condition on observed data to compute and summarize the posterior distribution.
-> Evaluate the fit, implications, and sensitivity of the model.

## Projects

| Project | Topic | Main methods | Status |
|---|---|---|---|
| 01 | Probability, Bayes' rule, and basic simulation | Conditional probability, prior, likelihood, posterior | In progress |
| 02 | Single-parameter Bayesian models | Beta-binomial, normal-normal, gamma-Poisson | Planned |
| 03 | Posterior simulation | Direct simulation, grid approximation, posterior summaries | Planned |
| 04 | Multiparameter models | Normal model, nuisance parameters, logistic bioassay | Planned |
| 05 | Hierarchical models | Partial pooling, exchangeability, hierarchical normal models | Planned |
| 06 | Model checking | Posterior predictive checks, sensitivity analysis | Planned |
| 07 | Bayesian computation | MCMC, Gibbs sampling, Metropolis-Hastings | Planned |

## Repository structure

```text
projects/
  01-probability-and-bayes-rule/
  02-single-parameter-models/
  03-posterior-simulation/
  04-multiparameter-models/
  05-hierarchical-models/
  06-model-checking/
  07-bayesian-computation/

shared/
  functions/
  plotting/

docs/
  index.md

references/
  README.md
  README.md
