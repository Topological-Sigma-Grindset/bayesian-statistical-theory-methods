# Project 02: Regression, Likelihood, and Mixture Models

This project develops a sequence of Bayesian and likelihood-based statistical analyses in R. The emphasis is on transparent model formulation, reproducible computation, posterior simulation, maximum likelihood estimation, and finite mixture modeling.

The project combines two applied settings:

1. **Radon measurements**: Bayesian linear regression for log radon levels using county and floor indicators.
2. **Baseball home-run rates**: likelihood-based beta modeling and finite mixture models for player-season home-run rates using the Lahman baseball database.

The project is organized into small standalone analyses. Each analysis has its own script, figures, results, and short README.

---

## Statistical methods demonstrated

This project demonstrates:

- Bayesian linear regression with a noninformative prior
- Posterior summaries for regression coefficients and residual scale
- Posterior predictive simulation
- Beta likelihood modeling for continuous rates
- Maximum likelihood estimation by grid search
- Maximum likelihood estimation by numerical optimization with `optim`
- EM algorithms for finite mixture models
- Multi-start optimization for mixture models
- Mixture-based posterior classification probabilities
- Reproducible R scripting with analysis-level outputs

---

## Project structure

```text
02-regression-likelihood-and-mixture-models/
├── README.md
├── run_all.R
├── MANIFEST.txt
└── analyses/
    ├── 01-radon-regression-posterior-prediction/
    ├── 02-beta-model-home-run-rates/
    ├── 03-beta-mle-grid-and-optim/
    ├── 04-normal-halfnormal-mixture-em/
    ├── 05-ortiz-elite-probabilities-halfnormal-model/
    ├── 06-normal-beta-mixture-em/
    └── 07-ortiz-elite-probabilities-beta-mixture/