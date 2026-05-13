# Analysis 03: Beta MLE by Grid Search and Optim

This analysis estimates beta parameters for home-run rates using both a grid search and `optim()` on the log-parameter scale. It compares estimates, log-likelihoods, and model-implied variance against empirical variance.

## Dependencies

```r
install.packages("Lahman")
```

## Run

```r
source("projects/02-regression-likelihood-and-mixture-models/analyses/03-beta-mle-grid-and-optim/code/03-beta-mle-grid-and-optim.R")
```
