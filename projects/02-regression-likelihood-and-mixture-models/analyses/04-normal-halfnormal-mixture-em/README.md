# Analysis 04: Normal/Half-Normal Mixture EM

This analysis fits a two-component finite mixture to Lahman home-run rates. The elite component is normal; the non-elite component is half-normal with fixed mode at zero. The script runs multiple EM starts, keeps the best observed log-likelihood, and saves fitted density curves.

## Dependencies

```r
install.packages("Lahman")
```

## Run

```r
source("projects/02-regression-likelihood-and-mixture-models/analyses/04-normal-halfnormal-mixture-em/code/04-normal-halfnormal-mixture-em.R")
```
