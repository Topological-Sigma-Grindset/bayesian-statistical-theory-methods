# Analysis 06: Normal/Beta Mixture EM

This analysis fits a two-component finite mixture to Lahman home-run rates. The elite component is normal; the non-elite component is beta. The beta M-step uses `optim()` on the log scale for alpha and beta.

## Dependencies

```r
install.packages("Lahman")
```

## Run

```r
source("projects/02-regression-likelihood-and-mixture-models/analyses/06-normal-beta-mixture-em/code/06-normal-beta-mixture-em.R")
```
