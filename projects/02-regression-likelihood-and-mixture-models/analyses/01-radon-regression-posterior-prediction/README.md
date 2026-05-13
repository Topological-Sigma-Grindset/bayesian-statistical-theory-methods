# Analysis 01: Radon Regression and Posterior Prediction

This analysis fits a Bayesian linear regression model for log radon using county and floor indicators. It reports posterior summaries for regression coefficients, the residual scale, and posterior predictive intervals for a new Clay County house on the first floor and in the basement.

## Data

Place the local radon CSV in this folder:

```text
data/
```

The script looks for common file names including `radon.table.7.3.csv`, `radon_table_7_3.csv`, `radon.csv`, and `table_7_3_radon.csv`. If exactly one CSV is present in `data/`, it uses that file.

Expected variables are radon measurement, county, and floor. If no local radon CSV is present, the script writes `results/radon_data_requirement.txt` and skips the calculation.

## Dependencies

Base R only.

## Run

```r
source("projects/02-regression-likelihood-and-mixture-models/analyses/01-radon-regression-posterior-prediction/code/01-radon-regression-posterior-prediction.R")
```
