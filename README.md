# CFCF
Cross-fitted Causal Forests 

This package implements a causal forest (Wager and Athey, 2018) learned (cross-fitted) with training and holdout data, fit over K folds, as in Hatamyar and Kreif (2023). 
Forests are estimated with all K-1 folds, then estimates are predicted for the holdout (kth) fold. This is done K times with each fold used as the holdout fold once.
This is repeated a specified number of times and the results are averaged across folds and repetitions. See Hatamyar and Kreif (2023) for further details. 


## Installation


```R
# Install the development version from GitHub:
devtools::install_github("jhatamyar/CFCF")
```

## Usage 

```R
library(CFCF)

# Generate data.
n <- 2000
p <- 10 

## matrix of covariates, vector of treatment and outcomes 
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)


## run the function: NOTE THIS WILL TAKE A FEW MINUTES 
crossfitted_CF <- CFCF(X, Y, W)
```

## Outputs 

```R
## CATES (Conditional Average Treatment Effects)
crossfitted_CF$cates

## ATE (Average Treatment Effect) 
crossfitted_CF$ATE

## Double-Robust Scores
crossfitted_CF$scores
```

The function also outputs standard errors of CATEs and ATE, estimated propensity scores, and the separate double-robust scores for treatment and control. 

### References
Julia Hatamyar and Noemi Kreif. **Policy Learning with Rare Outcomes**, 2023, [arxiv](https://arxiv.org/abs/2302.05260)

Stefan Wager and Susan Athey. **Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.** *Journal of the American Statistical Association*, 113(523), 2018.



