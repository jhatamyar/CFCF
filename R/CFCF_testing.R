# testing
#devtools::install_github("jhatamyar/CFCF")


#library(CFCF)

# Generate data.
#n <- 2000
#p <- 10

## matrix of covariates, vector of treatment and outcomes
#X <- matrix(rnorm(n * p), n, p)
#W <- rbinom(n, 1, 0.4 + 0.2 * (X[, 1] > 0))
#Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)


## run the function: NOTE THIS WILL TAKE A FEW MINUTES
#crossfitted_CF <- CFCF(X, Y, W)
