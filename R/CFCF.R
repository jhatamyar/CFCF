#' Causal forest learned (cross-fitted) with training and holdout data, fit
#' over K folds. Forests are estimated with all but i folds,
#' then estimates are predicted for the holdout (ith) fold. This is repeated
#' a specified number of times and the results are averaged across folds.
#' @importFrom stats predict quantile runif
#' @importFrom foreach foreach %do%
#' @importFrom grf causal_forest regression_forest
#'
#' @param X a matrix of covariates, dimension n rows by m columns
#' @param Y a vector of outcomes, dimension 1xn
#' @param W a vector of treatment status, dimension 1Xn
#' @param K the number of folds to be used in cross-fitting, default is 5
#' @param reps the number of times to repeat the cross-fitting, default is 4
#' @param num.trees the number of trees to be fitted in the causal forest, default is 2000
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item{\code{cates}}{A numeric vector containing the conditional average treatment effects for each unit or group analyzed, averaged over repetitions.}
#'     \item{\code{errors}}{A numeric vector of errors for each CATE estimation.}
#'     \item{\code{scores}}{A numeric vector of double-robust scores calculated within the function, according to "get_scores" in the grf package.}
#'     \item{\code{pscores}}{A numeric vector with the propensity scores estimating the probability of assignment to treatment given covariates.}
#'     \item{\code{gamma0}}{A numeric value or vector representing the double-robust score for the control group, according to "double_robust_scores" in the grf package.}
#'     \item{\code{gamma1}}{A numeric value or vector representing the double-robust score for the treated group, according to "double_robust_scores" in the grf package}
#'     \item{\code{ATE}}{A numeric value indicating the Average Treatment Effect estimated by the function, averaged over repetitions.}
#'     \item{\code{SE}}{A numeric vector containing the standard errors associated with ethe overall ATE.}
#'   }
#'
#' Note that the 'alpha' and 'honesty.fraction' parameters within the causal forest function
#' will be tuned within the training step, see documentation for grf package


CFCF <- function(X, Y, W,
                              K = 5,
                              reps = 4,
                              num.trees = 2000){

  r <- foreach(t = 1:reps, .combine='cbind') %do% {

    print(paste("Repetition", toString(t)))

    ## initialize results for this repetition
    res      <- array(NA, dim = c(length(Y),1,6))
    dimnames(res)[[3]] <- c("cates", "errors", "scores", "pscores", "gamma0", "gamma1")

    cates    <- array(NA, dim = length(Y))
    errors   <- array(NA, dim = length(Y))
    scores   <- array(NA, dim = length(Y))
    pscores   <- array(NA, dim = length(Y))
    gamma0   <- array(NA, dim = length(Y))
    gamma1   <- array(NA, dim = length(Y))

    ## create fold indicators
    split             <- runif(nrow(X))
    cvgroup           <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))

    ## loop over folds
    for(i in 1:K){

      ## train using training set !i
      cf.train <- causal_forest(X[cvgroup != i,], Y[cvgroup != i], W[cvgroup != i],
                                num.trees = num.trees, # need more for accurate CIs ?
                                tune.parameters = c("alpha", "honesty.fraction"))

      # predict using holdout fold i
      print(paste("CF fold", toString(i)))
      test_predictions <- predict(cf.train, newdata = X[cvgroup == i,], estimate.variance = TRUE)

      # We need estimates of Y and W for the score : OK to do this with train/test set?
      Y.forest.train <- regression_forest(X[cvgroup != i,], Y[cvgroup != i]) ## change this later to include other args

      W.forest.train <- regression_forest(X[cvgroup != i,], W[cvgroup != i])

      ## get nuisances from test set
      Y.hat.test    <- predict(Y.forest.train, X[cvgroup == i,])$predictions
      W.hat.test    <- predict(W.forest.train, X[cvgroup == i,])$predictions


      # extract cates and variance
      tau.hat <- test_predictions$predictions
      se.est  <- sqrt(test_predictions$variance.estimates)

      ## save for this OOS fold index
      cates[cvgroup == i]  <- tau.hat
      errors[cvgroup == i] <- se.est
      pscores[cvgroup == i] <- W.hat.test

      ## get the scores (this is why I estimated W.hat and Y.hat again)
      ## NOTE this is adapted from "get_scores" grf source code
      debiasing.weights    <- (W[cvgroup == i] - W.hat.test) / (W.hat.test * (1 - W.hat.test))
      Y.residual           <- Y[cvgroup == i] - (Y.hat.test + tau.hat * (W[cvgroup == i] - W.hat.test))
      scores.raw            <- tau.hat + debiasing.weights * Y.residual

      scores[cvgroup == i] <- scores.raw #

      ### now I need to get the gamma scores for policytree
      # conditional means
      Y.hat.0 <- Y.hat.test - W.hat.test * tau.hat
      Y.hat.1 <- Y.hat.test + (1 - W.hat.test) * tau.hat

      mu.matrix <- cbind("control" = Y.hat.0, "treated" = Y.hat.1)

      ## this is taken from double_robust_scores
      W.hat.matrix <- cbind(1 - W.hat.test, W.hat.test) # [control, treated]
      n.obs <- nrow(W.hat.matrix)
      observed.treatment.idx <- cbind(1:n.obs, W[cvgroup == i] + 1)

      YY <- matrix(0, n.obs, 2)
      IPW <- matrix(0, n.obs, 2)
      YY[observed.treatment.idx] <- Y[cvgroup == i]
      IPW[observed.treatment.idx] <- 1 / W.hat.matrix[observed.treatment.idx]
      Gamma.matrix <- (YY - mu.matrix) * IPW + mu.matrix

      gamma0[cvgroup == i] <- Gamma.matrix[,1]
      gamma1[cvgroup == i] <- Gamma.matrix[,2]

    }

    #get ATE and SE for this repetition
    ATE <- mean(scores)
    SE  <- sqrt(mean((scores - ATE)^2) / length(W))

    # save to output arrays
    res[,,1] <- cates
    res[,,2] <- errors
    res[,,3] <- scores
    res[,,4] <- pscores
    res[,,5] <- gamma0
    res[,,6] <- gamma1

    res <- c(as.vector(res), as.vector(ATE), as.vector(SE))

    r <- data.frame(res) #
  }


  ### post estimation (ATE, etc)

  ## extract from large matrix
  cates_all   <- r[1:nrow(X),]
  errors_all  <- r[(nrow(X)+1):(2*nrow(X)),]
  scores_all  <- r[((2*nrow(X))+1):(3*nrow(X)),]
  pscores_all <- r[((3*nrow(X))+1):(4*nrow(X)),]
  gamma0_all  <- r[((4*nrow(X))+1):(5*nrow(X)),]
  gamma1_all  <- r[((5*nrow(X))+1):(6*nrow(X)),]

  ATE_all    <- r[(length(r[,1])-1),]
  SE_all     <- r[(length(r[,1])),]

  ## get column means (average over the repetitions)
  res_cates     <- rowMeans(cates_all)
  res_errors    <- rowMeans(errors_all)
  res_scores    <- rowMeans(scores_all)
  res_pscores   <- rowMeans(pscores_all)
  res_gamma0    <- rowMeans(gamma0_all)
  res_gamma1    <- rowMeans(gamma1_all)

  res_ATE      <- rowMeans(ATE_all)  # -0.003
  res_SE       <- rowMeans(SE_all)

  res_all      <- list(res_cates, res_errors, res_scores, res_pscores, res_gamma0, res_gamma1, res_ATE, res_SE)
  names(res_all) <- c("cates", "errors", "scores", "pscores", "gamma0", "gamma1", "ATE", "SE")

  return(res_all)


}

