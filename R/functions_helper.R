library(nnet)
library(mvtnorm)
library(doSNOW)
library(doParallel)

### helper functions


comb = function(x, ...) {
  lapply(seq_along(x),
         function(ii) c(x[[ii]], lapply(list(...), function(y) y[[ii]])))
}

#' Logsumexp trick
#'
#' `logsumexp()` computes the softmax of a vector, using the logsumexp trick
#'
#' @param log_vect Input numerical vector.
#' @returns A probability vector, each element in [0,1] and all sum to 1
#' @export
logsumexp = function(log_vect)
{
  c = max(log_vect)
  return (as.vector(c + log(sum(exp(log_vect - c)))))
}

#' Predict class probabilities w_k(x) for each obs
#'
#' `predict_probs()` retrieves w_k(x) for each obs
#'
#' @param beta_matrix (K-1) x (p+1) matrix of beta (weight) parameters
#' @param cov_matrix n x p matrix of covariates for all obs
#' @param log True/False flag to determine whether to return probs in log form
#' @returns An n x K matrix of (log) probabilities for each class and obs
#' @export
predict_probs <- function(beta_matrix, cov_matrix, log=FALSE)
{
  n = nrow(cov_matrix)

  prod = cbind(rep(0,n),
               t(beta_matrix %*% t(cbind(rep(1, n),
                                         cov_matrix))))

  logW = prod - apply(prod, 1, logsumexp)

  if (log == FALSE)
  {
    W = exp(logW)

    return (W)
  } else {
    return (logW)
  }
}

#' Obtain the missing pattern matrix R
#'
#' `get_R_matrix()` returns the n x d missing pattern matrix R
#'
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param s_max vector of length d, containing the maximum score for each test
#' @returns An n x d binary matrix (1 indicates observed, 0 is missing)
#' @export
get_R_matrix <- function(data_outcomes, s_max)
{
  R = NULL
  n = nrow(data_outcomes)
  d = ncol(data_outcomes)

  for (j_var in 1:d)
  {
    idx = (data_outcomes[, j_var] <= s_max[j_var])
    idx = idx & (data_outcomes[, j_var] >= 0)
    R = cbind(R, idx)
  }

  R[is.na(R)] = 0

  R = R * 1

  return (R)
}

#' Shifts the theta and beta parameters to account for changing class reference
#'
#' `changeRef()` returns new theta and beta parameters, shifted, so that the
#' new reference class is such that `theta[,1]' is maximized
#'
#' @param theta K x d matrix of binomial parameters
#' @param beta (K-1) x (p+1) matrix of beta (weight) parameters
#' @returns A list containing the updated theta and beta matrices
#' @export
changeRef <- function(theta, beta)
{

  if (sum(is.na(theta)) + sum(is.na(beta)) != 0)
  {
    return (list("theta" = theta, "beta" = beta))
  }

  K = nrow(theta)
  p = ncol(beta) - 1
  d = ncol(theta)

  curr_ref = rank(theta[,1])[1]
  desired_ref = which.max(theta[,1])

  if (curr_ref != dim(theta)[1])
  {
    if (K > 2)
    {
      beta_temp = matrix(NA, K-1, p+1)
      for (ii in 1:K-1)
      {
        beta_temp[ii, ] = beta[ii, ] - beta[desired_ref-1, ]
      }
      beta_temp[desired_ref-1, ] = -1 * beta[desired_ref-1, ]

      idx = order(theta[,1], decreasing=TRUE)
      theta = theta[idx, ]
      beta = beta_temp
    } else
    {
      beta = -beta
      idx = order(theta[,1], decreasing=TRUE)
      theta = theta[idx, ]
    }
  } else
  {
    idx = order(theta[,1], decreasing=TRUE)
    theta = theta[idx, ]
    if (K > 2)
    {
      beta = beta[idx[-1]-1, ]
    }
  }

  return (list("theta" = theta, "beta" = beta))
}

#' Computes the observed log-likelihood
#'
#' `get_loglikelihood()` computes the observed log-likelihood at (theta, beta)
#' using the data
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param s_max vector of length d, containing the maximum score for each test
#' @param theta K x d matrix of binomial parameters
#' @param beta (K-1) x (p+1) matrix of beta (weight) parameters
#' @returns The observed log-likelihood value
#' @export
get_loglikelihood <- function(data_cov, data_outcomes, s_max, theta, beta)
{
  n = nrow(data_cov)
  K = nrow(theta)
  d = ncol(theta)

  W = predict_probs(beta, data_cov)
  if (sum(is.na(W)) != 0)
  {
    return (NA)
  }
  R = get_R_matrix(data_outcomes, s_max)
  cc_idx = ((R  %*% rep(1, d)) == d)


  for (k_iter in 1:K)
  {
    for (d_iter in 1:d)
    {
      W[cc_idx, k_iter] = W[cc_idx, k_iter] *
        dbinom(data_outcomes[cc_idx, d_iter], size=s_max[d_iter],
               prob = theta[k_iter, d_iter])
    }
  }

  curr_miss = !cc_idx


  while (sum(curr_miss) > 0)
  {
    pattern_idx = which(curr_miss)[1]

    idx_im = (colSums(abs(t(R) - R[pattern_idx,])) == 0)

    for (k_iter in 1:K)
    {
      for (d_im in which(!!R[pattern_idx,]))
      {
        W[idx_im, k_iter] = W[idx_im, k_iter] *
          dbinom(data_outcomes[idx_im, d_im], size=s_max[d_im],
                 prob = theta[k_iter, d_im])
      }
    }

    R[idx_im,] = 1
    # this denotes that we have completed the imputation for this pattern
    curr_miss = (R %*% rep(1,d) != d)
  }
  return (sum(log(rowSums(W))))
}
