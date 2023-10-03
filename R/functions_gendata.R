#' Generates a complete case data set without missing data
#'
#' `generate_cc_data()` returns a generated data set
#'
#' @param n number of observations
#' @returns A list of X (n x p) covariate matrix, Y (n x d) discrete outcome
#' matrix, Z (n x 1) vector of class labels
#' @export
generate_cc_data <- function(n)
{
  n1 = 5
  n2 = 10
  n3 = 10
  n4 = 20

  rho = 0.2
  d = 4
  X = rmvnorm(n, mean=c(2,3), sigma = (rho*(1-diag(2)) + diag(2)))

  beta2 = c(-1.5, 0.3, 0.4)
  beta3 = c(-2, 0.5, 0.25)

  beta_matrix = rbind(beta2, beta3)
  probs = predict_probs(beta_matrix, X)

  theta = matrix(c(0.8, 0.5, 0.2), 3, 1)
  theta = matrix(rep(theta, d), nrow=length(theta))


  Y = matrix(0, nrow=n, ncol=4)


  class = rep(NA, n)

  for (ii in 1:n)
  {
    class[ii] = sample(c(1,2,3), size=1, replace=TRUE,
                       prob=c(probs[ii,1], probs[ii,2], probs[ii,3]))

    if (class[ii] == 1)
    {
      Y[ii, 1] = rbinom(n=1, size=n1, prob=theta[1,1])
      Y[ii, 2] = rbinom(n=1, size=n2, prob=theta[1,2])
      Y[ii, 3] = rbinom(n=1, size=n3, prob=theta[1,3])
      Y[ii, 4] = rbinom(n=1, size=n4, prob=theta[1,4])
    } else if (class[ii] == 2)
    {
      Y[ii, 1] = rbinom(n=1, size=n1, prob=theta[2,1])
      Y[ii, 2] = rbinom(n=1, size=n2, prob=theta[2,2])
      Y[ii, 3] = rbinom(n=1, size=n3, prob=theta[2,3])
      Y[ii, 4] = rbinom(n=1, size=n4, prob=theta[2,4])
    } else if (class[ii] == 3)
    {
      Y[ii, 1] = rbinom(n=1, size=n1, prob=theta[3,1])
      Y[ii, 2] = rbinom(n=1, size=n2, prob=theta[3,2])
      Y[ii, 3] = rbinom(n=1, size=n3, prob=theta[3,3])
      Y[ii, 4] = rbinom(n=1, size=n4, prob=theta[3,4])
    }
  }

  X = data.frame(X)
  Y = data.frame(Y)

  return (list("X"=X, "Y"=Y, "Z"=class))
}

#' Computes expit of a vector
#'
#' `expit()` returns expit of a vector
#'
#' @param input numerical vector
#' @returns A vector containing expit(input)
#' @export
expit = function(input)
{
  return (exp(input) / (1+exp(input)))
}

#' Adds missingness to a generated complete case data set
#'
#' `makeMiss()` returns a generated data set with missing data
#'
#' @param dataset list of generated data set containing X, Y, Z
#' @param eta parameter that controls missingness, eta -> inf => less missingness,
#' eta -> 0 => most missingness
#' @returns A list of X (n x p) covariate matrix, Y (n x d) discrete outcome
#' matrix, Y_obs (n x d) discrete outcome matrix with NaNs,
#' Z (n x 1) vector of class labels
#' @export
makeMiss = function(dataset, eta=0)
{

  p = 2
  d = 4

  X = dataset$X
  Y = dataset$Y
  class = dataset$Z

  n1 = 5
  n2 = 10
  n3 = 10
  n4 = 20

  for (r in c(0001, 0110, 1010, 1110))
  {
    if (r == 0001)
    {
      lower_bound = 0
      upper_bound = 1

      intercept = -2 - eta
      beta_X = c(-0.25, 0.3)
      beta_Y = c(1/n4*3)

      beta = c(intercept, beta_X, beta_Y)

      input = t(t(cbind(rep(1, n), X, Y[,4]))) %*% as.matrix(beta)

      prob_R0001 = expit(input)
      prob_R0001 = pmin(pmax(prob_R0001, lower_bound), upper_bound)
    } else if (r == 0110)
    {
      lower_bound = 0
      upper_bound = 1

      intercept = -1 - eta
      beta_X = c(0.3, -0.7)
      beta_Y = c(1/n2*-1, 1/n3*1.5)

      beta = c(intercept, beta_X, beta_Y)

      input = t(t(cbind(rep(1, n), X, Y[,c(2,3)]))) %*% as.matrix(beta)

      prob_R0110 = expit(input)
      prob_R0110 = pmin(pmax(prob_R0110, lower_bound), upper_bound)
    } else if (r == 1010)
    {
      lower_bound = 0
      upper_bound = 1

      intercept = -2 - eta
      beta_X = c(0.7, -0.4)
      beta_Y = c(1/n1*1.2, 1/n3*-1.5)

      beta = c(intercept, beta_X, beta_Y)

      input = t(t(cbind(rep(1, n), X, Y[,c(1,3)]))) %*% as.matrix(beta)

      prob_R1010 = expit(input)
      prob_R1010 = pmin(pmax(prob_R1010, lower_bound), upper_bound)
    } else if (r == 1110)
    {
      lower_bound = 0
      upper_bound = 1

      intercept = -1 - eta
      beta_X = c(0.2, -0.15)
      beta_Y = c(1/n1*0.75, 1/n2*-1.4, 1/n3*-0.5)

      beta = c(intercept, beta_X, beta_Y)

      input = t(t(cbind(rep(1, n), X, Y[,c(1,2,3)]))) %*% as.matrix(beta)

      prob_R1110 = expit(input)
      prob_R1110 = pmin(pmax(prob_R1110, lower_bound), upper_bound)
    }
  }

  prob_R1111 = 1 - prob_R0001 - prob_R0110 - prob_R1010 - prob_R1110

  R = rep(NA, n)
  for (ii in 1:n)
  {
    R[ii] = sample(c(0001, 0110, 1010, 1110, 1111), 1,
                   prob=c(prob_R0001[ii], prob_R0110[ii], prob_R1010[ii], prob_R1110[ii],
                          prob_R1111[ii]))
  }

  Y_obs = Y

  Y_obs[R==0001, c(1,2,3)] = NA
  Y_obs[R==0110, c(1,4)] = NA
  Y_obs[R==1010, c(2,4)] = NA
  Y_obs[R==1110, 4] = NA

  return (list("X"=X, "Y"=Y, "Y_obs"=Y_obs, "Z"=class, "R"=R))
}
