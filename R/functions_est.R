#' Fit the mixture of binomial product experts model when there is no missingness
#' (Algorithm 1 from paper)
#'
#' `fit_latent_var()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param theta_init K x d matrix of binomial initial params
#' @param beta_init (K-1) x (p+1) matrix of beta (weight) initial params
#' @param K number of mixture components
#' @param M_em number of iterations of the EM algorithm
#' @param s_max vector of length d, containing the maximum score for each test
#' @param printFlag True/False flag to print status
#' @param stop_eps epsilon for checking convergence
#' @returns A list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#' @export
fit_latent_var <- function(data_cov, data_outcomes, theta_init, beta_init,
                           M_em=50, s_max, printFlag=FALSE, stop_eps=1e-4)
{
  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)
  K = dim(theta_init)[1]

  beta_t = beta_init
  theta_t = theta_init

  num_params = prod(dim(beta_t)) + prod(dim(theta_t))

  # get w_k(x)
  logW = predict_probs(beta_t, data_cov, log=TRUE)

  for(j_em in 1:M_em)
  {
    if (printFlag == TRUE)
      print(paste0("EM iteration: ", j_em))

    # get P(Z=k|X,Y)
    logZ = NULL
    # proportional vector
    for (k_iter in 1:K)
    {
      logZ_temp = logW[,k_iter]#rep(logW[k_iter], n)

      for (j_var in 1:d)
      {
        logZ_temp = logZ_temp + dbinom(data_outcomes[, j_var],
                                       s_max[j_var],
                                       theta_t[k_iter, j_var],
                                       log=TRUE)
      }
      logZ = cbind(logZ, logZ_temp)
    }

    ### logsumexp to prevent underflow
    Z = exp(logZ - apply(logZ, 1, logsumexp))

    ## M-step 1: update parameters of covariates
    logistic_fit = multinom(Z ~ .,
                            data=data_cov,
                            trace=FALSE,
                            maxit=1000)
    beta_t = summary(logistic_fit)$coeff
    W = logistic_fit$fitted.value
    logW = log(W)

    # M-step 2: update mixture proportion
    for (k_iter in 1:K)
    {
      for (j_var in 1:d)
      {
        theta_t[k_iter, j_var] = sum(Z[, k_iter] * data_outcomes[, j_var]) /
          (sum(Z[, k_iter]) * s_max[j_var])
      }
    }

    ### check convergence
    if (j_em == 1)
    {
      beta_old = beta_init
      theta_old = theta_init
      diff = sum((beta_t - beta_old)^2) + sum((theta_t - theta_old)^2)
      diff = diff / num_params

      if (printFlag == TRUE)
      {
        print(diff)
        print(get_loglikelihood(data_cov, data_outcomes, s_max, theta_t, beta_t))
      }
      if (diff < stop_eps)
      {
        break
      }

      beta_old = beta_t
      theta_old = theta_t
    }

    if (j_em > 1)
    {
      diff = sum((beta_t - beta_old)^2) + sum((theta_t - theta_old)^2)
      diff = diff / num_params
      beta_old = beta_t
      theta_old = theta_t

      if (printFlag == TRUE)
      {
        print(diff)
        print(get_loglikelihood(data_cov, data_outcomes, s_max, theta_t, beta_t))
      }
      if (diff < stop_eps)
      {
        break
      }
    }
  }

  # get log_likelihood
  ll = get_loglikelihood(data_cov, data_outcomes, s_max, theta_t, beta_t)


  theta_star = theta_t
  beta_star = beta_t

  ### make most healthy group reference category
  change_ref = changeRef(theta_star, beta_star)
  theta_star = change_ref$theta
  beta_star = change_ref$beta

  return (list("beta_star" = beta_star, "theta_star" = theta_star,
               "log_likelihood" = ll))
}

#' Fit the mixture of binomial product experts model when there is missingness
#' under MAR assumption
#' (Algorithm 3 from paper)
#'
#' `fit_mar()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param theta_init K x d matrix of binomial initial params
#' @param beta_init (K-1) x (p+1) matrix of beta (weight) initial params
#' @param K number of mixture components
#' @param M_mar number of iterations of the outer EM algorithm
#' @param s_max vector of length d, containing the maximum score for each test
#' @param printFlag True/False flag to print status
#' @param stop_eps_mar epsilon for checking convergence for each outer EM iter
#' @param stop_eps_latent epsilon for checking convergence for each inner EM iter
#' @param num_imps number of multiple imputations
#' @returns A list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#' @export
fit_mar <- function(data_cov, data_outcomes, theta_init, beta_init,
                    M_mar=50, s_max, printFlag=FALSE, stop_eps_mar=1e-4,
                    stop_eps_latent=1e-4, num_imps=20)
{
  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)
  K = dim(theta_init)[1]

  # get missing pattern matrix
  R = get_R_matrix(data_outcomes, s_max)
  theta_t = theta_init
  beta_t = beta_init

  num_params = prod(dim(beta_t)) + prod(dim(theta_t))

  for (j_mar in 1:M_mar)
  {
    if (printFlag == TRUE)
    {
      print(paste0("Iteration of MAR: ", j_mar))
    }
    ### E-step: impute
    data_outcomes_im = multiple_impute(data_cov, data_outcomes, theta_t, beta_t, s_max, M_imp=num_imps)

    ### M-step

    ### run Algorithm 1 on stacked imputed data
    latent_fit = fit_latent_var(do.call("rbind", replicate(num_imps, data_cov, simplify = FALSE)),
                                data_outcomes_im, theta_t, beta_t,
                                M_em=50, s_max, printFlag=printFlag, stop_eps=stop_eps_latent)

    beta_t = latent_fit$beta_star
    theta_t = latent_fit$theta_star

    ### check convergence
    if (j_mar == 1)
    {
      beta_mar_old = beta_init
      theta_mar_old = theta_init

      diff_mar = sum((beta_t - beta_mar_old)^2) + sum((theta_t - theta_mar_old)^2)
      diff_mar = diff_mar / num_params
      if (diff_mar < stop_eps_mar)
      {
        break
      }

      beta_mar_old = beta_t
      theta_mar_old = theta_t
    }

    if (j_mar > 1)
    {
      diff_mar = sum((beta_t - beta_mar_old)^2) + sum((theta_t - theta_mar_old)^2)

      diff_mar = diff_mar / num_params
      beta_mar_old = beta_t
      theta_mar_old = theta_t

      #print(get_loglikelihood(data_cov, data_outcomes, s_max, theta_t, beta_t))
      #print(diff_mar)
      if (diff_mar < stop_eps_mar)
      {
        break
      }
    }
  }

  # get log-likelihood
  ll = get_loglikelihood(data_cov, data_outcomes, s_max, theta_t, beta_t)

  ## change reference class
  change = changeRef(theta_t, beta_t)
  theta_star = change$theta
  beta_star = change$beta

  return (list("beta_star" = beta_t, "theta_star" = theta_t,
               "log_likelihood"=ll))
}


#' Fit the mixture of binomial product experts model when there is missingness
#' under MAR assumption with random initializations (serial implementation)
#' (Algorithm 3 from paper)
#'
#' `fit_mar_randinit()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param numRandInit number of random initializations
#' @param K number of mixture components
#' @param M_mar number of iterations of the outer EM algorithm
#' @param s_max vector of length d, containing the maximum score for each test
#' @param printFlag True/False flag to print status
#' @param stop_eps_mar epsilon for checking convergence for each outer EM iter
#' @param stop_eps_latent epsilon for checking convergence for each inner EM iter
#' @param num_imps number of multiple imputations
#' @returns A list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#' @export
fit_mar_randinit <- function(data_cov, data_outcomes, numRandInit,
                             K, M_mar, s_max, printFlag=FALSE, stop_eps_mar=1e-4,
                             stop_eps_latent=1e-4, num_imps=20)
{
  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)

  ll_vect = rep(NA, numRandInit)
  theta_all_iters = array(NA, c(K, d, numRandInit))
  beta_all_iters = array(NA, c(K-1, p+1, numRandInit))

  ### run fit_mar for every random init
  for (j_randinit in 1:numRandInit)
  {
    if (printFlag == TRUE)
      print(paste0("Current random init: ", j_randinit, " ---------"))

    ### random initializations
    beta_init = cbind(runif(K-1, -1, 1),
                      matrix(runif((K-1)*(p), -0.5, 0.5), K-1, p))
    theta_init = matrix(runif(K*d, min=0.1, max=0.9), K, d)
    theta_init = theta_init[order(theta_init[, 1]), ]

    beta_t = beta_init
    theta_t = theta_init

    res_mar = fit_mar(data_cov, data_outcomes, theta_t, beta_t,
                      M_mar, s_max, printFlag, stop_eps_mar=stop_eps_mar,
                      stop_eps_latent=stop_eps_latent,
                      num_imps=num_imps)

    theta_t = res_mar$theta_star
    beta_t = res_mar$beta_star
    ll = res_mar$log_likelihood

    ll_vect[j_randinit] = ll
    theta_all_iters[,,j_randinit] = theta_t
    beta_all_iters[,,j_randinit] = beta_t


    max_ll_idx = which.max(ll_vect)
  }

  ### pick beta, theta that maximizes log_likelihood
  beta_star = beta_all_iters[,,max_ll_idx]
  theta_star = theta_all_iters[,,max_ll_idx]

  ll = get_loglikelihood(data_cov, data_outcomes, s_max, theta_star, beta_star)

  ### change reference class
  change = changeRef(theta_star, beta_star)
  theta_star = change$theta
  beta_star = change$beta

  #### save all results
  # for (j_randinit in 1:numRandInit)
  # {
  #   change = changeRef(theta_all_iters[,,j_randinit], beta_all_iters[,,j_randinit])
  #   theta_all_iters[,,j_randinit] = change$theta
  #   beta_all_iters[,,j_randinit] = change$beta
  # }
  #
  # save(ll_vect, theta_all_iters, beta_all_iters,
  #      file=paste("all_results_K", K, ".Rdata", sep=""))

  return (list("beta_star" = beta_star, "theta_star" = theta_star,
               "log_likelihood"=ll_vect))
}

#' Fit the mixture of binomial product experts model when there is missingness
#' under MAR assumption with random initializations (parallel implementation)
#' (Algorithm 3 from paper)
#'
#' `fit_mar_randinit_parallel()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param numRandInit number of random initializations
#' @param K number of mixture components
#' @param M_mar number of iterations of the outer EM algorithm
#' @param s_max vector of length d, containing the maximum score for each test
#' @param printFlag True/False flag to print status
#' @param stop_eps_mar epsilon for checking convergence for each outer EM iter
#' @param stop_eps_latent epsilon for checking convergence for each inner EM iter
#' @param num_imps number of multiple imputations
#' @returns A list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#' @export
fit_mar_randinit_parallel <- function(data_cov, data_outcomes, numRandInit,
                                      K, M_mar, s_max, printFlag=FALSE,
                                      stop_eps_mar=1e-4, stop_eps_latent=1e-4,
                                      num_imps=20)
{
  ##### parallel part
  pb = txtProgressBar(max = numRandInit, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)


  ll_vect = rep(NA, numRandInit)
  theta_all_iters = array(NA, c(K, d, numRandInit))
  beta_all_iters = array(NA, c(K-1, p+1, numRandInit))
  beta_init_all_iters = array(NA, c(K-1, p+1, numRandInit))
  theta_init_all_iters = array(NA, c(K, d, numRandInit))

  #### run analysis in parallel
  output = foreach (jj = 1:numRandInit, .combine='comb', .packages=c('nnet'),
                    .export=c('fit_mar', 'get_R_matrix', 'fit_latent_var', 'predict_probs', 'logsumexp',
                              'get_loglikelihood', 'changeRef', 'multiple_impute'),
                    .multicombine=TRUE, .options.snow = opts) %dopar%
    {

      ######################### random initializations

      beta_init = cbind(runif(K-1, -1, 1),
                        matrix(runif((K-1)*(p), -0.5, 0.5), K-1, p))
      theta_init = matrix(runif(K*d, min=0.1, max=0.9), K, d)
      theta_init = theta_init[order(theta_init[, 1]), ]

      beta_t = beta_init
      theta_t = theta_init

      ########################## run analysis

      res_mar = fit_mar(data_cov, data_outcomes, theta_t, beta_t,
                        M_mar, s_max, printFlag,
                        stop_eps_mar=stop_eps_mar,
                        stop_eps_latent=stop_eps_latent,
                        num_imps=num_imps)

      ######################### get results

      as.list(c(c(res_mar$beta_star), c(res_mar$theta_star), c(beta_init), c(theta_init), res_mar$log_likelihood))
    }

  max_ll_idx = which.max(unlist(output[[length(output)]]))

  #### unpack output data
  index = 0
  for (ii in 1:(p+1))
  {
    for (jj in 1:(K-1))
    {
      index = index + 1
      beta_all_iters[jj, ii, ] = unlist(output[[index]])
    }
  }
  for (ii in 1:d)
  {
    for (jj in 1:K)
    {
      index = index + 1
      theta_all_iters[jj, ii, ] = unlist(output[[index]])
    }
  }
  for (ii in 1:(p+1))
  {
    for (jj in 1:(K-1))
    {
      index = index + 1
      beta_init_all_iters[jj, ii, ] = unlist(output[[index]])
    }
  }
  for (ii in 1:d)
  {
    for (jj in 1:K)
    {
      index = index + 1
      theta_init_all_iters[jj, ii, ] = unlist(output[[index]])
    }
  }

  index = index + 1
  ll_vect = unlist(output[[index]])

  ### pick beta, theta that maximizes log_likelihood
  beta_star = beta_all_iters[,,max_ll_idx]
  theta_star = theta_all_iters[,,max_ll_idx]

  #### best initialization point
  beta_best_init = beta_init_all_iters[,,max_ll_idx]
  theta_best_init = theta_init_all_iters[,,max_ll_idx]

  ll = ll_vect[max_ll_idx]

  # change reference class
  change = changeRef(theta_star, beta_star)
  theta_star = change$theta
  beta_star = change$beta

  ### save results
#   for (j_randinit in 1:numRandInit)
#   {
#     change = changeRef(theta_all_iters[,,j_randinit], beta_all_iters[,,j_randinit])
#     theta_all_iters[,,j_randinit] = change$theta
#     beta_all_iters[,,j_randinit] = change$beta
#   }
#
#   save(ll_vect, theta_all_iters, beta_all_iters,
#        file=paste("all_results_K", K, ".Rdata", sep=""))



  return (list("beta_star" = beta_star, "theta_star" = theta_star,
               "log_likelihood"=ll_vect))

}

#' Multiple imputation
#' (Algorithm 2 from paper)
#'
#' `multiple_impute()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param theta K x d matrix of binomial theta params
#' @param beta (K-1) x (p+1) matrix of beta (weight) initial params
#' @param M_imps number of multiple imputations
#' @param stacked True/False flag on whether to stack the imputed data sets
#' @returns All M imputed data sets
#' @export
multiple_impute = function(data_cov, data_outcomes,
                           theta, beta, s_max,
                           M_imp=50, stacked=TRUE)
{
  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)
  K = nrow(beta) + 1

  # get missing pattern matrix
  R = get_R_matrix(data_outcomes, s_max)

  if (stacked == TRUE)
  {
    data_outcomes_imputed = array(NA, dim=c(n*M_imp, d))
  } else
  {
    data_outcomes_imputed = array(NA, dim=c(n,d,M_imp))
  }

  num_params = dim(beta)[1] * dim(beta)[2]
  num_params = num_params + dim(theta)[1] * dim(theta)[2]

  ### impute data M times
  for (j_imp in 1:M_imp)
  {
    data_outcomes_temp = data_outcomes
    R_im = R

    ### impute data for each pattern
    curr_miss = (R_im %*% rep(1, d)) != d
    while (sum(curr_miss) > 0)
    {
      ## calculate weights
      pattern_idx = which(curr_miss)[1]

      idx_im = (colSums(abs(t(R_im) - R_im[pattern_idx,])) == 0)

      logW_im = predict_probs(beta, data_cov[idx_im,], log=TRUE)

      for (k_iter in 1:K)
      {
        for (d_im in which(!!R_im[pattern_idx,]))
        {
          logW_im[, k_iter] = logW_im[, k_iter] +
            dbinom(data_outcomes_temp[idx_im, d_im], size=s_max[d_im],
                   prob = theta[k_iter, d_im], log=TRUE)
        }
      }

      # logsumexp to prevent underflow
      W_im = exp(logW_im - apply(logW_im, 1, logsumexp))



      for(d_im in which(!R_im[pattern_idx,]))
      {

        #### sample Z
        class_sample = rep(NA, dim(W_im)[1])
        for (ii in 1:dim(W_im)[1])
        {
          class_sample[ii] = sample(K, 1, prob=W_im[ii,])
        }

        #### sample Y_miss
        data_outcomes_temp[idx_im,d_im] = rbinom(n=length(class_sample),
                                                 size=s_max[d_im],
                                                 prob=theta[class_sample,d_im])

      }

      R_im[idx_im,] = 1
      # this denotes that we have completed the imputation for this pattern
      curr_miss = (R_im %*% rep(1,d) != d)

    }

    if (stacked == TRUE)
    {
      data_outcomes_imputed[((j_imp-1)*n+1):(j_imp*n),] = as.matrix(data_outcomes_temp)
    } else
    {
      data_outcomes_imputed[,,j_imp] = as.matrix(data_outcomes_temp)
    }

  }

  return (data_outcomes_imputed)
}
