#' Runs bootstrap using a starting point estimate (parallel implementation)
#' (Algorithm 4 from paper)
#'
#' `bootstrap_parallel()` returns a list containing the beta_star,
#' theta_star, and log_likelihood evaluated at the convergence point
#'
#' @param data_cov n x p matrix of covariates
#' @param data_outcomes n x d matrix of discrete outcomes
#' @param theta_init K x d matrix of binomial initial params
#' @param beta_init (K-1) x (p+1) matrix of beta (weight) initial params
#' @param M_mar number of iterations of the outer EM algorithm
#' @param s_max vector of length d, containing the maximum score for each test
#' @param numBootstraps number of bootstrap samples
#' @param stop_eps_mar epsilon for checking convergence for each outer EM iter
#' @param stop_eps_latent epsilon for checking convergence for each inner EM iter
#' @param num_imps number of multiple imputations
#' @returns A list containing the bootstrap estimates for beta and theta
#' @export
bootstrap_parallel <- function(data_cov, data_outcomes, theta_init,
                               beta_init, M_mar, s_max, numBootstraps=1000,
                               stop_eps_mar=1e-4, stop_eps_latent=1e-4,
                               num_imps=20)
{
  ##### parallel part
  pb = txtProgressBar(max = numBootstraps, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  n = nrow(data_cov)
  p = ncol(data_cov)
  d = ncol(data_outcomes)
  K = dim(theta_init)[1]

  beta_boot = array(NA, c(K-1, p+1, numBootstraps))
  theta_boot = array(NA, c(K, d, numBootstraps))

  output = foreach (jj = 1:numBootstraps, .combine='comb', .packages=c('nnet'),
                    .export=c('fit_mar', 'get_R_matrix', 'fit_latent_var', 'predict_probs', 'logsumexp',
                              'get_loglikelihood', 'changeRef', 'multiple_impute'),
                    .multicombine=TRUE, .options.snow = opts) %dopar%
  {

      ########################## bootstrap

      boot_idx = sample(n, n, replace=TRUE)
      data_cov_boot = data_cov[boot_idx,]
      data_outcomes_boot = data_outcomes[boot_idx,]

      ########################## run analysis

      res_boot = fit_mar(data_cov_boot, data_outcomes_boot,
                         theta_init, beta_init,
                         M_mar=100, s_max,
                         stop_eps_mar=stop_eps_mar,
                         stop_eps_latent=stop_eps_latent,
                         num_imps=20)

      ######################### get results

      as.list(c(c(res_boot$beta_star), c(res_boot$theta_star)))
  }

  #### unpack output data
  index = 0
  for (ii in 1:(p+1))
  {
    for (jj in 1:(K-1))
    {
      index = index + 1
      beta_boot[jj, ii, ] = unlist(output[[index]])
    }
  }
  for (ii in 1:d)
  {
    for (jj in 1:K)
    {
      index = index + 1
      theta_boot[jj, ii, ] = unlist(output[[index]])
    }
  }


  return (list("beta_boot" = beta_boot, "theta_boot" = theta_boot))

}
