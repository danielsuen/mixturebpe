# mixturebpe
R implementation of algorithms for fitting mixture of binomial product experts model with missing at random (MAR) data

## Installation 
This package currently only exists on Github.  This package can be installed using the [devtools](https://github.com/hadley/devtools) package.  Run the following code:

```R
install.packages("devtools")
devtools::install_github("danielsuen/mixturebpe") 
```

## Example Code

In the package, we provide working implementations both sequentially and in parallel.  In practice, we recommend running the code in parallel because there are multiple random initializations and many bootstrap samples.  The following example shows how to run code in parallel.

We first load the package and simulate some data.  The data is simulated according to a mixture of binomial product experts model with $d=4$ outcome variables, $p=2$ covariates, and $K=3$ mixtures.  The parameter $\eta$ controls the amount of missingness.

```R
library(doParallel)
library(doSNOW)
library(mixturebpe)

#### init parallel cores

numCores = detectCores()
cluster = makeCluster(numCores-2)
registerDoSNOW(cluster)

#### init bootstrap arrays

numBootstraps = 1000
K = 3
p = 2
d = 4
beta_boot = array(NA, c(K-1, p+1, numBootstraps))
theta_boot = array(NA, c(K, d, numBootstraps))

#### generate simulated data

n = 2000
eta = 2.5

s_max = c(5,10,10,20)

dataset = generate_cc_data(n)
observed = makeMiss(dataset, eta=eta)
X = observed$X
Y = observed$Y
R = observed$R
Y_obs = observed$Y_obs
```

Next, we fit the model using a nested EM approach (Algorithm 3) from the paper:

```R
#### get point estimate

stop_eps_mar = 1e-4
stop_eps_latent = 1e-4

res_mar = fit_mar_randinit_parallel(X, Y_obs, numRandInit=100, K=3,
                                    M_mar=100, s_max,
                                    stop_eps_mar=stop_eps_mar,
                                    stop_eps_latent=stop_eps_latent,
                                    num_imps=20)
theta_mar = res_mar$theta_star
beta_mar = res_mar$beta_star
```

Finally, we run the bootstrap:

```R
##### bootstrap

pb = txtProgressBar(max = numBootstraps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


output = foreach (jj = 1:numBootstraps, .combine='comb', .packages=c('nnet'),
                  .multicombine=TRUE, .options.snow = opts) %dopar%
{

    ########################## bootstrap

    boot_idx = sample(n, n, replace=TRUE)
    X_boot = X[boot_idx,]
    Y_boot = Y_obs[boot_idx,]

    ########################## run analysis

    res_boot = fit_mar(X_boot, Y_boot, theta_mar, beta_mar,
                       K, M_mar=100, s_max,
                       stop_eps_mar=stop_eps_mar,
                       stop_eps_latent=stop_eps_latent,
                       num_imps=20)

    ######################### get results

    as.list(c(c(res_boot$beta_star), c(res_boot$theta_star)))
  }

#### unpack bootstrap results

beta_boot[1, 1, ] = unlist(output[[1]])
beta_boot[2, 1, ] = unlist(output[[2]])
beta_boot[1, 2, ] = unlist(output[[3]])
beta_boot[2, 2, ] = unlist(output[[4]])
beta_boot[1, 3, ] = unlist(output[[5]])
beta_boot[2, 3, ] = unlist(output[[6]])

theta_boot[1, 1, ] = unlist(output[[7]])
theta_boot[2, 1, ] = unlist(output[[8]])
theta_boot[3, 1, ] = unlist(output[[9]])
theta_boot[1, 2, ] = unlist(output[[10]])
theta_boot[2, 2, ] = unlist(output[[11]])
theta_boot[3, 2, ] = unlist(output[[12]])
theta_boot[1, 3, ] = unlist(output[[13]])
theta_boot[2, 3, ] = unlist(output[[14]])
theta_boot[3, 3, ] = unlist(output[[15]])
theta_boot[1, 4, ] = unlist(output[[16]])
theta_boot[2, 4, ] = unlist(output[[17]])
theta_boot[3, 4, ] = unlist(output[[18]])
```

We can construct $95\\%$ confidence intervals for a given parameter (say $\theta_{1,1}$) via the following code:

```R
upper = theta_mar[1,1] + 1.96 * sd(theta_boot[1,1,])
lower = theta_mar[1,1] - 1.96 * sd(theta_boot[1,1,])
```

Thus, the $95\\%$ confidence interval for $\theta_{1,1}$ is given by $(\text{lower}, \text{upper})$.
