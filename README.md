# mixturebpe
R implementation of algorithms for fitting mixture of binomial product experts model with missing at random (MAR) data

## Installation 
This package currently only exists on Github.  This package can be installed using the [devtools](https://github.com/hadley/devtools) package.  Run the following code:

```R
install.packages("devtools")
devtools::install_github("danielsuen/mixturebpe") 
```

## Example Code

In the package, we provide working implementations for obtaining point estimates both sequentially and in parallel.  In practice, we recommend running the code in parallel because there are multiple random initializations and many bootstrap samples.  The following example shows how to run code in parallel.

We first load the package and set up R for parallelization.  We identify the number of available cores on the machine via **detectCores()**.  We default to using 2 fewer cores than available to not take up all of the computer's resources.  The current progress while model fitting and bootstrapping will displayed using a progress bar.

```R
library(doParallel)
library(doSNOW)
library(mvtnorm)
library(mixturebpe)

#### init parallel cores

numCores = detectCores()
cluster = makeCluster(numCores-2)
registerDoSNOW(cluster)
```

Then, we simulate some data.  The data is simulated according to a mixture of binomial product experts model with $d=4$ outcome variables, $p=2$ covariates, and $K=3$ mixtures.  The parameter $\eta$ controls the amount of missingness, which is generated to be missing at random (MAR).

```R
#### generate simulated data

n = 2000
eta = 2.5
K = 3
p = 2
d = 4

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

Finally, we run the bootstrap (Algorithm 4):

```R
##### bootstrap

boot = bootstrap_parallel(X, Y_obs, theta_mar,
                          beta_mar, M_mar=100, s_max, numBootstraps=1000,
                          stop_eps_mar=1e-4, stop_eps_latent=1e-4,
                          num_imps=20)

#### unpack bootstrap results

beta_boot = boot$beta_boot
theta_boot = boot$theta_boot
```

We can construct $95\\%$ confidence intervals for a given parameter (say $\theta_{1,1}$) via the following code:

```R
upper = theta_mar[1,1] + 1.96 * sd(theta_boot[1,1,])
lower = theta_mar[1,1] - 1.96 * sd(theta_boot[1,1,])
```

Thus, the $95\\%$ confidence interval for $\theta_{1,1}$ is given by $(\text{lower}, \text{upper})$.
