library(rstan)
library(tidyverse)
library(loo)
options(mc.cores = parallel::detectCores())

########################



adstock_nbinom2 <- function(x, r, max_lag) {
  out <- numeric(length(x))
  
  for (i in 1:length(x)) {
    if (i == 1) {
      out[i] <- x[i] * r[i]
    } else {
      if (i <= max_lag) {
        out[i] <- sum(x[1:i] * rev(r[1:i]))
      } else {
        out[i] <- sum(x[(i - max_lag + 1):i] * rev(r[1:max_lag]))
      }
    }
  }
  
  return(out)
}


time <- 500

spend <- matrix(ncol=2, nrow= time)
spend[,1] <- abs(cumsum(rnorm(time,)))

spend[,2] <- abs(cumsum(rnorm(time,)))

true_betas <- c(0.5, 1.5)
mean_shift = runif(2, 0.5, 10)
concentration = runif(2, 0.1, 2)
true_kappas <- c(0.4, 0.99)
true_intercept <- 3

sample_df <- map_dfc (1:2, function(i) {
  dnbinom(1:(15), mu = mean_shift[i], 
          size = concentration[i],)
}) 

rates = as.matrix(sample_df)


as_spend <- matrix(nrow=time, ncol=2)
rois <- as_spend

for (i in 1:ncol(spend)) {
  as_spend[,i] <- adstock_nbinom2(spend[,i], rates[,i],max_lag=15)
  rois[,i] <- true_betas[i] * true_kappas[i] / (as_spend[,i] + true_kappas[i])
}

dv <- rep(true_intercept, time)

for (i in 1:ncol(spend)) {
  dv <- dv + as_spend[,i] * rois[,i]
}

dv <- dv + rnorm(time, 0, 0.01)

stan_data <- list(n_timesteps = time, n_channels = 2, spend = spend, depvar = dv)


