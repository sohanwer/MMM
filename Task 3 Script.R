library(rstan)
library(tidyverse)
library(loo)
options(mc.cores = parallel::detectCores())

########################
multiplicative_model <- "
functions {
  vector adstock(vector x, int n_timesteps, real rate) {
  vector[n_timesteps] out = x;
    for (i in 1 : n_timesteps) {
     if (i == 1) {
       out[i] = out[i] * rate ;
     } else {
        // note rate is a 'decay' rather than carryover
       out[i] = out[i] * rate + (1-rate) * out[i-1];
     }
}
    return (out);
  }
}

data {
  int n_timesteps;
  int n_channels;
  matrix[n_timesteps, n_channels] spend;
  vector[n_timesteps] depvar;
}

parameters {
  real<lower=0> intercept;
  real<lower=0> sigma;
  vector<lower=0>[n_channels] betas;
  vector[n_channels] adstock_rates_raw;
  vector<lower=0>[n_channels] kappas; // half-saturation point. 
}

transformed parameters {
  vector[n_timesteps] log_predicted = rep_vector(0, n_timesteps);
  vector[n_channels] adstock_rates = inv_logit(adstock_rates_raw);
  matrix[n_timesteps, n_channels] rois;
  matrix[n_timesteps, n_channels] adstocked_spend;

  log_predicted = rep_vector(intercept, n_timesteps);
  for(i in 1:n_channels) {
    adstocked_spend[, i] = adstock(spend[, i], n_timesteps, adstock_rates[i]);
    for (j in 1:n_timesteps) {
      rois[j, i] = betas[i] * kappas[i] / (adstocked_spend[j, i] + kappas[i]);
    }
    log_predicted += log(rois[, i] .* adstocked_spend[,i]);
  }
}

model {
  intercept ~ normal(0,1);
  betas ~ normal(0,1);
  adstock_rates_raw ~ normal(0,1);
  kappas ~ normal(0,1);
  sigma ~ exponential(1);
  depvar ~ normal(log_predicted, sigma);
} 
generated quantities{
  vector[n_timesteps] log_lik;
  for ( i in 1:n_timesteps ) {
    log_lik[i] = normal_lpdf( depvar[i] | log_predicted[i], sigma );
  }
}"

additive_model <- "
functions {
  vector adstock(vector x, int n_timesteps, real rate) {
  vector[n_timesteps] out = x;
    for (i in 1 : n_timesteps) {
     if (i == 1) {
       out[i] = out[i] * rate ;
     } else {
        // note rate is a 'decay' rather than carryover
       out[i] = out[i] * rate + (1-rate) * out[i-1];
     }
}
    return (out);
  }
}

data {
  int n_timesteps;
  int n_channels;
  matrix[n_timesteps, n_channels] spend;
  vector[n_timesteps] depvar;
}

parameters {
  real<lower=0> intercept;
  real<lower=0> sigma;
  vector<lower=0>[n_channels] betas;
  vector[n_channels] adstock_rates_raw;
  vector<lower=0>[n_channels] kappas; // half-saturation point. 
}

transformed parameters {
  vector[n_timesteps] predicted = rep_vector(0, n_timesteps);
  vector[n_channels] adstock_rates = inv_logit(adstock_rates_raw);
  matrix[n_timesteps, n_channels] rois;
  matrix[n_timesteps, n_channels] adstocked_spend;

  predicted = rep_vector(intercept, n_timesteps);
  for(i in 1:n_channels) {
    adstocked_spend[, i] = adstock(spend[, i], n_timesteps, adstock_rates[i]);
    for (j in 1:n_timesteps) {
      rois[j, i] = betas[i] * kappas[i] / (adstocked_spend[j, i] + kappas[i]);
    }
    predicted += rois[, i] .* adstocked_spend[,i];
  }
}

model {
  intercept ~ normal(0,1);
  betas ~ normal(0,1);
  adstock_rates_raw ~ normal(0,1);
  kappas ~ normal(0,1);
  sigma ~ exponential(1);
  depvar ~ normal(predicted, sigma);
} 
generated quantities{
  vector[n_timesteps] log_lik;
  for ( i in 1:n_timesteps ) {
    log_lik[i] = normal_lpdf( depvar[i] | predicted[i], sigma );
  }
}"


adstock <- function(x, rate) {
  for (i in 1:length(x)) {
    if (i == 1) {
      x[i] <- x[i] * rate
    } else {
      x[i] <- x[i] * rate + (1-rate) * x[i-1];
    } }
  return (x) }

time <- 500

spend <- matrix(ncol=2, nrow= time)
spend[,1] <- abs(cumsum(rnorm(time,)))

spend[,2] <- abs(cumsum(rnorm(time,)))

true_betas <- c(0.5, 1.5)
true_adstocks <- c(0.1, 0.35)
true_kappas <- c(0.4, 0.99)
true_intercept <- 3

as_spend <- matrix(nrow=time, ncol=2)
rois <- as_spend

for (i in 1:ncol(spend)) {
  as_spend[,i] <- adstock(spend[,i], true_adstocks[i])
  rois[,i] <- true_betas[i] * true_kappas[i] / (as_spend[,i] + true_kappas[i])
}
dv <- rep(true_intercept, time)

for (i in 1:ncol(spend)) {
  dv <- dv + as_spend[,i] * rois[,i]
}

dv <- dv + rnorm(time, 0, 0.01)


stan_data <- list(n_timesteps = time, n_channels = 2, spend = spend, depvar = dv)




###################################################################
# Fit multiplicative model
m1 <- rstan::stan_model(model_code = multiplicative_model)
n_samples <- 2000
mult_fit <- sampling(m1, data=stan_data, chains = 2, cores = 2, iter=n_samples)
mult_extracted_fit <- rstan::extract(mult_fit)

## Fit additive model
m2 <- rstan::stan_model(model_code = additive_model)
add_fit <- sampling(m2, data=stan_data, chains = 2, cores = 2, iter=n_samples)
add_extracted_fit <- rstan::extract(add_fit)

#####################################################################
## Leave one-out CV


loo_mult <- loo(mult_fit, save_psis = TRUE)

loo_add <- loo(add_fit,save_psis = TRUE)





