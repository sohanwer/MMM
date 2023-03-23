library(rstan)
library(tidyverse)
library(loo)
library(ggplot2)
options(mc.cores = parallel::detectCores())

########################
model <- "
functions {
  vector adstock(vector x, int n_timesteps, vector r) {
  
    vector[n_timesteps] out;
    out = rep_vector(0, n_timesteps); // Initialize out to 0
    
    for (i in 1:n_timesteps) {
      if (i == 1) {
        out[i] = x[i] * r[i];
      } else {
        vector[i] tmp;
        for (j in 1:i) {
          tmp[j] = r[i-j+1];
        }
        out[i] = dot_product(x[1:i], tmp);
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
  vector<lower=0>[n_channels] mus;
  vector<lower=0>[n_channels] phis;
  vector<lower=0>[n_channels] betas;
  vector<lower=0>[n_channels] kappas; // half-saturation point. 
}

transformed parameters {
  vector[n_timesteps] predicted = rep_vector(0, n_timesteps);
  matrix[n_timesteps, n_channels] adstock_rates = rep_matrix(0,n_timesteps,n_channels);
  matrix[n_timesteps, n_channels] rois;
  matrix[n_timesteps, n_channels] adstocked_spend;

  
  for(i in 1:n_channels){
    for(j in 1:n_timesteps){
      adstock_rates[j,i] = inv_logit(neg_binomial_2_lpmf(j | mus[i], phis[i]));
    }
    adstocked_spend[, i] = adstock(spend[, i], n_timesteps, adstock_rates[,i]);

  }
  
  predicted = rep_vector(intercept, n_timesteps);
  
  for(i in 1:n_channels) {

    for (j in 1:n_timesteps) {
      rois[j, i] = betas[i] * kappas[i] / (adstocked_spend[j, i] + kappas[i]);
    }
    predicted += rois[, i] .* adstocked_spend[,i];
  }
  
}


model {
  
  
  intercept ~ normal(0,1);
  
  betas ~ normal(0,1);
  
  kappas ~ normal(0,1);
  
  sigma ~ exponential(1);
  
  for (i in 1:n_channels){
    mus[i] ~ cauchy(0,2);
    phis[i] ~ cauchy(0,2);

  };
  
  depvar ~ normal(predicted, sigma);
  
  
} 
generated quantities{
  vector[n_timesteps] log_lik;
  for ( i in 1:n_timesteps ) {
    log_lik[i] = normal_lpdf( depvar[i] | predicted[i], sigma );
  }
}"


adstock_nbinom <- function(x, r) {
  out <- numeric(length(x))
  
  for (i in 1:length(x)) {
    if (i == 1) {
      out[i] <- x[i] * r[i]
      
    } else{
      
      # Calculate dot product for i > 1
      out[i] <- sum(x[1:i] * rev(r[1:i]))
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
  dnbinom(1:(time), mu = mean_shift[i], 
          size = concentration[i])
}) 

rates = as.matrix(sample_df)


as_spend <- matrix(nrow=time, ncol=2)
rois <- as_spend

for (i in 1:ncol(spend)) {
  as_spend[,i] <- adstock_nbinom(spend[,i], rates[,i])
  rois[,i] <- true_betas[i] * true_kappas[i] / (as_spend[,i] + true_kappas[i])
}

dv <- rep(true_intercept, time)

for (i in 1:ncol(spend)) {
  dv <- dv + as_spend[,i] * rois[,i]
}

dv <- dv + rnorm(time, 0, 0.01)

stan_data <- list(n_timesteps = time, n_channels = 2, spend = spend, depvar = dv)


###################################################################
# Fit model
m <- rstan::stan_model(model_code = model)
n_samples <- 2000
fit <- sampling(m, data=stan_data, chains = 2, cores = 2, iter=n_samples)

####################################################################
# Inspect estimates
print(fit, pars=c("intercept", "betas", "kappas", "adstock_rates",'mus','phis'))


#####################################################################

## Diagnostics

## Check if chains have converged

stan_trace(fit)

####################################################################

####################################################################

# Visualise predicted parameter
post_df <- as.data.frame(fit)
extracted_fit <- rstan::extract(fit)


predicted_mean <- matrix(rep(NA, n_samples* time), nrow= n_samples)
predicted_observation <- matrix(rep(NA, n_samples* time), nrow= n_samples)

for(i in 1:time) {
  
  predicted_mean[,i] <- extracted_fit$predicted[,i]
  predicted_observation[,i] <- rnorm(n_samples, extracted_fit$predicted[,i], extracted_fit$sigma)
}

interval_mean <- HDInterval::hdi(predicted_mean, 0.9)
interval_prediction <- HDInterval::hdi(predicted_observation, 0.9)



model_data_df <- tibble(time = 1:time,
                        spend1 = spend[,1],
                        spend2 = spend[,2],
                        dv = dv) # dependent variable

model_data_df$low_interval_mean <- interval_mean[1,]
model_data_df$high_interval_mean <- interval_mean[2,]
model_data_df$low_interval_pred <- interval_prediction[1,]
model_data_df$high_interval_pred <- interval_prediction[2,]


p <- ggplot()
p2 <- p +
  
  geom_ribbon(data = model_data_df,
              aes(x = time, ymin = low_interval_mean, ymax = high_interval_mean),
              alpha = .5) +
  geom_ribbon(data = model_data_df,
              aes(x = time, ymin = low_interval_pred, ymax = high_interval_pred),
              alpha = .1) +
  
  geom_line(data = model_data_df,aes(x = time, y = dv), color = "blue") +
  
  
  labs(subtitle="HPDI Interval = 0.90")
p2


###################################################

## Task 2 Graph

###################################################
###################################################

# Extract posterior samples
    # Plot for media channel 1

## Get actual sales for channel 1
sales_ch1 <- rep(true_intercept, time)

sales_ch1 <- sales_ch1 + as_spend[,1] * rois[,1]
sales_ch1 <- sales_ch1 + rnorm(500,0,0.001)


extracted_fit <- rstan::extract(fit)

mu1_estimate <- mean(extracted_fit$mus[,1])
phi1_estimate <- mean(extracted_fit$phis[,1])
kappa1_estimate <- mean(extracted_fit$kappas[,1])
beta1_estimate <- mean(extracted_fit$betas[,1])
intercept_esimate = mean(extracted_fit$intercept)


estimated_adstocked_rates_1 = map_dfc (1:1,function(i) {
  dnbinom(1:(time), mu = mu1_estimate, 
          size = phi1_estimate)
}) 
estimated_adstocked_rates_1 = as.matrix(estimated_adstocked_rates_1)


as_spend_predicted_1 <- matrix(nrow=time, ncol=1)
rois_predicted_1 <- as_spend_predicted_1

as_spend_predicted_1 = adstock_nbinom(spend[,1], estimated_adstocked_rates_1)
rois_predicted_1 <- beta1_estimate* kappa1_estimate/ (as_spend_predicted_1 + 
                                                        kappa1_estimate)


dv_ch1_predicted <- rep(intercept_esimate, time)

dv_ch1_predicted <- dv_ch1_predicted + as_spend_predicted_1 * rois_predicted_1


dv_ch1_predicted <- dv_ch1_predicted + rnorm(time, 0, 0.01)

ch1_df = data.frame(spend[,1],sales_ch1,dv_ch1_predicted)
names(ch1_df) = c("Spend","Actual_Sales","Predicted_Sales")

ggplot(data = ch1_df, aes(x = Spend)) +
  geom_line(aes(y = Actual_Sales, color = "Actual Sales"), size = 0.5) +
  geom_smooth(aes(y = Actual_Sales, color = "Actual Sales"), method = "gam", se = FALSE) +
  geom_line(aes(y = Predicted_Sales, color = "Predicted Sales"), size = 0.5) +
  geom_smooth(aes(y = Predicted_Sales, color = "Predicted Sales"), method = "gam", se = FALSE) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Actual Sales", "Predicted Sales")) +
  xlab("Spend") +
  ylab("Sales") +
  ggtitle("Actual Sales vs. Predicted Sales - Channel 1")  


# Extract posterior samples
# Plot for media channel 1

## Get actual sales for channel 1
sales_ch2 <- rep(true_intercept, time)

sales_ch2 <- sales_ch2 + as_spend[,2] * rois[,2]
sales_ch2 <- sales_ch2 + rnorm(500,0,0.001)


mu2_estimate <- mean(extracted_fit$mus[,2])
phi2_estimate <- mean(extracted_fit$phis[,2])
kappa2_estimate <- mean(extracted_fit$kappas[,2])
beta2_estimate <- mean(extracted_fit$betas[,2])
intercept_esimate = mean(extracted_fit$intercept)


estimated_adstocked_rates_2 = map_dfc (1:1,function(i) {
  dnbinom(1:(time), mu = mu2_estimate, 
          size = phi2_estimate)
}) 
estimated_adstocked_rates_2 = as.matrix(estimated_adstocked_rates_2)


as_spend_predicted_2 <- matrix(nrow=time, ncol=1)
rois_predicted_2 <- as_spend_predicted_2

as_spend_predicted_2 = adstock_nbinom(spend[,2], estimated_adstocked_rates_2)
rois_predicted_2 <- beta2_estimate* kappa2_estimate/ (as_spend_predicted_2 + 
                                                        kappa2_estimate)


dv_ch2_predicted <- rep(intercept_esimate, time)

dv_ch2_predicted <- dv_ch2_predicted + as_spend_predicted_2 * rois_predicted_2


dv_ch2_predicted <- dv_ch2_predicted + rnorm(time, 0, 0.01)

ch2_df = data.frame(spend[,2],sales_ch2,dv_ch2_predicted)
names(ch2_df) = c("Spend","Actual_Sales","Predicted_Sales")

ggplot(data = ch2_df, aes(x = Spend)) +
  geom_line(aes(y = Actual_Sales, color = "Actual Sales"), size = 0.5) +
  geom_smooth(aes(y = Actual_Sales, color = "Actual Sales"), method = "gam", se = FALSE) +
  geom_line(aes(y = Predicted_Sales, color = "Predicted Sales"), size = 0.5) +
  geom_smooth(aes(y = Predicted_Sales, color = "Predicted Sales"), method = "gam", se = FALSE) +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Actual Sales", "Predicted Sales")) +
  xlab("Spend") +
  ylab("Sales") +
  ggtitle("Actual Sales vs. Predicted Sales - Channel 2")  

