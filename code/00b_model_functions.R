library(cmdstanr)
library(posterior)
library(arrow)

# This is when we have multiple maxima per sample
# for this we need to data in slightly different format
#
stan_code_efsmult <- function(
  gamma_shape,
  gamma_rate,
  mu_prior_mean,
  mu_prior_sd
) {
  paste0(
    "
    data {
      int<lower=1> k;                    // Number of competitions/samples
      int<lower=1> n_obs;                // Total number of observations
      array[n_obs] real x;                   // All recorded fish sizes
      array[k] int n_per_sample; // Number of fish per sample
      array[k] int start_idx;   // Starting index for each sample
    }
    
    parameters {
      real<lower=0> mu;                  // Mean of the fish size distribution
      real<lower=0.001> sigma;           // Standard deviation
      real<lower=1> lambda;              // Estimated mean number of fish caught per competition
    }
    
    model {
      // Priors
      mu ~ normal(",
    mu_prior_mean,
    ", ",
    mu_prior_sd,
    ") T[0.001,];
       sigma ~ normal(",
    mu_prior_mean / 3,
    ", ",
    mu_prior_sd / 3,
    ") T[0.001,];
      lambda ~ gamma(",
    gamma_shape,
    ",",
    gamma_rate,
    ");
      
      target += -lambda * k;
      
      {
        
      // area of normal CDF above zero (normalisation constant)
      real norm_const = 1 - normal_cdf(0 | mu, sigma);  
        
      // looping over each sample
        for(j in 1:k) {
          
        target += n_per_sample[j] * log(lambda);

          // we use the smallest value since we know all other values are smaller than x_m
          // we could use the start_indx if we know the x values are in increasing order
          // but this is the safest way
          real smallest_in_m = min(x[start_idx[j]:(start_idx[j] + n_per_sample[j] - 1)]);

          // Truncated normal CDF
          real trunc_cdf = (normal_cdf(smallest_in_m | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;

          target += lambda * trunc_cdf;
          
          // looping over each maxima within the sample j
          for(i in 1:n_per_sample[j]) {
            int obs_idx = start_idx[j] + i - 1;
            real trunc_lpdf = normal_lpdf(x[obs_idx] | mu, sigma) - log(norm_const);
            target += trunc_lpdf;
          }
        }
      }
    }"
  )
}

stan_code_efs <- function(gamma_shape, gamma_rate, mu_prior_mean, mu_prior_sd) {
  paste0(
    "

  data {
    int<lower=1> k;  // Number of competitions
    vector[k] x;  // Recorded largest fish sizes at each competition
  }
  
  parameters {
    real<lower=0> mu;    // Mean of the fish size distribution
    real<lower=0.001> sigma; // Standard deviation
    real<lower=1> lambda; // Estimated mean number of fish caught per competition
  }
  
  model {
    // Priors
    mu ~ normal(",
    mu_prior_mean,
    ", ",
    mu_prior_sd,
    ") T[0.001,];
    sigma ~ normal(",
    mu_prior_mean / 3,
    ", ",
    mu_prior_sd / 3,
    ") T[0.001,];
    lambda ~ gamma(",
    gamma_shape,
    ",",
    gamma_rate,
    ");
    
    
    target += k*log(lambda) - lambda*k;
    
    for (i in 1:k) {

    // Using truncated normal instead of normal
    real norm_const = 1 - normal_cdf(0 | mu, sigma);  // Normalisation constant
    real trunc_cdf = (normal_cdf(x[i] | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;
    real trunc_lpdf = normal_lpdf(x[i] | mu, sigma) - log(norm_const);
    target += lambda * trunc_cdf + trunc_lpdf;

    }


  }"
  )
}

stan_code_evt <- {
  "data {
  int<lower=0> k;           // number of observations
  vector[k] x;              // observed maxima values
}

parameters {
  real loc;                  // location parameter
  real<lower=0> scale;      // scale parameter (must be positive)
  real shape;                  // shape parameter
}

model {
  // Priors
