

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
    mu ~ normal(, ) T[0.001,];
    sigma ~ normal(, ) T[0.001,];
    lambda ~ gamma(5,0.1);
    
    
    target += k*log(lambda) - lambda*k;
    
    for (i in 1:k) {

    // Using truncated normal instead of normal
    real norm_const = 1 - normal_cdf(0 | mu, sigma);  // Normalisation constant
    real trunc_cdf = (normal_cdf(x[i] | mu, sigma) - normal_cdf(0 | mu, sigma)) / norm_const;
    real trunc_lpdf = normal_lpdf(x[i] | mu, sigma) - log(norm_const);
    target += lambda * trunc_cdf + trunc_lpdf;

    }


  }
